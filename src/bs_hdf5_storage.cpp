/**
 * @file bs_hdf5_storage.cpp
 * @brief plugin for using hdf5 library using
 * @author Sayfullin Ilshat
 * @date 2009-01-28
 */

#include "bs_hdf5_storage.h"
#include "bs_kernel.h"
#include "bs_kernel_tools.h"
#include <stdlib.h>
#include "py_bs_hdf5_storage.h"
#include "bos_report.h"

using namespace blue_sky;
using namespace blue_sky::python;
using namespace boost::python;

namespace blue_sky
  {

  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("bs_hdf5_storage", "1.0.0", "Blue Sky HDF5 storage plugin", "Blue Sky HDF5 storage plugin", "bs_hdf5_storage")

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    const plugin_descriptor & pd = *bs_init.pd_;
    bool res = BLUE_SKY_REGISTER_TYPE(pd, bs_hdf5_storage);

    ///res &= blue_sky::give_kernel::Instance().register_type(*bs_init.pd_, bs_hdf5_storage::bs_type());
    ///BS_ASSERT (res);

    return res;
  }

#ifdef BSPY_EXPORTING_PLUGIN
  BLUE_SKY_INIT_PY_FUN
  {
    py_export_bs_hdf5_storage ();
  }
#endif //BSPY_EXPORTING_PLUGIN

  herr_t
  hdf5_walk_handler (unsigned n, const H5E_error2_t *err, void *p)
  {
    //BOSOUT (section::save_data, level::critical) 
    //  << "walk: " << n
    //  << bs_end;
    

    size_t class_name_size = H5Eget_class_name (err->cls_id, 0, 0);
    size_t major_size      = H5Eget_msg (err->maj_num, 0, 0, 0);
    size_t minor_size      = H5Eget_msg (err->min_num, 0, 0, 0);

    if (class_name_size <= 0 || major_size <= 0 || minor_size <= 0)
      {
        bs_throw_exception (boost::format ("class_name_size: %d, major_size: %d, minor_size: %d") 
          % class_name_size % major_size % minor_size);
      }

    shared_vector <char> class_name (2 * class_name_size, 0);
    shared_vector <char> major_name (2 * major_size, 0);
    shared_vector <char> minor_name (2 * minor_size, 0);

    if (H5Eget_class_name (err->cls_id, &class_name[0], 2 * class_name_size - 1) <= 0)
      {
        bs_throw_exception ("Can't get class name");
      }
    if (H5Eget_msg (err->maj_num, NULL, &major_name[0], 2 * major_size - 1) <= 0)
      {
        bs_throw_exception (boost::format ("Can't get major message (class: %s)") % &class_name[0]);
      }
    if (H5Eget_msg (err->min_num, NULL, &minor_name[0], 2 * minor_size - 1) <= 0)
      {
        bs_throw_exception (boost::format ("Can't get minor message (class: %s, major: %s)") 
          % &class_name[0] % &major_name[0]);
      }

    std::string *str = static_cast <std::string *> (p);


    (*str) += (boost::format ("#%u: %s in %s(): line %u\n\t\tclass: %s\n\t\tmajor: %s\n\t\tminor: %s\n") 
          % n % err->file_name % err->func_name % err->line
          % &class_name[0]
          % &major_name[0]
          % &minor_name[0]).str ();

    return 0;
  }

  herr_t 
  hdf5_error_handler (hid_t, void *)
  {
    std::string hdf5_stack_walk;

    hid_t stack_id = H5Eget_current_stack ();
    if (stack_id < 0)
      {
        bs_throw_exception ("Can't get current stack");
      }

    hid_t status = H5Ewalk2 (stack_id, H5E_WALK_DOWNWARD, hdf5_walk_handler, &hdf5_stack_walk);
    if (status < 0)
      {
        H5E_type_t type;
        size_t s = H5Eget_msg (status, &type, NULL, 0);
        if (!s)
          {
            bs_throw_exception ("Can't get error message");
          }

        shared_vector <char> buffer (s + 1, 0);
        s = H5Eget_msg (status, &type, &buffer[0], s);

        bs_throw_exception (boost::format ("Can't set walk2 handler: %s") % &buffer[0]);
      }

    H5Eclose_stack (stack_id);


    BOSOUT (section::save_data, level::critical)
    << "Error in HDF5 occured: \n"
    << hdf5_stack_walk
    << "Call stack which prevent to error: "
    << kernel_tools::get_backtrace (128)
    << bs_end;


    return 0;
  }

  bs_hdf5_storage::bs_hdf5_storage (bs_type_ctor_param param /* = NULL */)
      : file_id (0)
      , open_flag (0)
      //, group_wells (0)
      //, group_fip (0)
  {
    H5Eset_auto2 (H5E_DEFAULT, hdf5_error_handler, 0);
  }

  bs_hdf5_storage::bs_hdf5_storage (const bs_hdf5_storage& src)
  : bs_refcounter (src)
  , objbase (src)
  {
    *this = src;
  }

//// Constructor
//  bs_hdf5_storage::bs_hdf5_storage()
//  {
//    file_id = 0;
//  }
//// Destuctor
//  bs_hdf5_storage::~bs_hdf5_storage()
//  {
//    this->close();
//  }

  int bs_hdf5_storage::open (const char* filename, int flag)
  {
    // Close previous opened file (if it is opened)
    close ();

    open_flag = H5F_ACC_RDWR;
    if (flag & ST_ACC_RO)
      open_flag = H5F_ACC_RDONLY;

    if (flag & ST_ACC_CREATE)
      {
        file_id = H5Fcreate (filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (file_id < 0)
          throw bs_exception ("bs_hdf5_storage", "file_create_error");
        H5Fclose (file_id);
      }

    file_id = H5Fopen(filename, open_flag, H5P_DEFAULT);
    if (file_id < 0)
      {
        throw bs_exception ("bs_hdf5_storage", "file_open_error");
      }

    if ((size_t)open_flag == (size_t)H5F_ACC_RDONLY)    //read
      {
        //group_wells = open_group (file_id, "wells");
        //group_fip = open_group (file_id, "fip");
      }
    else    //write
      {
        hid_t group = H5Gcreate (file_id, "/results", NULL);
        H5Gclose (group);
      }
    return 0;
  }

#ifdef _MPI
// create Hdf5 file for writing
  hid_t* bs_hdf5_storage::mpi_create_indiv_file(const char* filename)
  {
    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
    indiv_file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist);
    if (file_out < 0)
      {
        rep->print (LOG_INIT_SECTION, LOG_ERR,
                    GET_TEXT ("Error: can't create HDF5 file %s\n"), filename);
        return 0;
      }
    H5Pclose(plist);
    return &indiv_file;
  }
#endif

// Close Hdf5 file
  int bs_hdf5_storage::close ()
  {
    if ((size_t)open_flag == (size_t)H5F_ACC_RDONLY)
      {
        //if (group_wells > 0)
        //  H5Gclose(group_wells);
        //if (group_fip > 0)
        //  H5Gclose(group_fip);
      }
    else // ST_CREATE, ST_RDWR: output file
      {
#ifdef _MPI
        if (H5Fclose (indiv_file))
          return -2;
        indiv_file = 0;
#endif
      }

    /*
      fprintf (stderr, "at H5Fclose: there are %d objs_all still opened!\n",H5Fget_obj_count(file_out, H5F_OBJ_ALL));
      fprintf (stderr, "at H5Fclose: there are %d files still opened!\n",H5Fget_obj_count(file_out, H5F_OBJ_FILE));
      fprintf (stderr, "at H5Fclose: there are %d groups still opened!\n",H5Fget_obj_count(file_out, H5F_OBJ_GROUP));
      fprintf (stderr, "at H5Fclose: there are %d datatypes still opened!\n",H5Fget_obj_count(file_out, H5F_OBJ_DATATYPE));
      fprintf (stderr, "at H5Fclose: there are %d attrs still opened!\n",H5Fget_obj_count(file_out, H5F_OBJ_ATTR));
      fprintf (stderr, "at H5Fclose: there are %d datasets still opened!\n",H5Fget_obj_count(file_out, H5F_OBJ_DATASET));
    */

    if (file_id > 0)
      {
        if (H5Fclose (file_id))
          bs_throw_exception ("Error occured while file closing");
      }

    file_id = 0;
    return 0;
  }

  bool bs_hdf5_storage::is_opened () const
  {
    if (file_id)
    return true;
    return false;
  }

// Creates a new group for object
  hid_t bs_hdf5_storage::begin_object (const hid_t location, const std::string& object_name) const
    {
      hid_t new_group = -1;
      // Check group already exists
      htri_t status = H5Lexists (location, object_name.c_str(), NULL);
      if (status == 0)  // Create group.
      {
        new_group = H5Gcreate (location, object_name.c_str(), H5P_DEFAULT);
        }
      else // Open group.
        {
          new_group = H5Gopen (location, object_name.c_str());
        }

    if (new_group < 0)
    {
      bs_throw_exception ("obj_opening_error");
      }

    return new_group;
  }

  hid_t bs_hdf5_storage::begin_object (const std::string &name) const
  {
    return begin_object (file_id, name);
  }

// Get rank of data set. Returns number of dimension
  int bs_hdf5_storage::get_rank (const hid_t location, const std::string& name) const
  {
    hid_t dataset = H5Dopen (location, name.c_str());
    hid_t dataspace = H5Dget_space (dataset);
    int rank;
    if ( (rank = H5Sget_simple_extent_ndims (dataspace)) < 0 )
    {
      throw bs_exception("bs_hdf5_storage", "rank_error");
      }
    H5Dclose (dataset);
    H5Sclose (dataspace);
    return rank;
  }


// Get dimensions.
  int bs_hdf5_storage::get_dimensions (const hid_t location, const std::string& name, int *dimensions) const
  {
    hid_t dataset;
    if ( (dataset = H5Dopen (location, name.c_str())) < 0 )
    {
      throw bs_exception("bs_hdf5_storage", "dataset_error");
        return -1;
      }
    hid_t dataspace = H5Dget_space (dataset);
    int rank = H5Sget_simple_extent_ndims (dataspace);

    hsize_t *t_dimensions = new hsize_t[rank];

    H5Sget_simple_extent_dims (dataspace, t_dimensions, NULL);

    for (int i=0; i<rank; ++i)
    dimensions[i] = t_dimensions[i];

    delete[] t_dimensions;
    H5Dclose (dataset);
    H5Sclose (dataspace);
    return 0;
  }

  /**
   * \brief Create char-dataset in specified location in hdf5 file and write str to it.
   *        Skip, if already exists.
   * \return 0
   */
  int bs_hdf5_storage::write_string_to_hdf5 (const hid_t location, const char *dataset_name, const std::string &str) const
    {
      //check if not already exists
      htri_t status = H5Lexists (location, dataset_name, NULL);
      if (status == 0)
      {
        // write data for the current well
        hid_t dataspace_group = H5Screate (H5S_SCALAR);
          // create string type of given length
          hid_t datatype_group = H5Tcopy (H5T_C_S1);
          H5Tset_size (datatype_group, str.length ());
          hid_t dataset_group = H5Dcreate (location, dataset_name, datatype_group, dataspace_group, H5P_DEFAULT);
          if (dataset_group >= 0)//succesfully created
            {
              hid_t fspace = H5Dget_space (dataset_group);
              hid_t plist = H5Pcreate (H5P_DATASET_XFER);
              status = H5Dwrite (dataset_group, datatype_group, dataspace_group, fspace, plist, str.c_str ());
              H5Pclose (plist);
              H5Sclose (fspace);
              H5Dclose (dataset_group);
              if (status < 0)
                return -3;
            }
          H5Tclose (datatype_group); // close type
        }
    return 0;
  }

//bs stuff
  BLUE_SKY_TYPE_STD_CREATE(bs_hdf5_storage)
  BLUE_SKY_TYPE_STD_COPY(bs_hdf5_storage)
  BLUE_SKY_TYPE_IMPL_SHORT(bs_hdf5_storage, objbase, "bs_hdf5_storage")

}
