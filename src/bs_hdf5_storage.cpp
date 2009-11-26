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

#include <boost/tokenizer.hpp>

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

  hdf5_group_v2
  hdf5_file::operator[] (const std::string &name)
  {
    return hdf5_group_v2 (*this, name);
  }

  namespace detail {

    enum hid_type {
      hdf5_property_type
      , hdf5_dataspace_type
      , hdf5_dataset_type
      , hdf5_group_type
    };

    template <hid_type type>
    struct hid_closer
    {
    };

    template <>
    struct hid_closer <hdf5_property_type>
    {
      static void
      close (hid_t v)
      {
        H5Pclose (v);
      }
    };

    template <>
    struct hid_closer <hdf5_dataspace_type>
    {
      static void
      close (hid_t v)
      {
        H5Sclose (v);
      }
    };

    template <>
    struct hid_closer <hdf5_dataset_type>
    {
      static void
      close (hid_t v)
      {
        H5Dclose (v);
      }
    };

    template <>
    struct hid_closer <hdf5_group_type>
    {
      static void
      close (hid_t v)
      {
        H5Gclose (v);
      }
    };

    template <hid_type hid_type_t>
    struct hid_holder
    {
      hid_holder (hid_t h)
      : h_ (h)
      {
      }

      ~hid_holder ()
      {
        if (h_ >= 0)
          hid_closer <hid_type_t>::close (h_);
      }

      hid_holder &
      operator= (hid_t h)
      {
        if (h_ >= 0)
          hid_closer <hid_type_t>::close (h_);

        h_ = h;
        return *this;
      }

      operator hid_t () 
      {
        return h_;
      }

      hid_t h_;
    };

    hid_holder <hdf5_property_type>
    create_dataset_property ()
    {
      hid_t plist = H5Pcreate (H5P_DATASET_CREATE);
      if (plist < 0)
        {
          bs_throw_exception ("Can't create property for create dataset");
        }

      return hid_holder <hdf5_property_type> (plist);
    }

    hid_holder <hdf5_property_type>
    raw_transfer_property ()
    {
      hid_t plist = H5Pcreate (H5P_DATASET_XFER);
      if (plist < 0)
        {
          bs_throw_exception ("Can't create property for raw data transfer");
        }

      return plist;
    }
  }

  typedef detail::hid_holder <detail::hdf5_group_type>      hid_group_t;
  typedef detail::hid_holder <detail::hdf5_dataspace_type>  hid_dspace_t;
  typedef detail::hid_holder <detail::hdf5_dataset_type>    hid_dset_t;
  typedef detail::hid_holder <detail::hdf5_property_type>   hid_property_t;

  typedef hdf5::private_::hdf5_buffer__ hdf5_buffer_t;

  hdf5_file::~hdf5_file ()
  {
  }

  struct hdf5_storage_v2::impl 
  {
    typedef std::map <std::string, hid_t>   file_id_map_t;

    impl ()
    : grub_call_stack_ (false)
    , error_context_ ("")
    {
    }

    ~impl ()
    {
      file_id_map_t::iterator it = file_map_.begin (), e = file_map_.end ();
      for (; it != e; ++it)
        {
          if (it->second >= 0)
            {
              H5Fclose (it->second);
            }
        }
    }

    void
    register_file (const std::string &name, hid_t id)
    {
      file_id_map_t::iterator it = file_map_.find (name);
      if (it == file_map_.end ())
        {
          file_map_.insert (std::make_pair (name, id));
        }
      else
        {
          if (it->second != -1)
            {
              BOSWARN (section::app_info, level::warning) 
                << boost::format ("Register already opened file %s, (old id: %d, new id: %d)")
                  % name % it->second % id
                << bs_end;
            }
          else
            {
              it->second = id;
            }
        }
    }

    hid_t
    get_file_id (const std::string &name) const
    {
      file_id_map_t::const_iterator it = file_map_.find (name);
      if (it == file_map_.end ())
        {
          return -1;
        }

      return it->second;
    }

    hid_t
    get_file_id (hdf5_file &file)
    {
      if (file.file_id_ < 0)
        {
          file.file_id_ = get_file_id (file.file_name_);
        }

      if (file.file_id_ < 0)
        {
          file.file_id_ = H5Fcreate (file.file_name_.c_str (), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
          if (file.file_id_ < 0)
            {
              bs_throw_exception (boost::format ("Can't open file %s") % file.file_name_);
            }

          hid_group_t group = H5Gcreate (file.file_id_, "/results", NULL);
          if (group < 0)
            {
              bs_throw_exception (boost::format ("Can't create group '/results' in file %s") % file.file_name_);
            }

          register_file (file.file_name_, file.file_id_);
        }

      return file.file_id_;
    }

    void
    create_group_hierarchy (hid_t file_id, const std::string &group_name)
    {
      typedef boost::tokenizer <boost::char_separator <char> > tokenizer_t;
      boost::char_separator <char> s ("/");
      tokenizer_t t (group_name, s);
      std::string gn = "";
      for (tokenizer_t::iterator it = t.begin (); it != t.end (); ++it)
        {
          gn += "/" + *it;
          set_error_context (gn);
          if (!detail::is_object_exists (file_id, gn.c_str ()))
            {
              detail::hid_holder <detail::hdf5_group_type> group = H5Gcreate (file_id, gn.c_str (), H5P_DEFAULT);
              if (group < 0)
                {
                  bs_throw_exception (boost::format ("Can't create group %s") % gn);
                }
            }
        }
    }

    template <typename data_t>
    void
    write_to_hdf5 (hdf5_file &file, 
                   const std::string &group_name, 
                   const std::string &dataset_name, 
                   data_t &buffer)
    {
      hid_t file_id               = get_file_id (file);

      const int rank              = 2;
      hsize_t chunk_dimens[rank]  = {buffer.size, 1};
      hsize_t fspace_dimens[rank] = {buffer.size, 0};
      hsize_t mspace_dimens[rank] = {buffer.size, 1};
      hsize_t max_dimens[rank]    = {buffer.size, H5S_UNLIMITED};

      if (set_error_context (group_name) && !detail::is_object_exists (file_id, group_name.c_str ()))
        {
          hid_group_t group = H5Gcreate (file_id, group_name.c_str (), H5P_DEFAULT);
          if (group < 0)
            {
              create_group_hierarchy (file_id, group_name);
            }
        }
      hid_group_t group = detail::open_group (file_id, group_name.c_str ());

      if (set_error_context (dataset_name) && !detail::is_object_exists (group, dataset_name.c_str ()))
        {
          hid_property_t plist = detail::create_dataset_property ();
          if (H5Pset_chunk (plist, rank, chunk_dimens) < 0)
            {
              bs_throw_exception (boost::format ("Can't set chunk for dataset %s in group %s") % dataset_name % group_name);
            }

          hid_dspace_t fspace = H5Screate_simple (rank, fspace_dimens, max_dimens);
          if (fspace < 0)
            {
              bs_throw_exception (boost::format ("Can't create simple dataspace for dataset %s in group %s") % dataset_name % group_name);
            }

          hid_dset_t ds_value = H5Dcreate (group, dataset_name.c_str (), buffer.type, fspace, plist);
          if (ds_value < 0)
            {
              bs_throw_exception (boost::format ("Can't create dataset %s in group %s") % dataset_name % group_name);
            }
        }

      hid_dset_t dataset = detail::open_dataset (group, dataset_name.c_str ());
      hid_dspace_t space = H5Dget_space (dataset);
      if (space < 0)
        {
          bs_throw_exception (boost::format ("Can't get space for dataset %s in group %s") % dataset_name % group_name);
        }

      hsize_t size_value[] = {0, 0};
      if (H5Sget_simple_extent_dims (space, size_value, NULL) < 0)
        {
          bs_throw_exception (boost::format ("Can't get dataspace dims for dataset %s in group %s") % dataset_name % group_name);
        }

      hsize_t size[] = {size_value[0], size_value[1] + 1};
      if (H5Dextend (dataset, size) < 0)
        {
          bs_throw_exception (boost::format ("Can't extend dataset %s in group %s") % dataset_name % group_name);
        }

      hid_dspace_t fspace = H5Dget_space (dataset);
      if (fspace < 0)
        {
          bs_throw_exception (boost::format ("Can't get dataspace for dataset %s in group %s") % dataset_name % group_name);
        }

      hsize_t offset[] = {0, size_value[1]};
      if (H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL, mspace_dimens, NULL) < 0)
        {
          bs_throw_exception (boost::format ("Can't select hyperslab for dataset %s in group %s") % dataset_name % group_name);
        }

      hid_dspace_t mspace = H5Screate_simple (rank, mspace_dimens, max_dimens);
      if (mspace < 0)
        {
          bs_throw_exception (boost::format ("Can't create simple dataspace for dataset %s in group %s") % dataset_name % group_name);
        }

      hid_t status = H5Dwrite (dataset, buffer.type, mspace, fspace, detail::raw_transfer_property (), buffer.data ());
      if (status < 0)
        {
          bs_throw_exception (boost::format ("Can't write data (dataset %s, group %s)") % dataset_name % group_name);
        }
    }

    static herr_t
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

    static herr_t 
    hdf5_error_handler (hid_t, void *p)
    {
      impl *impl_ = static_cast <impl *> (p);

      hid_t stack_id = H5Eget_current_stack ();
      if (stack_id < 0)
        {
          bs_throw_exception ("Can't get current stack");
        }

      std::string hdf5_stack_walk;
      hid_t status = H5Ewalk2 (stack_id, H5E_WALK_DOWNWARD, impl::hdf5_walk_handler, &hdf5_stack_walk);
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

      if (impl_->grub_call_stack_)
        {
          BOSOUT (section::save_data, level::critical)
            << boost::format ("HDF5 error occured (context: %s): \n%sCall stack which prevents to error: %s")
                % impl_->error_context_ % hdf5_stack_walk % kernel_tools::get_backtrace (128)
            << bs_end;
        }
      else
        {
          BOSOUT (section::save_data, level::critical)
            << boost::format ("HDF5 error occured (context: %s): \n%s")
                % impl_->error_context_ % hdf5_stack_walk 
            << bs_end;
        }

      return 0;
    }

    bool 
    set_error_context (const std::string &s)
    {
      error_context_ = s;
      return true;
    }

    bool            grub_call_stack_;
    std::string     error_context_;
    file_id_map_t   file_map_;
  };

  hdf5_storage_v2::hdf5_storage_v2 ()
  : impl_ (new impl ())
  {
    H5Eset_auto2 (H5E_DEFAULT, impl::hdf5_error_handler, impl_);
  }

  hdf5_storage_v2::~hdf5_storage_v2 ()
  {
    delete impl_;
  }

  hdf5_storage_v2 *
  hdf5_storage_v2::instance ()
  {
    static hdf5_storage_v2 instance_;
    return &instance_;
  }

  hdf5_group_v2 &
  hdf5_group_v2::write_buffer (const char *dataset, const hdf5_buffer_t &buffer)
  {
    hdf5_storage_v2::instance ()->impl_->write_to_hdf5 (file_, name_, dataset, buffer);
    return *this;
  }

  hdf5_group_v2 &
  hdf5_group_v2::write_pod (const char *dataset, const hdf5_pod_t &pod)
  {
    hdf5_storage_v2::instance ()->impl_->write_to_hdf5 (file_, name_, dataset, pod);
    return *this;
  }

  hdf5_group_v2 &
  hdf5_group_v2::write_struct (const char *dataset, const hdf5_struct_t &s)
  {
    hdf5_storage_v2::instance ()->impl_->write_to_hdf5 (file_, name_, dataset, s);
    return *this;
  }

//bs stuff
  BLUE_SKY_TYPE_STD_CREATE(bs_hdf5_storage)
  BLUE_SKY_TYPE_STD_COPY(bs_hdf5_storage)
  BLUE_SKY_TYPE_IMPL_SHORT(bs_hdf5_storage, objbase, "bs_hdf5_storage")

}
