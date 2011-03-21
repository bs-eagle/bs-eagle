/** 
 * @file h5_pool.cpp
 * @brief implementation of #h5_pool_iface
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-03-05
*/
#include <hdf5.h>
#include <stdio.h>
#include <iomanip>

#include "h5_pool.h"
#include "h5_helper.h"

using namespace std;
using namespace boost::python;


namespace blue_sky
{
  h5_pool::h5_pool (bs_type_ctor_param) 
        : h5_pool_iface ()
    {
      init_all ();
      mute_flag = 0;
      
    }
  h5_pool::h5_pool (const h5_pool & /*src*/) : bs_refcounter ()
    {
      init_all ();
      mute_flag = 0;
    }

  herr_t 
  h5_pool::it_group (hid_t g_id, const char *name, const H5L_info_t * /*info*/, void * m)
    {
      hsize_t dims_[10];
      int n_dims;
      hid_t dset;
      hid_t dspace;
      hid_t dtype;
      
      dset = H5Dopen (g_id, name);
      if (dset < 0)
        {
          // TODO: report error
          throw;
        }
      dspace = H5Dget_space (dset);
      dtype = H5Dget_type (dset);
      n_dims = H5Sget_simple_extent_dims (dspace, dims_, NULL);
      ((h5_pool *)m)->add_node (std::string (name), dset, dspace, dtype, n_dims, dims_);
      return 0; 
    }

  template <class T> h5_pool::map_t::iterator 
  h5_pool::add_node (const std::string &name, const hid_t dset, const hid_t dspace, 
                     const hid_t dtype, const int n_dims, const T *dims)
    {
      pair_t p;
      //printf ("Add node");

      p.first = name;
      p.second.dset = dset;
      p.second.dspace = dspace;
      p.second.dtype = dtype;
      p.second.n_dims = n_dims;
      for (int i = 0; i < n_dims; ++i)
        {
          p.second.py_dims[i] = (npy_intp)dims[i];
          p.second.h5_dims[i] = (hsize_t)dims[i];
        }
      return h5_map.insert (p).first;
    }

  void 
  h5_pool::fill_map ()
    {
      clear_map ();
      if (group_id <= 0)
        {
          // TODO: report error
          throw;
        }
      H5Literate (group_id, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, &(h5_pool::it_group), this); 

    }

  void 
  h5_pool::close_node (h5_pair &p)
    {
      H5Dclose (p.dset);
      H5Sclose (p.dspace);
      H5Tclose (p.dtype);
    }

  void 
  h5_pool::clear_map ()
    {
      map_t::iterator i, e;
      e = h5_map.end ();
      for (i = h5_map.begin (); i != e; ++i)
        {
          close_node (i->second);
        }
      h5_map.clear ();

    }
  void 
  h5_pool::open_file (const std::string &fname_, const std::string &path_)
    {
      if (file_id)
        {
          close_file ();
        }
      if (file_exists (fname_))
        {
          htri_t ret = H5Fis_hdf5 (fname_.c_str ());
          if (ret > 0)
            {
              file_id = H5Fopen (fname_.c_str (), H5F_ACC_RDWR, H5P_DEFAULT);
              if (file_id < 0)
                {
                  // TODO: report error
                  throw;
                }
            }
          else
            {
              file_id = H5Fcreate (fname_.c_str (), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
              if (file_id < 0)
                {
                  // TODO: report error
                  throw;
                }
            }
        }
      else 
        // file dose not exist already, create it
        {
          file_id = H5Fcreate (fname_.c_str (), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
          if (file_id < 0)
            {
              // TODO: report error
              throw;
            }
        }
      fname = fname_;

      mute ();
      // create or open group
      group_id = H5Gopen (file_id, path_.c_str ());
      un_mute ();
      if (group_id < 0)
        {
          group_id = H5Gcreate (file_id, path_.c_str (), -1);
          if (group_id < 0)
            {
              // TODO: report error
              throw;
            }
        }
      path = path_;
      fill_map ();
    }
  
  void 
  h5_pool::close_file ()
    {
      clear_map ();
      if (group_id > 0)
        {
          H5Gclose (group_id);
        }
      if (file_id > 0)
        {
          H5Fclose (file_id);
        }
      init_all ();
    }

  spv_float 
  h5_pool::get_fp_data (const std::string &name)
    {
      map_t::iterator it;
      t_long n;
      spv_float a;

      if (group_id <= 0)
        {
          // TODO: report error
          return spv_float ();
        }
      it = h5_map.find (name);
      if (it == h5_map.end ())
        {
          // TODO: report error
          return spv_float ();
        }
      n = H5Sget_simple_extent_npoints (it->second.dspace);

      a = BS_KERNEL.create_object (v_float::bs_type ());
      a->resize (n);
      a->reshape (it->second.n_dims, it->second.py_dims); 
      if (sizeof (t_float) == sizeof (float))
        {
          H5Dread (it->second.dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(*a)[0]);
        }
      else if (sizeof (t_float) == sizeof (double))
        {
          H5Dread (it->second.dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(*a)[0]);
        }
      else
        {
          //TODO: print error message
          return spv_float ();
        }
      return a;
    }

  spv_int 
  h5_pool::get_i_data (const std::string &name)
    {
      map_t::iterator it;
      t_long n;
      spv_int a;

      if (group_id <= 0)
        {
          // TODO: report error
          return spv_int ();
        }
      it = h5_map.find (name);
      if (it == h5_map.end ())
        {
          // TODO: report error
          return spv_int ();
        }
      n = H5Sget_simple_extent_npoints (it->second.dspace);

      a = BS_KERNEL.create_object (v_int::bs_type ());
      a->resize (n);
      a->reshape (it->second.n_dims, it->second.py_dims); 
      if (sizeof (t_int) == sizeof (int))
        {
          H5Dread (it->second.dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(*a)[0]);
        }
      else if (sizeof (t_int) == sizeof (long))
        {
          H5Dread (it->second.dset, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(*a)[0]);
        }
      else
        {
          //TODO: print error message
          return spv_int ();
        }
      return a;
    }


  int 
  h5_pool::set_i_data (const std::string &name, spv_int data)
    {
      t_int ndims;
      const npy_intp *dims;
      map_t::iterator it, e;

      if (group_id <= 0)
        {
          // TODO: report error
          throw;
        }

      ndims = data->ndim ();
      dims = data->dims ();

      it = h5_map.find (name);
      if (it != h5_map.end ())
        {
          if (it->second.n_dims == ndims)
            {
              int flg = 0;
              for (int i = 0; i < ndims; ++i)
                {
                  if (it->second.py_dims[i] != dims[i])
                    {
                      flg = 1;
                      break;
                    }
                }
              if (!flg)
                {
                  //printf ("array OK\n");
                  if (sizeof (t_int) == sizeof (int))
                    {
                      H5Dwrite (it->second.dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, 
                                H5P_DEFAULT, &(*data)[0]);
                    }
                  else if (sizeof (t_int) == sizeof (long))
                    {
                      H5Dwrite (it->second.dset, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, 
                                H5P_DEFAULT, &(*data)[0]);
                    }
                  else
                    {
                      //TODO: print error message
                      throw;
                    }
                  return 0;
                }
            }
          // erase element 
          H5Ldelete (group_id, it->first.c_str (), H5P_DEFAULT);
          close_node (it->second);
          h5_map.erase (it);
        }
      // create new
      {
        printf ("Create %s\n", name.c_str ());

        map_t::iterator new_it = add_node (name, 0, 0, 0, ndims, dims);
        
        new_it->second.dspace = H5Screate_simple (ndims, new_it->second.h5_dims, NULL);
        
        if (sizeof (t_int) == sizeof (int))
          {
            new_it->second.dtype = H5T_NATIVE_INT;
          }
        else if (sizeof (t_int) == sizeof (long))
          {
            new_it->second.dtype = H5T_NATIVE_LONG;
          }
        else
          {
            //TODO: print error message
            throw;
          }
        new_it->second.dset = H5Dcreate (group_id, name.c_str (), 
                                         new_it->second.dtype, 
                                         new_it->second.dspace,
                                         H5P_DEFAULT);
        H5Dwrite (new_it->second.dset, new_it->second.dtype, 
                  H5S_ALL, H5S_ALL, H5P_DEFAULT, &(*data)[0]);
      }
      return 0;
    }

  int 
  h5_pool::set_fp_data (const std::string &name, spv_float data)
    {
      t_int ndims;
      const npy_intp *dims;
      map_t::iterator it, e;

      if (group_id <= 0)
        {
          // TODO: report error
          throw;
        }

      ndims = data->ndim ();
      dims = data->dims ();

      it = h5_map.find (name);
      if (it != h5_map.end ())
        {
          if (it->second.n_dims == ndims)
            {
              int flg = 0;
              for (int i = 0; i < ndims; ++i)
                {
                  if (it->second.py_dims[i] != dims[i])
                    {
                      flg = 1;
                      break;
                    }
                }
              if (!flg)
                {
                  //printf ("array OK\n");
                  if (sizeof (t_float) == sizeof (float))
                    {
                      H5Dwrite (it->second.dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, 
                                H5P_DEFAULT, &(*data)[0]);
                    }
                  else if (sizeof (t_float) == sizeof (double))
                    {
                      H5Dwrite (it->second.dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
                                H5P_DEFAULT, &(*data)[0]);
                    }
                  else
                    {
                      //TODO: print error message
                      throw;
                    }
                  return 0;
                }
            }
          // erase element 
          H5Ldelete (group_id, it->first.c_str (), H5P_DEFAULT);
          close_node (it->second);
          h5_map.erase (it);
        }
      // create new
      {
        printf ("Create %s\n", name.c_str ());

        map_t::iterator new_it = add_node (name, 0, 0, 0, ndims, dims);

        new_it->second.dspace = H5Screate_simple (ndims, new_it->second.h5_dims, NULL);
        
        if (sizeof (t_float) == sizeof (float))
          {
            new_it->second.dtype = H5T_NATIVE_FLOAT;
          }
        else if (sizeof (t_float) == sizeof (double))
          {
            new_it->second.dtype = H5T_NATIVE_DOUBLE;
          }
        else
          {
            //TODO: print error message
            throw;
          }
        new_it->second.dset = H5Dcreate (group_id, name.c_str (), 
                                         new_it->second.dtype, 
                                         new_it->second.dspace,
                                         H5P_DEFAULT);
        H5Dwrite (new_it->second.dset, new_it->second.dtype, 
                  H5S_ALL, H5S_ALL, H5P_DEFAULT, &(*data)[0]);
      }
      return 0;
    }

  void
  h5_pool::mute ()
    {
      if (!mute_flag)
        {
          /* Save old error handler */
          H5Eget_auto2 (H5E_DEFAULT, &old_func, &old_client_data);

          /* Turn off error handling */
          H5Eset_auto2 (H5E_DEFAULT, NULL, NULL);

          mute_flag = 1;
        }
    }

  void
  h5_pool::un_mute ()
    {
      if (mute_flag)
        {
          /* Restore previous error handler */
          H5Eset_auto2 (H5E_DEFAULT, old_func, old_client_data);
          mute_flag = 0;
        }
    }

#ifdef BSPY_EXPORTING_PLUGIN
  std::string 
  h5_pool::py_str () const
    {
      std::stringstream s;
      std::string ss; 
      map_t::const_iterator i, e;
      H5T_class_t dt;

      e = h5_map.end ();
      s << "H5_pool: " << h5_map.size () << std::endl;
      s << "-----------------------------------------------------------\n";
      for (i = h5_map.begin (); i != e; ++i)
        {
          std::stringstream sdim;
          s << std::setw (15) <<i->first << "\t [";
          sdim  << i->second.py_dims[0];
          for (int j = 1; j < i->second.n_dims; ++j)
            sdim << ", " <<i->second.py_dims[j];
          s << std::setw (15) << sdim.str () << "]\t";
          dt =  H5Tget_class (i->second.dtype);
          if (dt == H5T_INTEGER)
            ss = "Integer ";
          else if (dt == H5T_FLOAT)
            ss = "Float ";
          else
            ss = "Unknown";

          s << std::setw (8) << ss << H5Tget_precision (i->second.dtype);
          s << std::endl;
        }
      s << "-----------------------------------------------------------\n";
      return s.str ();
    }
#endif //BSPY_EXPORTING_PLUGIN
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (h5_pool);
  BLUE_SKY_TYPE_STD_COPY (h5_pool);

  BLUE_SKY_TYPE_IMPL (h5_pool, h5_pool_iface, "h5_pool", "Array pool stored in hdf5 file", "Array pool stored in HDF5 file");
}  // blue_sky namespace
