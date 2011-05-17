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

#include "h5_pool.hpp"
#include "h5_helper.h"

#include "hdf5_type.h"
#include "hdf5_hid_holder.hpp"

#include "bs_hdf5_storage.h"

using namespace std;
using namespace boost::python;


namespace blue_sky
{
  error_h5_no_array::error_h5_no_array (std::string const &array_name)
  : bs_exception ("error_h5_no_array", array_name)
  {
  }

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
      hid_t dset = detail::open_dataset (g_id, name);
      hid_t dspace = H5Dget_space (dset);
      if (dspace < 0)
        {
          bs_throw_exception (boost::format ("Can't get dataspace for dataset %s in group %d") % name % g_id);
        }

      hid_t dtype = H5Dget_type (dset);
      if (dtype < 0)
        {
          bs_throw_exception (boost::format ("Can't get datatype for dataset %s in group %d") % name % g_id);
        }

      hsize_t dims_[10] = {0};
      int n_dims = H5Sget_simple_extent_dims (dspace, dims_, NULL);
      ((h5_pool *)m)->add_node (name, dset, dspace, dtype, n_dims, dims_);
      return 0; 
    }

  template <class T> h5_pool::map_t::iterator 
  h5_pool::add_node (const std::string &name, const hid_t dset, const hid_t dspace, 
                     const hid_t dtype, const int n_dims, const T *dims, const bool var_dims)
    {
      pair_t p;
      //printf ("Add node");

      p.first = name;
      p.second.dset = dset;
      p.second.dspace = dspace;
      p.second.dtype = dtype;
      p.second.plist = 0;
      p.second.n_dims = n_dims;
      p.second.var_dims = var_dims;
      p.second.size = 0;
      if (var_dims)
        {
          for (int i = 0; i < 6; ++i)
            {
              p.second.src_dims[i] = (hsize_t)dims[i];
            }
        }
      else
        {
          for (int i = 0; i < n_dims; ++i)
            {
              p.second.py_dims[i] = (npy_intp)dims[i];
              p.second.h5_dims[i] = (hsize_t)dims[i];
            }
        }
      return h5_map.insert (p).first;
    }

  void 
  h5_pool::fill_map ()
  {
    BS_ASSERT (group_id);
    clear_map ();
    herr_t r = H5Literate (group_id, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, &(h5_pool::it_group), this); 
    if (r < 0)
      {
        bs_throw_exception (boost::format ("iterate fails"));
      }
  }

  void 
  h5_pool::close_node (h5_pair &p)
    {
      H5Dclose (p.dset);
      H5Sclose (p.dspace);
      H5Tclose (p.dtype); 
      // FIXME: close plist
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
      if (!file_exists (fname_) || !detail::is_hdf5_file (fname_))
        {
          file_id = H5Fcreate (fname_.c_str (), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
          if (file_id < 0)
            {
              bs_throw_exception (boost::format ("Can't create h5 file: %s") % fname);
            }
        }
      else
        {
          file_id = H5Fopen (fname_.c_str (), H5F_ACC_RDWR, H5P_DEFAULT);
          if (file_id < 0)
            {
              bs_throw_exception (boost::format ("Can't open h5 file: %s") % fname);
            }
        }

      fname = fname_;

      if (detail::is_object_exists (file_id, path_))
        {
          group_id = H5Gopen (file_id, path_.c_str ());
          if (group_id < 0)
            {
              bs_throw_exception (boost::format ("Can't open existing group: %s") % path_);
            }
        }
      else
        {
          group_id = H5Gcreate (file_id, path_.c_str (), -1);
          if (group_id < 0)
            {
              bs_throw_exception (boost::format ("Can't create group: %s") % path_);
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

  t_long 
  h5_pool::calc_data_dims (map_t::iterator it)
  {
    t_long n, d;
    t_int old_n_dims;
    int i, flg;
    
    n = 1;
     
    if (it->second.var_dims)
      {
        flg = 0;
        old_n_dims = it->second.n_dims;
        it->second.n_dims = 0;
        for (i = 0; i < 3; i++)
          {
            d = it->second.src_dims[2 * i] * pool_dims[i] + it->second.src_dims[2 * i + 1];
            if (d)
              {
                if (d != it->second.h5_dims[it->second.n_dims])
                  flg = 1;
                it->second.h5_dims[it->second.n_dims] = d;
                it->second.py_dims[it->second.n_dims] = d;
                it->second.n_dims++;
                n *= d;
              }
          }
        if (it->second.dset && (flg || (old_n_dims != it->second.n_dims)))
          {
            H5Ldelete (group_id, it->first.c_str (), H5P_DEFAULT);
            close_node (it->second);
          }
       }
     else
       {
         for (i = 0; i < it->second.n_dims; i++)
           n *= it->second.py_dims[i];
       }

     it->second.size = n;
     return n;   
  }
  
  t_long 
  h5_pool::calc_data_dims (const std::string &name)
    {
      map_t::iterator it;
      it = h5_map.find (name);
      return calc_data_dims(it);
  }

  spv_float
  h5_pool::get_fp_data (std::string const &name)
  {
    spv_float a = get_fp_data_unsafe (name);
    if (!a || a->empty ())
      throw error_h5_no_array (name);

    return a;
  }

  spv_int
  h5_pool::get_i_data (std::string const &name)
  {
    spv_int a = get_i_data_unsafe (name);
    if (!a || a->empty ())
      throw error_h5_no_array (name);

    return a;
  }

  bool
  h5_pool::is_opened (std::string const &name)
  {
    BS_ASSERT (group_id >= 0);
    map_t::iterator it = h5_map.find (name);
    if (it == h5_map.end ())
      {
        bs_throw_exception (boost::format ("No array %s in pool") % name);
      }

    return it->second.dspace != 0;
  }

  template <typename T>
  BS_SP (T)
  get_data (std::string const &name, h5_pair const &p)
  {
    t_long n;
    BS_SP (T) a = BS_KERNEL.create_object (T::bs_type ());
    typename T::value_type def_value;
    
    if (p.dset)
      {
        n = H5Sget_simple_extent_npoints (p.dspace);
      }
    else
      {
        n = 1;
        for (int i = 0; i < p.n_dims; ++i)
          n *= p.py_dims[i];
      }   
    a->resize (n);
    a->reshape (p.n_dims, p.py_dims);
    
    if (p.dset)
      {
        herr_t r = H5Dread (p.dset, get_hdf5_type <typename T::value_type> (), H5S_ALL, H5S_ALL, H5P_DEFAULT, a->data ());
        if (r < 0)
          {
            bs_throw_exception (boost::format ("Can't read data: %s") % name);
          }
      }
    else
      {
        H5Pget_fill_value (p.plist, p.dtype, &def_value);
        a->assign(def_value);
      }

    return a;
  }

  spv_float 
  h5_pool::get_fp_data_unsafe (const std::string &name)
  {
    // FIXME: was <=
    if (group_id < 0)
      {
        bs_throw_exception (boost::format ("Get data %s: group not opened") % name);
      }
    if (h5_map.find (name) == h5_map.end ())
      {
        return spv_float ();
      }

    return get_data <v_float> (name, open_data (name));
  }

  spv_int 
  h5_pool::get_i_data_unsafe (const std::string &name)
  {
    // FIXME: was <=
    if (group_id < 0)
      {
        bs_throw_exception (boost::format ("Get data %s: group not opened") % name);
      }
    if (h5_map.find (name) == h5_map.end ())
      {
        return spv_int ();
      }

    return get_data <v_int> (name, open_data (name));
  }
    
    
  void
  h5_pool::declare_data (std::string const &name, hid_t dtype, void *value, int n_dims, npy_intp *dims, int var_dims)
  {
    if (group_id <= 0)
      {
        bs_throw_exception (boost::format ("Declare: %s, group not opened.") % name);
      }

    map_t::iterator it = h5_map.find (name);
    if (it != h5_map.end ())
      {
        h5_map.erase (it);
      }

    it = add_node (name, 0, 0, 0, n_dims, dims, var_dims);
    // FIXME: should we do copy of dtype?
    it->second.dtype = dtype;

    if (!detail::is_object_exists (group_id, name.c_str ()))
      {
        hid_t plist = H5Pcreate (H5P_DATASET_CREATE);
        if (plist < 0)
          {
            bs_throw_exception (boost::format ("Can't create property for dataset %s in group %d") % name % group_id);
          }

        herr_t r = H5Pset_fill_value (plist, it->second.dtype, value);
        if (r < 0)
          {
            bs_throw_exception (boost::format ("Can't set fill for property, dataset %s in group %d") % name % group_id);
          }

        it->second.plist = plist;
      }
  }

  h5_pair
  h5_pool::open_data (std::string const &name, int h5_write)
  {
    map_t::iterator it = h5_map.find (name);
    if (it == h5_map.end ())
      {
        bs_throw_exception (boost::format ("Can't open data: %s") % name);
      }

    if (!it->second.dset)
      {
        if (!detail::is_object_exists (group_id, name.c_str ()))
          {
            const h5_pair &p = it->second;
            if (h5_write)
              {
                hid_t dspace = H5Screate_simple (p.n_dims, p.h5_dims, NULL);
                if (dspace < 0)
                  {
                    bs_throw_exception (boost::format ("Can't create simple dataspace for dataset %s in group %d") % name % group_id);
                  }

                hid_t dset = H5Dcreate (group_id, name.c_str (), p.dtype, dspace, p.plist);
                if (dset < 0)
                  {
                    bs_throw_exception (boost::format ("Can't create dataset %s in group %d") % name % group_id);
                  }

                it->second.dspace = dspace;
                it->second.dset = dset;
              }
            calc_data_dims (it);
          }
        else 
          {
            hid_t dset = detail::open_dataset (group_id, name.c_str ());
            hid_t dspace = H5Dget_space (dset);
            if (dspace < 0)
              {
                bs_throw_exception (boost::format ("Can't get dataspace for dataset %s in group %d") % name % group_id);
              }

            it->second.dspace = dspace;
            it->second.dset = dset;
            calc_data_dims (it);
          }
      }

    return it->second;
  }

  int 
  h5_pool::declare_i_data (const std::string &name, t_int def_value, int n_dims, npy_intp *dims, int var_dims)
  {
    declare_data (name, get_hdf5_type <t_int> (), &def_value, n_dims, dims, var_dims);
    return 0;
  }

  int 
  h5_pool::declare_fp_data (const std::string &name, t_float def_value, int n_dims, npy_intp *dims, int var_dims)
  {
    declare_data (name, get_hdf5_type <t_float> (), &def_value, n_dims, dims, var_dims);
    return 0;
  }

  void
  h5_pool::set_data (std::string const &name, hid_t dtype, void *data, t_long data_size, int n_dims, npy_intp *dims, void *def_value)
  {
    if (h5_map.find (name) == h5_map.end ())
      {
        declare_data (name, dtype, def_value, n_dims, dims);
      }

    const h5_pair &p = open_data (name, 1);
    
    if (p.size > data_size)
      {
        bs_throw_exception (boost::format ("Size mismatch for %s: %d > %d") % name % p.size % data_size);
      }

    herr_t r = H5Dwrite (p.dset, p.dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    if (r < 0)
      {
        bs_throw_exception (boost::format ("Can't write data: %s") % name);
      }
  }
    
  int 
  h5_pool::set_i_data (const std::string &name, spv_int data, t_int def_value)
  {
    set_data (name, 
              get_hdf5_type <t_int> (), 
              data->data (), 
              data->size (), 
              data->ndim (),
              data->dims (),
              &def_value);
    return 0;
  }

  int 
  h5_pool::set_fp_data (const std::string &name, spv_float data, t_float def_value)
  {
    set_data (name, 
              get_hdf5_type <t_float> (), 
              data->data (), 
              data->size (), 
              data->ndim (), 
              data->dims (),
              &def_value);
    return 0;
  }
    
  void 
  h5_pool::set_pool_dims (t_long *dims, int ndims)
  {
    map_t::iterator it, e;
    int i;
    e = h5_map.end ();
    
    for (i = 0; i < 3; i++)
      pool_dims[i] = 0;
      
    n_pool_dims = ndims;
    for (i = 0; i < ndims; i++)
      pool_dims[i] = dims[i];
    
    for (it = h5_map.begin (); it != e; ++it)
      {
        calc_data_dims (it);
      }
  }
  

  // FIXME: obsolete
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

  // FIXME: obsolete
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
    
  std::string
  h5_pool::get_data_type(const std::string &name) const
  {
    map_t::const_iterator i;
    H5T_class_t dt;
    
    // FIXME: was <=
    if (group_id < 0)
      {
        bs_throw_exception (boost::format ("Get data type %s: group not opened") % name);
      }
      
    i = h5_map.find (name);  
    if (i == h5_map.end ())
      {
        bs_throw_exception (boost::format ("Get data type %s: array not found") % name);
      }
    dt =  H5Tget_class (i->second.dtype);
    if (dt == H5T_INTEGER)
      return "int";
    else if (dt == H5T_FLOAT)
      return "float";
    else
      return "unknown";
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
      s << "H5_pool: " << h5_map.size ()  << "\t [" << pool_dims[0];
      for (int j = 1; j < n_pool_dims; ++j)
        s << ", " <<pool_dims[j];
      s <<  "]" << std::endl;
      s << "-----------------------------------------------------------------------------------------\n";
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
           
          if (i->second.dset)
            s << std::setw (10) << " Created ";
          else
            s << std::setw (10) << " Declared ";
          
          if (i->second.var_dims)
            {
              s << " Variable Dims ";
              s <<  "[ " << i->second.src_dims[0];
              for (int j = 1; j < 6; ++j)
                s << ", " <<i->second.src_dims[j];
              s << "]\t";
           }
          
          s << std::endl;
        }
      s << "-----------------------------------------------------------------------------------------\n";
      return s.str ();
    }
    
     
  
  boost::python::list 
  h5_pool::py_list_data() const
  {
    map_t::const_iterator i, e;
    boost::python::list items;
    
    e = h5_map.end ();
    for (i = h5_map.begin (); i != e; ++i)
      items.append(i->first);
    return items;
  }
  
  void 
  h5_pool::py_set_pool_dims (boost::python::list &dims)
  {
    t_long arr_dims[10];
    t_int i, n_dims = len(dims);
    
    if (n_dims > 10)
      bs_throw_exception (boost::format ("py_set_pool_dims: too big dims number %d!") % n_dims);
    
    for (int i = 0; i < n_dims; ++i)
      {
        arr_dims[i] = boost::python::extract<int>(dims[i]);
      }
    set_pool_dims (arr_dims, n_dims);
  }
  
  boost::python::list 
  h5_pool::py_get_pool_dims ()
  {
    boost::python::list dims;
    
    for (int i = 0; i < n_pool_dims; ++i)
      {
        dims.append(pool_dims[i]);
      }
    return dims;
  }
    
#endif //BSPY_EXPORTING_PLUGIN
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (h5_pool);
  BLUE_SKY_TYPE_STD_COPY (h5_pool);

  BLUE_SKY_TYPE_IMPL (h5_pool, h5_pool_iface, "h5_pool", "Array pool stored in hdf5 file", "Array pool stored in HDF5 file");
}  // blue_sky namespace

