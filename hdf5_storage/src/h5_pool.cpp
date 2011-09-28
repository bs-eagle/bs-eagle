/** 
 * @file h5_pool.cpp
 * @brief implementation of #h5_pool_iface
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-03-05
*/
#ifdef BSPY_EXPORTING_PLUGIN
#include <boost/python.hpp>
#endif

#include <hdf5.h>
#include <stdio.h>
#include <iomanip>

#include "h5_pool.hpp"
#include "h5_helper.h"

#include "hdf5_type.h"
#include "hdf5_type_to_hid.hpp"
#include "hdf5_hid_holder.hpp"
#include "hdf5_functions.h"

#include "bos_report.h"

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

      hid_t plist = H5Dget_create_plist (dset);
      if (plist < 0)
        {
          bs_throw_exception (boost::format ("Can't get property list for dataset %s in group %d") % name % g_id);
        }

      hsize_t dims_[10] = {0};
      // FIXME: check 
      int n_dims = H5Sget_simple_extent_dims (dspace, dims_, NULL);
      map_t::iterator it = ((h5_pool *)m)->add_node (name, dset, dspace, dtype, plist, n_dims, dims_, 0);

      BOSOUT (section::h5, level::low) 
        << "Node from file: " << name << bs_end;
      return 0; 
    }

  template <class T> h5_pool::map_t::iterator 
  h5_pool::add_node (const std::string &name, const hid_t dset, const hid_t dspace, 
                     const hid_t dtype, const hid_t plist, 
                     const int n_dims, const T *dims, const bool var_dims)
    {
      BS_ASSERT (dtype >= 0) (name);
      BS_ASSERT (plist >= 0) (name);
      pair_t p;

      p.first = name;
      p.second.dset = dset;
      p.second.dspace = dspace;
      p.second.dtype = dtype;
      p.second.plist = plist;
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
          p.second.size = 1;
          // FIXME: check n_dims
          for (int i = 0; i < n_dims; ++i)
            {
              p.second.py_dims[i] = (npy_intp)dims[i];
              p.second.h5_dims[i] = (hsize_t)dims[i];

              p.second.size *= dims[i];
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
      if (p.dset >= 0)
        {
          H5Dclose (p.dset);
          p.dset = -1;
        }
      if (p.dspace >= 0)
        {
          H5Sclose (p.dspace);
          p.dspace = -1;
        }
      if (p.dtype >= 0)
        {
          H5Tclose (p.dtype); 
          p.dtype = -1;
        }
      if (p.plist >= 0)
        {
          H5Pclose (p.plist);
          p.plist = -1;
        }
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
      // FIXME:
      if (file_id)
        {
          close_file ();
        }

      file_id = H5Fcreate (fname_.c_str (), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      if (file_id < 0)
        {
          bs_throw_exception (boost::format ("Can't create h5 file: %s") % fname);
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

  typedef hsize_t pool_dims_t[3];
  // calcs dims only for var_dims
  t_long 
  calc_data_dims (std::string const &name, h5_pair &p, t_int n_dims, pool_dims_t const &pool_dims)
  {
    BOOST_STATIC_ASSERT (sizeof (t_long) >= sizeof (npy_intp));
    //BOOST_STATIC_ASSERT (sizeof (t_long) >= sizeof (hsize_t));

    p.n_dims = 0;
    p.size = 1;
    for (int i = 0; i < 3; i++)
      {
        t_long d = p.src_dims[2 * i] * pool_dims[i] + p.src_dims[2 * i + 1];
        if (d)
          {
            p.h5_dims[p.n_dims] = d;
            p.py_dims[p.n_dims] = d;
            p.n_dims++;
            p.size *= d;
          }
      }

    if (p.dset > 0)
      {
        t_long size = p.size;
        p.size = -1; // set wrong size to push dataset recreation when next set_data is called
        return size;
      }
      
    return p.size;   
  }

  typedef hsize_t src_dims_t[6];
  t_long
  calc_size (t_int n_pool_dims, pool_dims_t const &pool_dims, src_dims_t const &src_dims)
  {
    t_long size = 1;
    for (t_int i = 0; i < n_pool_dims; i++)
      {
        t_long d = src_dims[2 * i] * pool_dims[i] + src_dims[2 * i + 1];
        if (d)
          {
            size *= d;
          }
      }

    BS_ASSERT (size > 0) (size);
    return size;
  }
  
  t_long 
  h5_pool::calc_data_dims (const std::string &name)
  {
    map_t::iterator it = h5_map.find (name);
    if (it == h5_map.end ())
      bs_throw_exception (boost::format ("No array %s in pool") % name);

    h5_pair &p = it->second;
    if (p.var_dims)
      {
        blue_sky::calc_data_dims(name, p, n_pool_dims, pool_dims);
      }

    t_long size = calc_size (n_pool_dims, pool_dims, p.src_dims);
    if (p.size != size)
      {
        BOSOUT (section::h5, level::warning)
          << boost::format ("Size mismatch for %s in pool: %d == %d") % name % size % p.size
          << ". Calculated size will be used."
          << bs_end;
      }

    return size;
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
    BS_ASSERT (group_id >= 0) (name);
    map_t::iterator it = h5_map.find (name);
    if (it == h5_map.end ())
      {
        bs_throw_exception (boost::format ("No array %s in pool") % name);
      }

    return it->second.dspace >= 0;
  }

  template <typename T>
  BS_SP (T)
  get_data (std::string const &name, h5_pair const &p)
  {
    t_long n = H5Sget_simple_extent_npoints (p.dspace);
    BS_SP (T) a = BS_KERNEL.create_object (T::bs_type ());

    a->resize (n);
    a->reshape (p.n_dims, p.py_dims);

    herr_t r = H5Dread (p.dset, get_hdf5_type <hdf5_type_helper <typename T::value_type>::type> (), H5S_ALL, H5S_ALL, H5P_DEFAULT, a->data ());
    if (r < 0)
      {
        bs_throw_exception (boost::format ("Can't read data: %s") % name);
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
    BOOST_STATIC_ASSERT (sizeof (t_long) >= sizeof (npy_intp));
    //BOOST_STATIC_ASSERT (sizeof (t_long) >= sizeof (hsize_t));

    if (group_id <= 0)
      {
        bs_throw_exception (boost::format ("Declare: %s, group not opened.") % name);
      }

    map_t::iterator it = h5_map.find (name);
    if (it != h5_map.end ())
      {
        BOSOUT (section::h5, level::warning) << "Declared data " << name << " already exists in pool" << bs_end;
        h5_pair &p = it->second;

        for (int i = 0; i < 6; ++i)
          {
            p.src_dims[i] = dims[i];
          }
      }
    else
      {
        hid_t plist = H5Pcreate (H5P_DATASET_CREATE);
        if (plist < 0)
          {
            bs_throw_exception (boost::format ("Can't create property for dataset %s in group %d") % name % group_id);
          }

        herr_t r = H5Pset_fill_value (plist, dtype, value);
        if (r < 0)
          {
            bs_throw_exception (boost::format ("Can't set fill for property, dataset %s in group %d") % name % group_id);
          }

        hid_t dtype_copy = H5Tcopy (dtype);
        if (dtype_copy < 0)
          {
            bs_throw_exception (boost::format ("Can't copy datatype for %s in group %d") % name % group_id);
          }

        add_node (name, -1, -1, dtype_copy, plist, n_dims, dims, var_dims);
      }
  }

  h5_pair
  h5_pool::open_data (std::string const &name, const int)
  {
    map_t::iterator it = h5_map.find (name);
    if (it == h5_map.end ())
      {
        bs_throw_exception (boost::format ("Can't open data: %s") % name);
      }

    h5_pair &p = it->second;
    if (p.var_dims && p.size != calc_size (n_pool_dims, pool_dims, p.src_dims))
      {
        if (p.dset >= 0)
          {
            H5Ldelete (group_id, name.c_str (), H5P_DEFAULT);
          }

        hid_t dtype = p.dtype;
        hid_t plist = p.plist;

        p.dtype = p.plist = -1; // prevent closing of dtype and plist
        close_node (p);

        p.dtype = dtype;
        p.plist = plist;

        blue_sky::calc_data_dims (name, p, n_pool_dims, pool_dims);
      }

    if (p.dset < 0)
      {
        if (!detail::is_object_exists (group_id, name.c_str ()))
          {
            hid_t dspace = H5Screate_simple (p.n_dims, p.h5_dims, NULL);
            if (dspace < 0)
              {
                bs_throw_exception (boost::format ("Can't create simple dataspace for dataset %s in group %d") % name % group_id);
              }

            BS_ASSERT (p.dtype >= 0) (name);
            BS_ASSERT (p.plist) (name);

            hid_t dset = H5Dcreate (group_id, name.c_str (), p.dtype, dspace, p.plist);
            if (dset < 0)
              {
                bs_throw_exception (boost::format ("Can't create dataset %s in group %d") % name % group_id);
              }

            p.dspace = dspace;
            p.dset = dset;
            BS_ASSERT (p.dset >= 0) (name);
          }
        else 
          {
            hid_t dset = detail::open_dataset (group_id, name.c_str ());
            hid_t dspace = H5Dget_space (dset);
            if (dspace < 0)
              {
                bs_throw_exception (boost::format ("Can't get dataspace for dataset %s in group %d") % name % group_id);
              }

            p.dspace = dspace;
            p.dset = dset;
            BS_ASSERT (p.dset >= 0) (name);
          }
      }

    return p;
  }

  int 
  h5_pool::declare_i_data (const std::string &name, t_int def_value, int n_dims, npy_intp *dims, int var_dims)
  {
    declare_data (name, get_hdf5_type <hdf5_type_helper <t_int>::type> (), &def_value, n_dims, dims, var_dims);
    return 0;
  }

  int 
  h5_pool::declare_fp_data (const std::string &name, t_float def_value, int n_dims, npy_intp *dims, int var_dims)
  {
    declare_data (name, get_hdf5_type <hdf5_type_helper <t_float>::type> (), &def_value, n_dims, dims, var_dims);
    return 0;
  }

  void
  h5_pool::set_data (std::string const &name, hid_t dtype, void *data, t_long data_size, int n_dims, npy_intp *dims, void *def_value)
  {
    if (h5_map.find (name) == h5_map.end ())
      {
        declare_data (name, dtype, def_value, n_dims, dims);
      }

    const h5_pair &p = open_data (name);
    if (p.size > data_size)
      {
        bs_throw_exception (boost::format ("Size mismatch for %s: %d > %d") % name % p.size % data_size);
      }

    herr_t r = H5Dwrite (p.dset, p.dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    if (r < 0)
      {
        bs_throw_exception (boost::format ("Can't write data (%s)") % name);
      }
  }
    
  int 
  h5_pool::set_i_data (const std::string &name, spv_int data, t_int def_value)
  {
    set_data (name, 
              get_hdf5_type <hdf5_type_helper <t_int>::type> (), 
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
              get_hdf5_type <hdf5_type_helper <t_float>::type> (), 
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
    BS_ASSERT (ndims == 3) (ndims);
    for (int i = 0; i < 3; i++)
      pool_dims[i] = 0;
      
    n_pool_dims = ndims;
    for (int i = 0; i < ndims; i++)
      {
        pool_dims[i] = dims[i];
      }
    
    map_t::iterator it = h5_map.begin (), e = h5_map.end ();
    for (; it != e; ++it)
      {
        std::string const &name = it->first;
        h5_pair &p = it->second;
        if (p.var_dims)
          blue_sky::calc_data_dims (name, p, n_pool_dims, pool_dims);
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
           
          if (i->second.dset >= 0)
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
    
    for (i = 0; i < n_dims; ++i)
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

