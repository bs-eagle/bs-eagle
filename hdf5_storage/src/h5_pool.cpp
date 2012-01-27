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
#include "date_helper.h"
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
      p.second.diff_from_base = true;
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
  h5_pool::fill_map (hid_t group)
  {
    BS_ASSERT (group);
    clear_map ();
    herr_t r = H5Literate (group, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, &(h5_pool::it_group), this);
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
  h5_pool::open_file (const std::string &fname_)
    {
	  //this function is not really open, it creates file
	  //TODO: set edit_base = false; at open_file
      edit_base = true;

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

      //TODO: Open existing file: remove H5F_ACC_TRUNC, parse and add existing groups
      group_id.clear ();
      //std::string base_name = "base_" + get_date_time_str();
      //group_id.insert (pair<std::string, hid_t>(base_name.c_str(), -1));
      group_id.insert (pair<std::string, hid_t>("actual", -1));

      map_hid_t::iterator i;
      for (i = group_id.begin(); i != group_id.end(); ++i)
        {
          std::string path_ = i->first;
          if (detail::is_object_exists (file_id, path_))
            {
        	  i->second = H5Gopen (file_id, path_.c_str ());
              if (i->second < 0)
                {
                  bs_throw_exception (boost::format ("Can't open existing group: %s") % path_);
                }
            }
          else
            {
        	  i->second = H5Gcreate (file_id, path_.c_str (), -1);
              if (i->second < 0)
                {
                  bs_throw_exception (boost::format ("Can't create group: %s") % path_);
                }
            }
          fill_map (i->second);
        }
    }
  
  void 
  h5_pool::close_file ()
    {
      clear_map ();
      map_hid_t::iterator i;
      for (i = group_id.begin(); i != group_id.end(); ++i)
        {
          if (i->second > 0)
            {
              H5Gclose (i->second);
            }
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
	std::string group_name = "";
	//std::cout<<"get_fp_data ( " << name.c_str () <<" )\n";
	//std::cout<<"try actual\n";
	map_hid_t::iterator group = group_id.find ("actual");
    if (detail::is_object_exists (group->second, name.c_str ()))
      {
		group_name = group->first;
    	//std::cout<<"found at "<< group_name.c_str () <<"\n";
        return get_fp_data_group (name, group_name);
      }
    map_hid_t::iterator b = group_id.begin();
	map_hid_t::reverse_iterator rev_end (++b), rev_i (group_id.end());
	for (;rev_i != rev_end; ++rev_i)
	  {
		group_name = rev_i->first;
	    //std::cout<<"try "<< group_name.c_str () <<"\n";
		group = group_id.find (group_name);
        if (detail::is_object_exists (group->second, name.c_str ()))
	      break;
	  }
	//if (rev_i != rev_end)
	//  std::cout<<"found at "<< group_name.c_str () <<"\n";
	//else
	//  std::cout<<"not found !!!\n";
    return get_fp_data_group (name, group_name);
  }

  spv_float
  h5_pool::get_fp_data_group (std::string const &name, const std::string &group_name)
  {
    spv_float a = get_fp_data_unsafe_group (name, group_name);
    if (!a || a->empty ())
      throw error_h5_no_array (name);

    return a;
  }

  spv_int
  h5_pool::get_i_data (std::string const &name)
  {
	std::string group_name = "";
	//std::cout<<"get_i_data ( " << name.c_str () <<" )\n";
	//std::cout<<"try actual\n";
	map_hid_t::iterator group = group_id.find ("actual");
    if (detail::is_object_exists (group->second, name.c_str ()))
      {
		group_name = group->first;
    	//std::cout<<"found at "<< group_name.c_str () <<"\n";
        return get_i_data_group (name, group_name);
      }
    map_hid_t::iterator b = group_id.begin();
	map_hid_t::reverse_iterator rev_end (++b), rev_i (group_id.end());
	for (;rev_i != rev_end; ++rev_i)
	  {
		group_name = rev_i->first;
	    group = group_id.find (group_name);
	    //std::cout<<"try "<< group_name.c_str () <<" id "<<group->second<< "\n";
        if (detail::is_object_exists (group->second, name.c_str ()))
	      break;
	  }
	//if (rev_i != rev_end)
	//  std::cout<<"found at "<< group_name.c_str () <<"\n";
	//else
	//  std::cout<<"not found !!!\n";
    return get_i_data_group (name, group_name);
  }

  spv_int
  h5_pool::get_i_data_group (std::string const &name, const std::string &group_name)
  {
    spv_int a = get_i_data_unsafe_group (name, group_name);
    if (!a || a->empty ())
      throw error_h5_no_array (name);

    return a;
  }

  bool
  h5_pool::is_opened_group (std::string const &name, const std::string &group_name)
  {
	hid_t group = group_id.find (group_name)->second;
    BS_ASSERT (group >= 0) (name);
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
  h5_pool::get_fp_data_unsafe_group (const std::string &name, const std::string &group_name)
  {
    // FIXME: was <=
	hid_t group = group_id.find (group_name)->second;
    if (group < 0)
      {
        bs_throw_exception (boost::format ("Get data %s: group not opened") % name);
      }
    if (h5_map.find (name) == h5_map.end ())
      {
        return spv_float ();
      }

    return get_data <v_float> (name, open_data (name, group_name));
  }

  spv_int 
  h5_pool::get_i_data_unsafe_group (const std::string &name, const std::string &group_name)
  {
    // FIXME: was <=
	hid_t group = group_id.find (group_name)->second;
    if (group < 0)
      {
        bs_throw_exception (boost::format ("Get data %s: group not opened") % name);
      }
    if (h5_map.find (name) == h5_map.end ())
      {
        return spv_int ();
      }

    return get_data <v_int> (name, open_data (name, group_name));
  }
    
    
  void
  h5_pool::declare_data (std::string const &name, hid_t dtype, void *value, int n_dims, npy_intp *dims, int var_dims)
  {
    BOOST_STATIC_ASSERT (sizeof (t_long) >= sizeof (npy_intp));
    //BOOST_STATIC_ASSERT (sizeof (t_long) >= sizeof (hsize_t));

    hid_t group = group_id.find ("actual")->second;

    map_t::iterator it = h5_map.find (name);
    if (it != h5_map.end ())
      {
        BOSOUT (section::h5, level::warning) << "Declared data " << name << " already exists in pool" << bs_end;
        /*
        h5_pair &p = it->second;

        for (int i = 0; i < 6; ++i)
          {
            p.src_dims[i] = dims[i];
          }
        */
      }
    else
      {
        hid_t plist = H5Pcreate (H5P_DATASET_CREATE);
        if (plist < 0)
          {
            bs_throw_exception (boost::format ("Can't create property for dataset %s in group %d") % name % group);
          }

        herr_t r = H5Pset_fill_value (plist, dtype, value);
        if (r < 0)
          {
            bs_throw_exception (boost::format ("Can't set fill for property, dataset %s in group %d") % name % group);
          }

        hid_t dtype_copy = H5Tcopy (dtype);
        if (dtype_copy < 0)
          {
            bs_throw_exception (boost::format ("Can't copy datatype for %s in group %d") % name % group);
          }

        add_node (name, -1, -1, dtype_copy, plist, n_dims, dims, var_dims);
      }
  }

  h5_pair
  h5_pool::open_data (std::string const &name, const std::string &group_name)
  {
	hid_t group = group_id.find (group_name)->second;
	hid_t group_actual = group_id.find ("actual")->second;

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
            H5Ldelete (group, name.c_str (), H5P_DEFAULT);
          }

        hid_t dtype = p.dtype;
        hid_t plist = p.plist;

        p.dtype = p.plist = -1; // prevent closing of dtype and plist
        close_node (p);

        p.dtype = dtype;
        p.plist = plist;

        blue_sky::calc_data_dims (name, p, n_pool_dims, pool_dims);
      }

    //if (p.dset < 0)
      {
        if (!detail::is_object_exists (group, name.c_str ()))
          {
            hid_t dspace = H5Screate_simple (p.n_dims, p.h5_dims, NULL);
            if (dspace < 0)
              {
                bs_throw_exception (boost::format ("Can't create simple dataspace for dataset %s in group %d") % name % group);
              }

            BS_ASSERT (p.dtype >= 0) (name);
            BS_ASSERT (p.plist) (name);

            hid_t dset = H5Dcreate (group, name.c_str (), p.dtype, dspace, p.plist);
            if (dset < 0)
              {
                bs_throw_exception (boost::format ("Can't create dataset %s in group %d") % name % group);
              }

            p.dspace = dspace;
            p.dset = dset;
            BS_ASSERT (p.dset >= 0) (name);

            //H5Lcreate_soft(("/actual/" + name).c_str (), group_actual, name.c_str ());
          }
        else 
          {
            hid_t dset = detail::open_dataset (group, name.c_str ());
            hid_t dspace = H5Dget_space (dset);
            if (dspace < 0)
              {
                bs_throw_exception (boost::format ("Can't get dataspace for dataset %s in group %d") % name % group);
              }

            p.dspace = dspace;
            p.dset = dset;
            BS_ASSERT (p.dset >= 0) (name);
          }
      }

    return p;
  }

  const char *
  h5_pool::get_script ()
  {
	const hsize_t rank = 1;
	hid_t dset, dspace;

	// Create the datatype
	hid_t dtype = H5Tcopy (H5T_C_S1);
	H5Tset_size (dtype, BUF_SIZE);

    if (detail::is_object_exists (file_id, "script"))
      {
        dset = detail::open_dataset (file_id, "script");
        herr_t r = H5Dread (dset, dtype, NULL, NULL, H5P_DEFAULT, buf);
      }
    else
      return "";
    return buf;
  }

  int 
  h5_pool::add_script (const std::string &script)
  {
	const hsize_t rank = 1;
	hid_t dset, dspace;
    std::string script2 = "";

	// create the datatype
	hid_t dtype = H5Tcopy (H5T_C_S1);
	H5Tset_size (dtype, BUF_SIZE);

	// create "script" if it not exists
	bool is_exist = detail::is_object_exists (file_id, "script");
    if (!is_exist)
      {
	    const hsize_t dims = 1;
	    const hsize_t max_dims = 1;
        dspace = H5Screate_simple (1, &dims, &max_dims);
	    dset = H5Dcreate (file_id, "script", dtype, dspace, H5P_DEFAULT);
      }
    //else// read "script" and store to script2
    //  {
    //    dset = detail::open_dataset (file_id, "script");
    //    herr_t r = H5Dread (dset, dtype, NULL, NULL, H5P_DEFAULT, buf);
    //    script2 = std::string (buf);
    //  }

    //script2 = script2 + "\n# " + get_date_time_str () +"\n" + script;
    //script2 = "\n# " + get_date_time_str () +"\n" + script;
    script2 = script;

    sprintf (buf, "%s", script2.c_str ());
    herr_t r = H5Dwrite (dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);

    return 0;
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

#ifdef BSPY_EXPORTING_PLUGIN
  int 
  h5_pool::py_declare_i_data (const std::string &name, t_int def_value, int n_dims, boost::python::list &dims, int var_dims)
  {
    t_int i;
    npy_intp arr_dims[10];

    if (var_dims)
      for (i = 0; i < 6; ++i)
        {
          arr_dims[i] = boost::python::extract<int>(dims[i]);
        }
    else
      for (i = 0; i < n_dims; ++i)
        {
          arr_dims[i] = boost::python::extract<int>(dims[i]);
        }
      
    declare_data (name, get_hdf5_type <hdf5_type_helper <t_int>::type> (), &def_value, n_dims, arr_dims, var_dims);
    return 0;
  }

  int 
  h5_pool::py_declare_fp_data (const std::string &name, t_float def_value, int n_dims, boost::python::list &dims, int var_dims)
  {
    t_int i;
    npy_intp arr_dims[10];

    if (var_dims)
      for (i = 0; i < 6; ++i)
        {
          arr_dims[i] = boost::python::extract<int>(dims[i]);
        }
    else
      for (i = 0; i < n_dims; ++i)
        {
          arr_dims[i] = boost::python::extract<int>(dims[i]);
        }
      
    declare_data (name, get_hdf5_type <hdf5_type_helper <t_float>::type> (), &def_value, n_dims, arr_dims, var_dims);
    return 0;
  }
#endif //BSPY_EXPORTING_PLUGIN

  void
  h5_pool::set_data (std::string const &name, hid_t dtype, void *data, t_long data_size, int n_dims, npy_intp *dims, void *def_value, bool create_base)
  {
	std::string last_base_name = "actual";
	map_t::iterator map_it = h5_map.find (name);
    if (map_it == h5_map.end ()) // new array
      {
        declare_data (name, dtype, def_value, n_dims, dims);
      }
    else if (create_base)// already exist and create_base=true
      {
		if (map_it->second.diff_from_base) // array in last base different from prev base
    	  {
    	    //create new base group
    	    last_base_name = "base_" + get_date_time_str();

    	    hid_t base_id = H5Gcreate (file_id, last_base_name.c_str(), -1);
    	    group_id.insert (pair<std::string, hid_t>(last_base_name, base_id));

    	    std::string comment = "# Created base\n";
    	    add_script (comment);

    	    // set other array's flags difference from (new) base is false
    	    map_t::iterator it = h5_map.begin (), e = h5_map.end ();
            for (; it != e; ++it)
              {
            	if (it != map_it)
                  it->second.diff_from_base = false;
              }
          }
    	else
    	  {
    	    //set array path to last base instead of "actual"
    	    map_hid_t::iterator e = group_id.end();
    	    if (!group_id.empty ())
    	      last_base_name = (--e)->first;
    	    else
    	      bs_throw_exception (boost::format ("set data %s error: No groups in pool!") % name);
    	    //set flag current array is changed
        	if (!edit_base)
        	  map_it->second.diff_from_base = true;
    	  }
      }

    const h5_pair &p = open_data (name, last_base_name);
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

  void
  h5_pool::finish_base ()
  {
	edit_base = false;
    // set flag difference from (new) base is true for all arrays
	// next call of set_data will create new base
    map_t::iterator it = h5_map.begin (), e = h5_map.end ();
    for (; it != e; ++it)
      {
        it->second.diff_from_base = true;
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
              &def_value, true);
    return 0;
  }

  int
  h5_pool::set_i_data_script (const std::string &name, spv_int data, t_int def_value)
  {
    set_data (name,
              get_hdf5_type <hdf5_type_helper <t_int>::type> (),
              data->data (),
              data->size (),
              data->ndim (),
              data->dims (),
              &def_value, false);
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
              &def_value, true);
    return 0;
  }

  int
  h5_pool::set_fp_data_script (const std::string &name, spv_float data, t_float def_value)
  {
    set_data (name,
              get_hdf5_type <hdf5_type_helper <t_float>::type> (),
              data->data (),
              data->size (),
              data->ndim (),
              data->dims (),
              &def_value, false);
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
  h5_pool::get_data_type_group (const std::string &name, const std::string &group_name) const
  {
	hid_t group = group_id.find (group_name)->second;

    map_t::const_iterator i;
    H5T_class_t dt;

    // FIXME: was <=
    if (group < 0)
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
  h5_pool::py_get_pool_dims () const
  {
    boost::python::list dims;
    
    for (int i = 0; i < n_pool_dims; ++i)
      {
        dims.append(pool_dims[i]);
      }
    return dims;
  }

  boost::python::list 
  h5_pool::py_get_data_dims (const std::string &name) const
  {
    boost::python::list dims;
    map_t::const_iterator ti, te; 
    ti = h5_map.find (name);
    if (ti != h5_map.end ())
      {
        int n = ti->second.n_dims;
        for (int i = 0; i < n; ++i)
          {
            dims.append(ti->second.py_dims[i]);
          }

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

