/** 
 * @file h5_pool_iface.h
 * @brief interface class for data pool (storing int and float arrays)
 * @author Oleg Borschuk
 * @version 0.1
 * @date 2011-02-23
 */
#ifndef H5_POOL_IFACE_XO0QK1SS

#define H5_POOL_IFACE_XO0QK1SS

//#ifdef USE_H5

#include <string>
#include <hdf5.h>

#include "bs_assert.h"
#include "bs_tree.h"
#include "smart_ptr.h"
#include "bs_array.h"
#include "conf.h"


namespace blue_sky
{
  /**
   * @brief raise if no array in H5 pool
   */
  class error_h5_no_array : public bs_exception
  {
  public:
    error_h5_no_array (std::string const &array_name);
  };

/** 
 * @brief 
 */
class h5_pool_iface : public objbase
  {
    // --------------------------------
    // METHODS
    // --------------------------------
    public:
      /** 
       * @brief destructor
       */
      virtual ~h5_pool_iface ()
        {}
      
      /** 
       * @brief create file for reading and writing 
       * 
       * @param fname   -- <INPUT> file name
       * @param path    -- <INPUT> path to the pool
       */
      virtual void open_file (const std::string &fname) = 0;

      /** 
       * @brief close file
       */
      virtual void close_file () = 0;

      /** 
       * @brief flush data buffers
       */
      virtual void flush () const = 0;
      
      /** 
       * @brief set pool dims to calculate each pool array size
       * 
       * @param dims    -- <INPUT> dimensions
       * @param n_dim    -- <INPUT> dimension count (1, 2 or 3)
       */
      virtual void set_pool_dims (t_long *dims, int n_dim) = 0;

      /** 
       * @brief read data by name from file and return it as vector
       * 
       * @param name    -- <INPUT> given name of data
       * 
       * @return smart pointer to the data vector or
       * raise exception if no array with name 'name'
       */
      virtual spv_float get_fp_data (const std::string &name) = 0;

      /** 
       * @brief read data by name from file and return it as vector
       * 
       * @param name    -- <INPUT> given name of data
       * 
       * @return smart pointer to the data vector or
       * null pointer if no array with name 'name'
       */
      virtual spv_float get_fp_data_unsafe (const std::string &name) = 0;

      /** 
       * @brief read data by name from file and return it as vector
       * 
       * @param name    -- <INPUT> given name of data
       * 
       * @return smart pointer to the data vector or
       * raise exception if no array with name 'name'
       */
      virtual spv_int get_i_data (const std::string &name) = 0;

      /** 
       * @brief read data by name from file and return it as vector
       * 
       * @param name    -- <INPUT> given name of data
       * 
       * @return smart pointer to the data vector or
       * null pointer if no array with name 'name'
       */
      virtual spv_int get_i_data_unsafe (const std::string &name) = 0;

      /** 
       * @brief rewrite existing array in file
       * 
       * @param name    -- <INPUT> name of the array
       * @param data    -- <INPUT> smart pointer to the array 
       * @param src_dims  -- <INPUT> source dimensions to calculate array dimensions using pool_dims 
       * 
       * @return 0 if success
       */
       
      virtual int declare_fp_data (const std::string &name, t_float def_value, int n_dims, npy_intp *dims, int var_dims = 0) = 0;
      
        /** 
       * @brief create new or replace existing array in file
       * 
       * @param name      -- <INPUT> name of the array
       * @param def_value -- <INPUT> default value to fill the array 
       * @param n_dims    -- <INPUT> number of dims 
       * @param dims      -- <INPUT/OUTPUT> if var_dims = 1, source dimensions / calculated dimensions. if var_dims = 0, then just dimensions
       * @param var_dims  -- <INPUT> if var_dims = 1, array dimensions depend on pool dimensions. 
       * 
       * @return 0 if success
       */
       
      
      virtual int declare_i_data (const std::string &name, t_int def_value, int n_dims, npy_intp *dims, int var_dims = 0) = 0;

     /** 
       * @brief rewrite existing array in file
       * 
       * @param name    -- <INPUT> name of the array
       * @param data    -- <INPUT> smart pointer to the array 
       * @param src_dims  -- <INPUT> source dimensions to calculate array dimensions using pool_dims 
       * 
       * @return 0 if success
       */
      
      
      /** 
       * @brief calculate current data dimensions
       * 
       * @param name    -- <INPUT> given name of data
       * 
       * @return data size 
       */
      virtual t_long calc_data_dims (const std::string &name) = 0;
       
      virtual int set_fp_data (const std::string &name, spv_float data, t_float def_value = 0) = 0;
      virtual int set_fp_data_script (const std::string &name, spv_float data, t_float def_value = 0) = 0;

      /** 
       * @brief rewrite existing array in file
       * 
       * @param name      -- <INPUT> name of the array
       * @param data      -- <INPUT> smart pointer to the array 
       * 
       * @return 0 if success
       */
      virtual int set_i_data (const std::string &name, spv_int data, t_int def_value = 0) = 0;
      virtual int set_i_data_script (const std::string &name, spv_int data, t_int def_value = 0) = 0;

      virtual void finish_base () = 0;
      virtual int add_script (const std::string &script, bool replace) = 0;
      virtual const char * get_script () = 0;
      /**
       * @brief returns is data array opened or not
       *
       * @param name -- <INPUT> name of the array
       * @return true if opened or false if not, throws exception 
       * if no array in pool
       * */
      virtual bool is_opened (const std::string &name) = 0;
      
      /**
       * @brief returns data array type
       *
       * @param name -- <INPUT> name of the array
       * @return string name of type, throws exception 
       * if no array in pool
       * */
      virtual std::string get_data_type(const std::string &name) const = 0;


#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const = 0;
      virtual boost::python::list py_list_data () const = 0;
      virtual void py_set_pool_dims (boost::python::list &dims) = 0;
      virtual boost::python::list py_get_pool_dims () const = 0;
      virtual boost::python::list py_get_data_dims (const std::string &name) const = 0;

      virtual int py_declare_i_data (const std::string &name, t_int def_value, int n_dims, boost::python::list &dims, int var_dims) = 0;
      virtual int py_declare_fp_data (const std::string &name, t_float def_value, int n_dims, boost::python::list &dims, int var_dims) = 0;
#endif //BSPY_EXPORTING_PLUGIN
  };
} // end of blue_sky namespace

//#endif //USE_H%
#endif /* end of include guard: H_5_POOL_IFACE_XO0QK1SS */

