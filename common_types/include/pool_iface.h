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
  class error_h5_no_error : public bs_exception
  {
  public:
    error_h5_no_error (std::string const &array_name);
  };

/** 
 * @brief 
 */
class h5_pool_iface : public bs_node
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
      virtual void open_file (const std::string &fname, const std::string &path) = 0;

      /** 
       * @brief close file
       */
      virtual void close_file () = 0;

      /** 
       * @brief flush data buffers
       */
      virtual void flush () const = 0;

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
       * 
       * @return 0 if success
       */
      virtual int set_fp_data (const std::string &name, spv_float data) = 0;

      /** 
       * @brief rewrite existing array in file
       * 
       * @param name    -- <INPUT> name of the array
       * @param data    -- <INPUT> smart pointer to the array 
       * 
       * @return 0 if success
       */
      virtual int set_i_data (const std::string &name, spv_int data) = 0;

#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const = 0;
#endif //BSPY_EXPORTING_PLUGIN
  };
} // end of blue_sky namespace

//#endif //USE_H%
#endif /* end of include guard: H_5_POOL_IFACE_XO0QK1SS */

