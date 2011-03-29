/** 
 * @file h5_pool.h
 * @brief implementation of #h5_pool_iface
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-03-05
 */
#ifndef H5_POOL_G42LCHNB

#define H5_POOL_G42LCHNB

#include <map>
#include <string>
#include "pool_iface.h"

namespace blue_sky
{
  struct h5_pair
    {
      hid_t dset;
      hid_t dspace;
      hid_t dtype;
      int n_dims;
      npy_intp py_dims[10];
      hsize_t h5_dims[10];
    };

  /** 
   * @brief interface class for block CSR matrix storage and manipulation
   */
  class BS_API_PLUGIN h5_pool: public h5_pool_iface
    {

    public:

      typedef h5_pool_iface                                     base_t;
      typedef h5_pool                                           this_t;
      typedef std::map<std::string, h5_pair>                    map_t;
      typedef std::pair<std::string, h5_pair>                   pair_t;


    public:
      
      //! destructor
      virtual ~h5_pool ()
        {};

      /** 
       * @brief create file for reading and writing 
       * 
       * @param fname   -- <INPUT> file name
       * @param path    -- <INPUT> path to the pool
       */
      virtual void open_file (const std::string &fname, const std::string &path);

      /** 
       * @brief close file
       */
      virtual void close_file ();

      /** 
       * @brief flush data buffers
       */
      virtual void flush () const
        {
          if (file_id > 0)
            {
              H5Fflush (file_id, H5F_SCOPE_LOCAL);
            }
        }

      /** 
       * @brief read data by name from file and return it as vector
       * 
       * @param name    -- <INPUT> given name of data
       * 
       * @return smart pointer to the data vector
       */
      virtual spv_float get_fp_data (const std::string &name);

      /** 
       * @brief read data by name from file and return it as vector
       * 
       * @param name    -- <INPUT> given name of data
       * 
       * @return smart pointer to the data vector
       */
      virtual spv_int get_i_data (const std::string &name);

      /** 
       * @brief rewrite existing array in file
       * 
       * @param name    -- <INPUT> name of the array
       * @param data    -- <INPUT> smart pointer to the array 
       * 
       * @return 0 if success
       */
      virtual int set_fp_data (const std::string &name, spv_float data);

      /** 
       * @brief rewrite existing array in file
       * 
       * @param name    -- <INPUT> name of the array
       * @param data    -- <INPUT> smart pointer to the array 
       * 
       * @return 0 if success
       */
      virtual int set_i_data (const std::string &name, spv_int data);

#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const;
#endif //BSPY_EXPORTING_PLUGIN
    public:
    
    protected:

      /** 
       * @brief should be used only in constructor
       */
      void init_all ()
        {
          fname.clear ();
          path.clear ();
          file_id = 0;
          group_id = 0;

        }

      void fill_map ();
      void clear_map ();
      void close_node (h5_pair &p);
      template <class T>
      map_t::iterator add_node (const std::string &name, const hid_t dset, const hid_t dspace, 
                     const hid_t dtype, const int n_dims, const T *dims);
      static herr_t it_group (hid_t g_id, const char *name, const H5L_info_t *info, 
                              void *op_data); 

      void mute ();
      void un_mute ();

    //______________________________________
    //  VARIABLES
    //______________________________________
    public:
      std::string       fname;          //!< h5 file name
      hid_t             file_id;        //!< h5 file identifier
      std::string       path;           //!< path of the pool in h5 file
      hid_t             group_id;       //!< base group in file

   protected:
      map_t             h5_map;

      herr_t (*old_func)(hid_t, void*);
      void *old_client_data;
      int mute_flag;

      //blue-sky class declaration
      BLUE_SKY_TYPE_DECL (h5_pool);
    };

}//namespace blue_sky

#endif /* end of include guard: H5_POOL_G42LCHNB */
