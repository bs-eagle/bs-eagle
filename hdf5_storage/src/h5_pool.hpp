/** 
 * @file h5_pool.hpp
 * @brief implementation of #h5_pool_iface
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-03-05
 */
#ifndef H5_POOL_G42LCHNB

#define H5_POOL_G42LCHNB

#include <map>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <string>
#include "pool_iface.h"

namespace blue_sky
{
  typedef boost::archive::text_iarchive           tia_t;
  typedef boost::archive::text_oarchive           toa_t;
  
  struct h5_pair
    {
      hid_t     dset;
      hid_t     dspace;
      hid_t     dtype;
      hid_t     plist;
      int       n_dims;
      t_long    size;
      npy_intp  py_dims[10];
      hsize_t   h5_dims[10];
      hsize_t   src_dims[6];  // dims for array are calculated from pool_dims and src_dims: dims[i] = src_dims[2 * i] * pool_dims[i] + src_dims[i + 1]. i = 0..2

      bool      var_dims;         // flag is true if array has variable dimensins depend on pool dimensions.
     /* 
      void save (toa_t &ar) const
        {
          ar & dset;
          ar & dspace;
          ar & dtype;
          ar & plist;
          ar & n_dims;
          ar & size;
          ar & py_dims;
          ar & h5_dims;
          ar & src_dims;
          ar & var_dims;
        }
        
      void load (tia_t &ar) const
        {
          ar & dset;
          ar & dspace;
          ar & dtype;
          ar & plist;
          ar & n_dims;
          ar & size;
          ar & py_dims;
          ar & h5_dims;
          ar & src_dims;
          ar & var_dims;
        }
        */
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
       * @brief set pool dims to calculate each pool array size
       * 
       * @param dims    -- <INPUT> dimensions
       * @param n_dim    -- <INPUT> dimension count (1, 2 or 3)
       */
      virtual void set_pool_dims (t_long *dims, int n_dim);

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
       * @brief calculate current data dimensions
       * 
       * @param name    -- <INPUT> given name of data
       * 
       * @return data size 
       */
      virtual t_long calc_data_dims (const std::string &name);

      /** 
       * @brief read data by name from file and return it as vector
       * 
       * @param name    -- <INPUT> given name of data
       * 
       * @return smart pointer to the data vector or
       * raise exception if no array with name 'name'
       */
      virtual spv_float get_fp_data (const std::string &name);

      /** 
       * @brief read data by name from file and return it as vector
       * 
       * @param name    -- <INPUT> given name of data
       * 
       * @return smart pointer to the data vector or
       * null pointer if no array with name 'name'
       */
      virtual spv_float get_fp_data_unsafe (const std::string &name);

      /** 
       * @brief read data by name from file and return it as vector
       * 
       * @param name    -- <INPUT> given name of data
       * 
       * @return smart pointer to the data vector or 
       * raise exception if no array with name 'name'
       */
      virtual spv_int get_i_data (const std::string &name);
      
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
       
      
      virtual int declare_fp_data (const std::string &name, t_float def_value, int n_dims, npy_intp *dims, int var_dims = 0);
      
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
       
      
      virtual int declare_i_data (const std::string &name, t_int def_value, int n_dims, npy_intp *dims, int var_dims = 0);

      /** 
       * @brief read data by name from file and return it as vector
       * 
       * @param name    -- <INPUT> given name of data
       * 
       * @return smart pointer to the data vector or 
       * null pointer if no array with name 'name'
       */
      virtual spv_int get_i_data_unsafe (const std::string &name);

      /** 
       * @brief rewrite existing array in file
       * 
       * @param name    -- <INPUT> name of the array
       * @param data    -- <INPUT> smart pointer to the array 
       * @param def_value  -- <INPUT> default value 
       * 
       * @return 0 if success
       */
       
      virtual int set_fp_data (const std::string &name, spv_float data, t_float def_value = 0);

      /** 
       * @brief rewrite existing array in file
       * 
       * @param name      -- <INPUT> name of the array
       * @param data      -- <INPUT> smart pointer to the array 
       * @param def_value  -- <INPUT> default value 
       * 
       * @return 0 if success
       */
      virtual int set_i_data (const std::string &name, spv_int data, t_int def_value = 0);

      /**
       * @brief returns is data array opened or not
       *
       * @param name -- <INPUT> name of the array
       * @return true if opened or false if not, throws exception 
       * if no array in pool
       * */
      virtual bool is_opened (const std::string &name);
      virtual std::string get_data_type (const std::string &name) const;

#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const;
      virtual boost::python::list py_list_data () const;
      virtual void py_set_pool_dims (boost::python::list &dims);
      virtual boost::python::list py_get_pool_dims () const;
      virtual boost::python::list py_get_data_dims (const std::string &name) const;
#endif //BSPY_EXPORTING_PLUGIN
    protected:

      /** 
       * @brief create new or replace existing array in file
       * 
       * @param name      -- <INPUT> name of the array
       * @param dtype     -- hdf5 type of data
       * @param def_value -- <INPUT> raw pointer to default value to fill the array 
       * @param n_dims    -- <INPUT> number of dims 
       * @param dims      -- <INPUT/OUTPUT> if var_dims = 1, source dimensions / calculated dimensions. if var_dims = 0, then just dimensions
       * @param var_dims  -- <INPUT> if var_dims = 1, array dimensions depend on pool dimensions. 
       * 
       * @return nothing, throws bs_exception on errors
       */
      void
      declare_data (std::string const &name, hid_t dtype, void *def_value, int n_dims, npy_intp *dims, int var_dims = 0);

      /** 
       * @brief rewrite existing array in file
       * 
       * @param name    -- <INPUT> name of the array
       * @param dtype   -- hdf5 type of data
       * @param data    -- <INPUT> raw pointer to data
       * @param data_size -- <INPUT> length of data
       * @param def_value  -- <INPUT> raw pointer to default value
       * 
       * @return 0 if success
       */
      void
      set_data (std::string const &name, 
                hid_t dtype, 
                void *data, 
                t_long data_size, 
                int n_dims,
                npy_intp *dims,
                void *def_value);

      /**
       * @brief opens data set
       *
       * @param name -- <INPUT> name of the array
       *
       * @return h5_pair, throws bs_exception on error
       * */
      h5_pair
      open_data (std::string const &name, const int h5_write = 0);


      /** 
       * @brief should be used only in constructor
       */
      void init_all ()
        {
          fname.clear ();
          path.clear ();
          file_id = 0; // FIXME: 0 is valid identifier
          group_id = 0; // FIXME: 0 is valid identifier

        }

      void fill_map ();
      void clear_map ();
      void close_node (h5_pair &p);

      template <class T>
      map_t::iterator add_node (const std::string &name, 
                                const hid_t dset, 
                                const hid_t dspace, 
                                const hid_t dtype, 
                                const hid_t plist, 
                                const int n_dims, 
                                const T *dims, 
                                const bool var_dims = 0);

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

      hsize_t           pool_dims[3];
      t_int             n_pool_dims;

      herr_t (*old_func)(hid_t, void*);
      void *old_client_data;
      int mute_flag;

      //blue-sky class declaration
      BLUE_SKY_TYPE_DECL (h5_pool);
    };

}//namespace blue_sky

#endif /* end of include guard: H5_POOL_G42LCHNB */
