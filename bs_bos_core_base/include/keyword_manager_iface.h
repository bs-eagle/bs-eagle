 /**
 * @file keywords.h
 * @brief interface for keyword_manager_iface class
 * @author Mark Khait
 * @date 2009-08-03
 * */
#ifndef KEYWORD_MANAGER_IFACE_H_
#define KEYWORD_MANAGER_IFACE_H_

#include "conf.h"

#include <boost/shared_ptr.hpp>

namespace blue_sky
{
  struct keyword_params //: boost::noncopyable
  {
    typedef smart_ptr <objbase, true>     sp_objbase;

    keyword_params ()
    {
    }

    keyword_params (const sp_objbase &km,
                    const sp_objbase &reader, 
                    const sp_objbase &data, 
                    const sp_objbase &em, 
                    const sp_objbase &mesh, 
                    const sp_objbase &fi_params,
                    const sp_objbase &scal_3p)
    : km (km),
    reader (reader),
    data (data),
    em (em),
    mesh (mesh),
    fi_params (fi_params),
    scal_3p (scal_3p)
    {
    }

    sp_objbase          km;
    sp_objbase          reader;
    sp_objbase          data;
    sp_objbase          em;
    sp_objbase          mesh;
    sp_objbase          fi_params;
    sp_objbase          scal_3p;
  };
      
  struct keyword_handler_iface
  {
    typedef keyword_params keyword_params_t;

    virtual ~keyword_handler_iface () {}

    virtual void
    handler (const std::string &, keyword_params_t &) = 0;
  };

  //! keyword_manager_iface - class-factory which contain a set of handlers for different keywords
  
  class BS_API_PLUGIN keyword_manager_iface: public objbase
    {
    public:
      //-----------------------------------------
      //  TYPES
      //-----------------------------------------
      
      typedef keyword_params            keyword_params_t;
      typedef keyword_handler_iface     keyword_handler_iface_t;

      typedef smart_ptr <objbase, true>               sp_objbase;
      
      //! type of pointer to function-handler keyword_manager_iface
      typedef void (*handler_t)(const std::string &, keyword_params_t &);
      typedef boost::shared_ptr <keyword_handler_iface_t> shared_handler_t;

      //! structure of keyword handler
      struct keyword_handler
      {
        keyword_handler ()
        : handle_function (0)
        , second_handle_function (0)
        , index_in_pool (-3)
        {
        }

        keyword_handler (handler_t handle_function)
        : handle_function (handle_function)
        , second_handle_function (0)
        , index_in_pool (-3)
        {
        }

        keyword_handler (handler_t handle_function, t_long index_in_pool)
        : handle_function (handle_function)
        , second_handle_function (0)
        , index_in_pool (index_in_pool)
        {
        }
        
        keyword_handler (handler_t handle_function, t_long def_value, int *new_dimens)
          : second_handle_function (handle_function)
          , index_in_pool (-1)
          , int_def_value (def_value)
          {
            dimens[0] = new_dimens[0];
            dimens[1] = new_dimens[1];
            dimens[2] = new_dimens[2];
            dimens[3] = new_dimens[3];
            dimens[4] = new_dimens[4];
            dimens[5] = new_dimens[5];
          }
        
        keyword_handler (handler_t handle_function, t_float def_value, t_long *new_dimens)
          : second_handle_function (handle_function)
          , index_in_pool (-2)
          , float_def_value (def_value)
          {
            dimens[0] = new_dimens[0];
            dimens[1] = new_dimens[1];
            dimens[2] = new_dimens[2];
            dimens[3] = new_dimens[3];
            dimens[4] = new_dimens[4];
            dimens[5] = new_dimens[5];
          }  

        keyword_handler (const shared_handler_t &handle_object)
        : handle_function (0)
        , handle_object (handle_object)
        , second_handle_function (0)
        , index_in_pool (-3)
        , int_def_value (0)
        , float_def_value (0)
        {
        }

        handler_t         handle_function;          //<! pointer to function
        shared_handler_t  handle_object;            //<! alternative for handle_function
        handler_t         second_handle_function;   //<! pointer to function
        
        t_long           index_in_pool;            //<! index in pool (for pooled keywords (pooled==handles by array handlers))
                                                    //<! -1: insert in imap (int) with available index
                                                    //<! -2: insert in dmap (float) with available index
                                                    //<! -3: invalid_value (because 0 is valid)
        t_long           int_def_value;            //<! default value (for int pooled keywords (pooled==handles by array handlers))
        t_float        float_def_value;          //<! default value (for float pooled keywords (pooled==handles by array handlers))
        int            dimens[6];                //<! dimensions of array (for pooled keywords (pooled==handles by array handlers))
      };

      

      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
    public:
      virtual ~keyword_manager_iface () {};
      
      //! registration of active keyword in factory
      virtual void register_keyword(const std::string &keyword, keyword_handler handler) = 0;
      
      //! registration of active integer pool keyword in factory
      virtual void register_i_pool_keyword(const std::string &keyword, int *dimens, t_long def_value, handler_t external_handler = 0) = 0;
      
      //! registration of active floating point pool keyword in factory
      virtual void register_fp_pool_keyword(const std::string &keyword, int *dimens, t_float def_value, handler_t external_handler = 0) = 0;
      
      //! registration of supported keywords in factory
      virtual void register_supported_keyword(const std::string &keyword, const std::string &provider) = 0;
    };
    
}//ns bs

#endif //KEYWORD_MANAGER_IFACE_H_

