 /**
 * @file keywords.h
 * @brief interface for keyword_manager_iface class
 * @author Mark Khait
 * @date 2009-08-03
 * */
#ifndef KEYWORD_MANAGER_IFACE_H_
#define KEYWORD_MANAGER_IFACE_H_

#include "conf.h"
#include "hydrodynamic_model_iface.h"

#include <boost/shared_ptr.hpp>

namespace blue_sky
{
  struct keyword_params //: boost::noncopyable
  {
    typedef smart_ptr <hydrodynamic_model_iface, true>     sp_hdm;

    keyword_params ()
    {
    }

    keyword_params (sp_hdm hdm)
    : hdm (hdm)
    {
    }

    sp_hdm          hdm;
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
      typedef std::vector<std::string>  prop_names_t;
      
      //! type of pointer to function-handler keyword_manager_iface
      typedef void (*handler_t)(const std::string &, keyword_params_t &);
      typedef boost::shared_ptr <keyword_handler_iface_t> shared_handler_t;

      //! structure of keyword handler
      struct keyword_handler
      {
        keyword_handler ()
        : read_handle_function (0)
        , react_handle_function (0)
        {
        }

        keyword_handler (handler_t reader, handler_t reactor)
        : read_handle_function (reader)
        , react_handle_function (reactor)
        {
        }

        keyword_handler (handler_t read_handle_function, t_int def_value, t_int *new_dimens)
          : read_handle_function (read_handle_function)
          , int_def_value (def_value)
          {
            dimens[0] = new_dimens[0];
            dimens[1] = new_dimens[1];
            dimens[2] = new_dimens[2];
            dimens[3] = new_dimens[3];
            dimens[4] = new_dimens[4];
            dimens[5] = new_dimens[5];
          }
        
        keyword_handler (handler_t read_handle_function, t_float def_value, t_int *new_dimens)
          : read_handle_function (read_handle_function)
          , float_def_value (def_value)
          {
            dimens[0] = new_dimens[0];
            dimens[1] = new_dimens[1];
            dimens[2] = new_dimens[2];
            dimens[3] = new_dimens[3];
            dimens[4] = new_dimens[4];
            dimens[5] = new_dimens[5];
          }  

        keyword_handler (handler_t read_handle_function, std::string format, prop_names_t names)
          : read_handle_function (read_handle_function)
          , prop_names (names)
          , prop_format (format)
          {
          }  

        /*       
        keyword_handler (const shared_handler_t &handle_object)
        : read_handle_function (0)
        , handle_object (handle_object)
        , react_handle_function (0)
        , int_def_value (0)
        , float_def_value (0)
        {
        }
        */

        handler_t         read_handle_function;          //<! pointer to function
        //shared_handler_t  handle_object;            //<! alternative for read_handle_function
        handler_t         react_handle_function;   //<! pointer to function
        
        t_int        int_def_value;            //<! default value (for int pooled keywords (pooled==handles by array handlers))
        t_float      float_def_value;          //<! default value (for float pooled keywords (pooled==handles by array handlers))
        t_int        dimens[6];                //<! dimensions of array (for pooled keywords (pooled==handles by array handlers))
        prop_names_t prop_names;               //<! names of properties for prop keyword
        std::string  prop_format;              //<! format for property keyword
      };

      

      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
    public:
      virtual ~keyword_manager_iface () {};

      //! register all plugins keywords 
      virtual void init() = 0;
      
      //! registration of active keyword in factory
      virtual void register_keyword(const std::string &keyword, keyword_handler handler) = 0;
      
      //! registration of active integer pool keyword in factory
      virtual void register_i_pool_keyword(const std::string &keyword, int *dimens, t_int def_value, handler_t external_handler = 0) = 0;
      
      //! registration of active floating point pool keyword in factory
      virtual void register_fp_pool_keyword(const std::string &keyword, int *dimens, t_float def_value, handler_t external_handler = 0) = 0;

      //! registration of active property(-ies) keyword in factory
      virtual void register_prop_keyword(const std::string &keyword, const std::string &format, prop_names_t &prop_names , handler_t external_handler = 0) = 0;

      //! registration of supported keywords in factory
      virtual void register_supported_keyword(const std::string &keyword, const std::string &provider) = 0;

      //! launch keyword handler
      virtual void handle_keyword (const std::string &keyword, keyword_params_t &params) = 0;
    };
    
}//ns bs

#endif //KEYWORD_MANAGER_IFACE_H_

