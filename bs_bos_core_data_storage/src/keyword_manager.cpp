/**
 *       \file  keyword_manager.cpp
 *      \brief  Implementation of keyword_manager
 *     \author  Morozov Andrey
 *       \date  02.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "bs_bos_core_data_storage_stdafx.h"
#include "keyword_manager.h"
#include "read_class.h"
#include "data_class.h"

namespace blue_sky
  {
  //constructors
  
  keyword_manager::keyword_manager(bs_type_ctor_param /*param*/)
  {
  }

  
  keyword_manager::keyword_manager(const keyword_manager & src)
  : bs_refcounter (src)
  {
    //*this = src;
  }

  
  bool 
  keyword_manager ::is_keyword_supported (const std::string &keyword, keyword_params_t &params) const
  {
    BS_SP (FRead) reader = params.hdm->get_reader ();

    supported_keywords_t::const_iterator sup_it = supported_keywords.find (keyword);
    if (sup_it == supported_keywords.end ())
      {  
        BOSWARN (section::keywords, level::warning) << (boost::format ("Keyword %s wasn't registered (%s)") % keyword % reader->get_prefix ()) << bs_end;
        return false;
      }
    else
      {
        std::list<std::string>::const_iterator i;
        BOSWARN (section::keywords, level::warning) << (boost::format ("Keyword %s is supported, but inactive now (%s). This keyword is supported by:") % keyword % reader->get_prefix ()) << bs_end;
        for (i = sup_it->second.begin(); i != sup_it->second.end(); ++i)
          BOSWARN (section::keywords, level::warning) << *i << bs_end;

        return true;
      }
  }

  
  bool
  keyword_manager ::is_keyword_activated (const std::string &keyword, keyword_params_t & /*params*/) const
  {
    return handlers.find(keyword) != handlers.end ();
  }

  
  void keyword_manager::handle_keyword (const std::string &keyword, keyword_params_t &params)
  {
    handlers_t::iterator it = handlers.find(keyword);
    BS_SP (FRead) reader = params.hdm->get_reader ();

    if (it == handlers.end ())
      {
        is_keyword_supported (keyword, params);
      }
    else
      {
        keyword_handler &handler = it->second;
        if (handler.read_handle_function)
          {
            (*(handler.read_handle_function)) (keyword, params);
          }
        if (handler.react_handle_function)
          {
            (*(handler.react_handle_function)) (keyword, params);
          }
        /*
        else if (handler.handle_object)
          {
            handler.handle_object->handler (keyword, params);
          }
        */
      }
  }
  
  void keyword_manager::handle_keyword_reactor (const std::string &keyword, keyword_params_t &params)
  {
    handlers_t::iterator it = handlers.find(keyword);
    BS_SP (FRead) reader = params.hdm->get_reader ();

    if (it == handlers.end ())
      {
        is_keyword_supported (keyword, params);
      }
    else
      {
        keyword_handler &handler = it->second;
        
        if (handler.react_handle_function)
          {
            (*(handler.react_handle_function)) (keyword, params);
          }
        /*
        else if (handler.handle_object)
          {
            handler.handle_object->handler (keyword, params);
          }
        */
      }
  }


  
  void keyword_manager::register_keyword(const std::string &keyword, keyword_handler handler)
  {
    handlers_t::iterator it = handlers.find(keyword);
    if (it != handlers.end())
      {
        bs_throw_exception(boost::format ("Keyword [%s] registration failed, keyword is already registered") % keyword);
      }
      
    handlers.insert (std::make_pair (keyword, handler));
    BOSOUT (section::keywords, level::low) << boost::format ("Keyword [%s] registered") % keyword << bs_end;
  }

  
/* 
  void
  keyword_manager ::register_keyword (const std::string &keyword, const shared_handler_t &handler, bool replace_existing)
  {
    handlers_t::iterator it = handlers.find (keyword);
    if (it != handlers.end ())
      {
        if (replace_existing)
          handlers.erase (it);
        else
          bs_throw_exception(boost::format ("Keyword [%s] registration failed, keyword is already registered") % keyword);
      }

    handlers.insert (std::make_pair (keyword, keyword_handler (handler, 0)));
    BOSOUT (section::keywords, level::low) << boost::format ("Keyword [%s] registered") % keyword << bs_end;
  }
  
  */

  //! registration of active integer pool keyword in factory
  
  void 
  keyword_manager ::register_i_pool_keyword(const std::string &keyword, npy_intp *dimens, t_int def_value, handler_t external_handler)
  {
    handlers_t::iterator it = handlers.find(keyword);
    if (it != handlers.end())
      {
        bs_throw_exception(boost::format ("Keyword [%s] registration failed, keyword is already registered") % keyword);
      }
    keyword_handler handler (&this_t::int_array_handler, def_value, dimens);
    handler.react_handle_function = external_handler;
      
    handlers.insert (std::make_pair (keyword, handler));
    hdm->get_pool()->declare_i_data (keyword, def_value, 3, dimens, 1);
    BOSOUT (section::keywords, level::low) << boost::format ("Keyword [%s] registered") % keyword << bs_end;
  }
  
  //! registration of active floating point pool keyword in factory
  
  void  keyword_manager::register_fp_pool_keyword (const std::string &keyword, npy_intp *dimens, t_float def_value, handler_t external_handler)
  {
    handlers_t::iterator it = handlers.find(keyword);
    if (it != handlers.end())
      {
        bs_throw_exception(boost::format ("Keyword [%s] registration failed, keyword is already registered") % keyword);
      }
    keyword_handler handler (&this_t::float_array_handler, def_value, dimens);
    handler.react_handle_function = external_handler;
    
    handlers.insert (std::make_pair (keyword, handler));
    hdm->get_pool()->declare_fp_data (keyword, def_value, 3, dimens, 1);
    BOSOUT (section::keywords, level::low) << boost::format ("Keyword [%s] registered") % keyword << bs_end;
  }

  void  keyword_manager::register_prop_keyword (const std::string &keyword, const std::string &format, prop_names_t &prop_names , handler_t external_handler)
  {
    handlers_t::iterator it = handlers.find(keyword);
    if (it != handlers.end())
      {
        bs_throw_exception(boost::format ("Keyword [%s] registration failed, keyword is already registered") % keyword);
      }
    
    keyword_handler handler (&this_t::prop_handler, format, prop_names);
    handler.react_handle_function = external_handler;
      
    handlers.insert (std::make_pair (keyword, handler));
  }
  
  void  keyword_manager::py_register_i_pool_keyword (const std::string keyword, boost::python::list dimens, t_int def_value)
    {
      npy_intp new_dimens[ARRAY_POOL_TOTAL];
      for (int i = 0; i < ARRAY_POOL_TOTAL; i++) 
        {
          new_dimens[i] = boost::python::extract<int>(dimens[i]);
        }
      register_i_pool_keyword (keyword, &new_dimens[0], def_value);
    }
  
  void  keyword_manager::py_register_fp_pool_keyword (const std::string keyword, boost::python::list dimens, t_float def_value)
    {
      npy_intp new_dimens[ARRAY_POOL_TOTAL];
      for (int i = 0; i < ARRAY_POOL_TOTAL; i++) 
        {
          new_dimens[i] = boost::python::extract<int>(dimens[i]);
        }
      register_fp_pool_keyword (keyword, &new_dimens[0], def_value);
    }

  
    
  
  boost::python::list 
  keyword_manager::py_list_active_keywords()
  {
    handlers_t::iterator it;
    boost::python::list items;
    for (it = handlers.begin(); it != handlers.end(); it++)
      items.append(it->first);
    return items;
  }
  
  
  boost::python::list 
  keyword_manager::py_list_supported_keywords()
  {
    supported_keywords_t::iterator it;
    boost::python::list items;
    for (it = supported_keywords.begin(); it != supported_keywords.end(); it++)
      items.append(it->first);
    return items;
  }
  
  
  void keyword_manager::register_supported_keyword(const std::string &keyword, const std::string &provider)
    {
      supported_keywords.insert (std::make_pair (keyword, std::list<std::string> ())).first->second.push_back(provider);
      BOSOUT (section::keywords, level::low) << boost::format ("Keyword [%s] supported") % keyword << bs_end;
    }

  //bs stuff
  BLUE_SKY_TYPE_STD_CREATE (keyword_manager)
  BLUE_SKY_TYPE_STD_COPY (keyword_manager)
  
  BLUE_SKY_TYPE_IMPL(keyword_manager, objbase, "keyword_manager", "BOS_Core keyword_manager class", "BOS_Core keyword_manager class");

  //!TODO: kill next string after debug
  //template class keyword_manager <base_strategy_di>;
  //template class keyword_manager <base_strategy_fi>;

}//ns_bs
