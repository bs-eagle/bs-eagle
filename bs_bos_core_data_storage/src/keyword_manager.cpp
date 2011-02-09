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

namespace blue_sky
  {
  //constructors
  template <class strategy_t>
  keyword_manager<strategy_t>::keyword_manager(bs_type_ctor_param /*param*/)
  {
  }

  template <class strategy_t>
  keyword_manager<strategy_t>::keyword_manager(const keyword_manager <strategy_t>& src)
  : bs_refcounter (src)
  {
    //*this = src;
  }

  template <typename strategy_t>
  bool 
  keyword_manager <strategy_t>::is_keyword_supported (const std::string &keyword, keyword_params_t &params) const
  {
    sp_reader_t reader (params.reader, bs_dynamic_cast ());

    typename supported_keywords_t::const_iterator sup_it = supported_keywords.find (keyword);
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

  template <class strategy_t>
  bool
  keyword_manager <strategy_t>::is_keyword_activated (const std::string &keyword, keyword_params_t &params) const
  {
    return handlers.find(keyword) != handlers.end ();
  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::handle_keyword (const std::string &keyword, keyword_params_t &params)
  {
    typename handlers_t::iterator it = handlers.find(keyword);
    sp_reader_t reader (params.reader, bs_dynamic_cast ());
    if (it == handlers.end ())
      {
        is_keyword_supported (keyword, params);
      }
    else
      {
        keyword_handler &handler = it->second;
        if (handler.handle_function)
          {
            (*(handler.handle_function)) (keyword, params);
          }
        else if (handler.handle_object)
          {
            handler.handle_object->handler (keyword, params);
          }
      }
  }

  template <class strategy_t>
  void keyword_manager<strategy_t>::register_keyword(const std::string &keyword, keyword_handler handler)
  {
    typename handlers_t::iterator it = handlers.find(keyword);
    if (it != handlers.end())
      {
        bs_throw_exception(boost::format ("Keyword [%s] registration failed, keyword is already registered") % keyword);
      }
    switch (handler.index_in_pool)
      {
        case (-1) : handler.handle_function = &this_t::int_array_handler; break;
        case (-2) : handler.handle_function = &this_t::float_array_handler; break;
        default : break;
      }
      
    handlers.insert (std::make_pair (keyword, handler));
    BOSOUT (section::keywords, level::low) << boost::format ("Keyword [%s] registered") % keyword << bs_end;
  }

  template <typename strategy_t>
  void
  keyword_manager <strategy_t>::register_keyword (const std::string &keyword, const shared_handler_t &handler, bool replace_existing)
  {
    typename handlers_t::iterator it = handlers.find (keyword);
    if (it != handlers.end ())
      {
        if (replace_existing)
          handlers.erase (it);
        else
          bs_throw_exception(boost::format ("Keyword [%s] registration failed, keyword is already registered") % keyword);
      }

    handlers.insert (std::make_pair (keyword, keyword_handler (handler)));
    BOSOUT (section::keywords, level::low) << boost::format ("Keyword [%s] registered") % keyword << bs_end;
  }
  
  
  //! registration of active integer pool keyword in factory
  template <class strategy_t>
  void 
  keyword_manager <strategy_t>::register_i_pool_keyword(const std::string &keyword, int *dimens, i_type_t def_value, handler_t external_handler)
  {
    typename handlers_t::iterator it = handlers.find(keyword);
    if (it != handlers.end())
      {
        bs_throw_exception(boost::format ("Keyword [%s] registration failed, keyword is already registered") % keyword);
      }
    keyword_handler handler (&this_t::int_array_handler, def_value, dimens);
    handler.second_handle_function = external_handler;
      
    handlers.insert (std::make_pair (keyword, handler));
    BOSOUT (section::keywords, level::low) << boost::format ("Keyword [%s] registered") % keyword << bs_end;
  }
  
  //! registration of active floating point pool keyword in factory
  template <class strategy_t>
  void  keyword_manager<strategy_t>::register_fp_pool_keyword(const std::string &keyword, int *dimens, fp_type_t def_value, handler_t external_handler)
  {
    typename handlers_t::iterator it = handlers.find(keyword);
    if (it != handlers.end())
      {
        bs_throw_exception(boost::format ("Keyword [%s] registration failed, keyword is already registered") % keyword);
      }
    keyword_handler handler (&this_t::int_array_handler, def_value, dimens);
    handler.second_handle_function = external_handler;
      
    handlers.insert (std::make_pair (keyword, handler));
    BOSOUT (section::keywords, level::low) << boost::format ("Keyword [%s] registered") % keyword << bs_end;
  }
  
  template <class strategy_t>
  void  keyword_manager<strategy_t>::py_register_fp_pool_keyword (const std::string keyword, boost::python::list dimens, fp_type_t def_value)
    {
      int new_dimens[ARRAY_POOL_TOTAL];
      for (int i = 0; i < ARRAY_POOL_TOTAL; i++) 
        {
          new_dimens[i] = boost::python::extract<int>(dimens[i]);
        }
      register_fp_pool_keyword (keyword, &new_dimens[0], def_value);
    }
    
  template <class strategy_t>
  boost::python::list 
  keyword_manager<strategy_t>::py_list_active_keywords()
  {
    typename handlers_t::iterator it;
    boost::python::list items;
    for (it = handlers.begin(); it != handlers.end(); it++)
      items.append(it->first);
    return items;
  }
  
  template <class strategy_t>
  boost::python::list 
  keyword_manager<strategy_t>::py_list_supported_keywords()
  {
    supported_keywords_t::iterator it;
    boost::python::list items;
    for (it = supported_keywords.begin(); it != supported_keywords.end(); it++)
      items.append(it->first);
    return items;
  }
  
  template <class strategy_t>
  void keyword_manager<strategy_t>::register_supported_keyword(const std::string &keyword, const std::string &provider)
    {
      supported_keywords.insert (std::make_pair (keyword, std::list<std::string> ())).first->second.push_back(provider);
      BOSOUT (section::keywords, level::low) << boost::format ("Keyword [%s] supported") % keyword << bs_end;
    }

  /*template <class strategy_t>
  const typename keyword_manager<strategy_t>::sp_event_manager_t & keyword_manager<strategy_t>::get_event_manager()
  {
    return top->em;
  }*/

  //bs stuff
  BLUE_SKY_TYPE_STD_CREATE_T_DEF(keyword_manager, (class))
  //BLUE_SKY_TYPE_STD_COPY_T_DEF(keyword_manager, (class))
  
  template <typename strategy_t>
  blue_sky::objbase* 
  keyword_manager <strategy_t>::bs_create_copy (bs_type_cpy_ctor_param src) 
    {
      const keyword_manager <strategy_t> *src_ptr = dynamic_cast <const keyword_manager <strategy_t> *> (src.get ());
      if (!src_ptr)
        bs_throw_exception ("Can't cast to keyword_manager");

      return new keyword_manager <strategy_t> (src_ptr);
    }
    
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (keyword_manager <base_strategy_fif>) , 1, (objbase), "keyword_manager_fif", "BOS_Core keyword_manager class", "BOS_Core keyword_manager class", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (keyword_manager <base_strategy_did>) , 1, (objbase), "keyword_manager_did", "BOS_Core keyword_manager class", "BOS_Core keyword_manager class", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (keyword_manager <base_strategy_dif>) , 1, (objbase), "keyword_manager_dif", "BOS_Core keyword_manager class", "BOS_Core keyword_manager class", false);
  
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (keyword_manager <base_strategy_flf>) , 1, (objbase), "keyword_manager_flf", "BOS_Core keyword_manager class", "BOS_Core keyword_manager class", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (keyword_manager <base_strategy_dld>) , 1, (objbase), "keyword_manager_dld", "BOS_Core keyword_manager class", "BOS_Core keyword_manager class", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (keyword_manager <base_strategy_dlf>) , 1, (objbase), "keyword_manager_dlf", "BOS_Core keyword_manager class", "BOS_Core keyword_manager class", false);

  //!TODO: kill next string after debug
  //template class keyword_manager <base_strategy_di>;
  //template class keyword_manager <base_strategy_fi>;


}//ns_bs
