/**
 *       \file  event_base.cpp
 *      \brief  Implementation of event_base
 *     \author  Morozov Andrey
 *       \date  07.06.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"
#include <boost/spirit.hpp>
#include <boost/spirit/core.hpp>

#include "event_base.h"
#include "str_functor.h"

using namespace boost::spirit;

namespace blue_sky {

  namespace detail {

    /**
     * \class str_functor
     * \brief Handles actions in boost::spirit parser
     * */
    template <typename T>
    struct str_functor
    {
      typedef smart_ptr <named_pbase, true> sp_named_pbase;
      typedef void (T::*method_t)(const sp_named_pbase &params_, const char*, const char*);
      typedef T self_t;

      str_functor (self_t *self, const sp_named_pbase &params, method_t method)
          : self (self)
          , params_ (params)
          , method (method)
      {
      }

      void operator() (const char *first, const char *last) const
        {
          (self->*method) (params_, first, last);
        }

      self_t          *self;
      sp_named_pbase  params_;
      method_t        method;
    };

    /**
     * \class char_functor
     * \brief Handles actions in boost::spirit parser
     * */
    template <typename T>
    struct char_functor
    {
      typedef smart_ptr <named_pbase, true> sp_named_pbase;
      typedef void (T::*method_t)(const sp_named_pbase &params_, const char);
      typedef T self_t;

      char_functor (self_t *self, const sp_named_pbase &params, method_t method)
          : self (self)
          , params_ (params)
          , method (method)
      {
      }

      void operator() (const char c) const
        {
          (self->*method) (params_, c);
        }

      self_t          *self;
      sp_named_pbase  params_;
      method_t        method;
    };
  }


  template <typename strategy_t>
  event_base<strategy_t>::~event_base ()
  {

  }

  /**
   * \brief  'default' ctor for event_base
   * \param  param Additional params for ctor
   * */
  template <typename strategy_t>
  event_base<strategy_t>::event_base(bs_type_ctor_param /*param*/)
  {
    num   = 1;
    index = 1;
  }
  /**
   * \brief  copy-ctor for event_base
   * \param  src Instance of event_base to be copied
   * */
  template <typename strategy_t>
  event_base<strategy_t>::event_base(const event_base& src)
  : bs_refcounter (src), objbase (src)
  {
    *this = src;
  }

  template <typename strategy_t>
  void
  event_base<strategy_t>::init (const std::string &event_name, const std::string & params)
  {
#ifdef _DEBUG
    debug_data_.params = params;
#endif

    parse (params, main_params ());
    if (!no_error_)
      {
        BOSWARN (section::read_data, level::warning) <<
          boost::format ("Warnings exist while processing keyword %s") % event_name
          << bs_end;
      }
  }

  template <typename strategy_t>
  void
  event_base <strategy_t>::add_next_line (const std::string &event_name, const std::string &params)
  {
    const sp_named_pbase &next_line_params = add_next_line_params ();
    if (!next_line_params)
      {
        bs_throw_exception ("Can't create next_line_params");
      }

    parse (params, next_line_params);
    if (!no_error_)
      {
        BOSWARN (section::read_data, level::warning) <<
          boost::format ("Warnings exist while processing keyword %s") % event_name
          << bs_end;
      }
  }

  template <typename strategy_t>
  void
  event_base<strategy_t>::save_value(const sp_named_pbase &params_, char const* first, char const* last) //< save a string value
  {
    std::string str (first, last);

    for (int c = 0; c < num; ++c, ++index)
      {
        no_error_ = params_->set_value (index, str) && no_error_;
      }
  }


  template <typename strategy_t>
  void
  event_base<strategy_t>::save_def (const sp_named_pbase &params_, const char /*c*/) //< paste a default value
  {
    if (index + num > (int)params_->size ())
      {
        BOSERR (section::schedule, level::error) << "Event params parser: too many parameters (" << index << ", " << main_params ()->size () << ")" << bs_end;
        return ;
      } 

    index += num;
  }

  template <typename strategy_t>
  void
  event_base <strategy_t>::clear_num (const sp_named_pbase &, const char *, const char *)
  {
    num = 1;
  }

  template <typename strategy_t>
  void
  event_base<strategy_t>::parse(const std::string &params_str, const sp_named_pbase &params)
  {
    BS_ASSERT (params);

    detail::char_functor <this_t> save_def_   (this, params, &self_t::save_def);
    detail::str_functor  <this_t> clear_num_  (this, params, &self_t::clear_num);
    detail::str_functor  <this_t> save_value_ (this, params, &self_t::save_value);

    // matches to tabs and/or spaces repeated one or more times
    rule<> delimeter_p = +(space_p); 

    // matches to any type value
    rule<> value_p = confix_p ("\"", (+anychar_p)[save_value_], "\"")
                    | confix_p ("'", (+anychar_p)[save_value_], "'")
                    | (+graph_p)[save_value_];

    // matches to * (default value)
    rule<> default_value_p = ch_p ('*') [save_def_];

    // handled N*A constructions - repeat A value N times
    rule<> wrap_p = int_p[assign_a (num)] >> ch_p ('*') >> value_p;

    // handled N* constructions - repeat default value N times
    rule<> wrap_default_p = int_p[assign_a (num)] >> ch_p ('*') [save_def_];

    // matches to parameter
    rule <> param_p = wrap_p
                    | wrap_default_p
                    | (*anychar_p)[clear_num_] >> nothing_p
                    | default_value_p
                    | value_p
                    ; 

    rule <> event_p = param_p >> *(delimeter_p >> param_p); //main parser with delimiters

    index = 1;
    num = 1;
    no_error_ = true;
    boost::spirit::parse (params_str.c_str(), event_p).full; //Go! Go! Go!

    num = (int) params->size()-index;
    //if number of values in parameter string is
    //less then the number of prescribed parameters
    //we skip rest as defaults
    save_def (params);
  }

  template <typename strategy_t>
  void
  event_base<strategy_t>::apply (const sp_top &/*top*/, const sp_mesh_iface_t &/*mesh*/, const sp_calc_model_t &/*calc_model*/) const
  {
    BS_ASSERT (false && "BASE METHOD CALL");
  }

  //bs stuff
  BLUE_SKY_TYPE_STD_CREATE_T_DEF (event_base, (class))
  BLUE_SKY_TYPE_STD_COPY_T_DEF (event_base, (class))
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (event_base<base_strategy_fi>), 1, (property_base), "BOS_Core event_base class (seq fi)", "BOS_Core event_base class", "BOS_Core event_base class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (event_base<base_strategy_di>), 1, (property_base), "BOS_Core event_base class (seq di)", "BOS_Core event_base class", "BOS_Core event_base class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (event_base<base_strategy_mixi>), 1, (property_base), "BOS_Core event_base class (seq mixi)", "BOS_Core event_base class", "BOS_Core event_base class", false)


  bool
  event_base_register_types (const plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, event_base <base_strategy_fi>::bs_type ());
    BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, event_base <base_strategy_di>::bs_type ());
    BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, event_base <base_strategy_mixi>::bs_type ());
    BS_ASSERT (res);

    return true;
  }
}//ns bs

