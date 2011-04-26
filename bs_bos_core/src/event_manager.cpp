/**
 *       \file  event_manager.cpp
 *      \brief  event_manager functions
 *     \author  Morozov Andrey
 *       \date  07.06.2008
 *  \copyright  This source code is released under the terms of
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"

#include "event_manager.h"
#include "well_events.h"

#include <boost/regex.hpp>

//using namespace boost::spirit;
//using namespace boost::gregorian;

namespace blue_sky
  {

  event_manager::~event_manager ()
  {

  }

  /**
   * \brief  'default' ctor for event_manager
   * \param  param additional params for event_manager
   * */
  event_manager::event_manager(bs_type_ctor_param /*param*/)
  {
  }

  /**
   * \brief  copy-ctor for event_manager
   * \param  src source copy of event_manager
   * */
  event_manager::event_manager(const event_manager& src)
  : bs_refcounter (src), objbase (src)
  {
    *this = src;
  }

  void
  event_manager ::process_event (const date_t &date, const std::string &event_name, const std::string &event_params)
  {
    if (current_event_)
      {
        current_event_->add_next_line (event_name, event_params);
      }
    else
      {
        current_event_ = create_event (date, event_name, event_params);
        if (!current_event_->is_multi_line ())
          {
            current_event_ = 0;
          }
      }
  }
  void
  event_manager ::end_event ()
  {
    current_event_ = 0;
  }

  event_manager::sp_event_base
  event_manager::create_event (const boost::posix_time::ptime &date, const std::string & event_name, const std::string & event_params)
  {
    //boost::regex re_check_type (event_name + "(.*)_" + tools::strategy_name ::name ());

    sp_obj event_object;

    const std::vector <type_tuple> &types = BS_KERNEL.registered_types ();
    for (size_t i = 0, cnt = types.size (); i < cnt; ++i)
      {
        const type_descriptor &td = types[i].td_;

        // FIXME:
        //if (boost::regex_match (td.stype_.c_str (), re_check_type))
        if (td.stype_.c_str () == event_name)
          {
            event_object = BS_KERNEL.create_object (td);
            break;
          }
      }
    if (!event_object)
      {
        bs_throw_exception (boost::format ("Type (%s) not registered (params: %s)") % event_name % event_params);
      }

    sp_event_base event (event_object, bs_static_cast ());
    if (!event)
      {
        bs_throw_exception (boost::format ("Created object for type (%s) can't be casted to event_base (params: %s)") % event_name % event_params);
      }

    event_list[date].push_back (event);  // TODO: posix_time::ptime
    event->init (event_name, event_params);
    return event;
  }


  //bs stuff
  BLUE_SKY_TYPE_STD_CREATE (event_manager)
  BLUE_SKY_TYPE_STD_COPY (event_manager)
  BLUE_SKY_TYPE_IMPL (event_manager, objbase, "event_manager_seq", "BOS_Core event_manager class", "BOS_Core event_manager class")

  bool
  well_events_register_type (const blue_sky::plugin_descriptor &pd)
  {
    bool res  = true;
    res &= BS_KERNEL.register_type (pd, WELSPECS_event::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, WELLCON_event::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, COMPDAT_event::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, WCONPROD_event::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, WCONHIST_event::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, WCONINJE_event::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, WECON_event::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, WECONINJ_event::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, WEFAC_event::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, WELTARG_event::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, WPIMULT_event::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, COMPENSATION_event::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, PERMFRAC_event::bs_type ()); BS_ASSERT (res);
    return res;
  }

  bool
  event_manager_register_types (const plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, event_manager::bs_type ()); BS_ASSERT (res);

    return true;
  }

}//ns bs

