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
  event_manager::process_event (const date_t &date, const std::string &event_name, const std::string &event_params)
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
  event_manager::end_event ()
  {
    current_event_ = 0;
  }

  event_manager::sp_event_base
  event_manager::create_event (const boost::posix_time::ptime &date, const std::string & event_name, const std::string & event_params)
  {
    sp_obj event_object = BS_KERNEL.create_object (event_name);
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


#define REG_EVENT_IN_KERNEL(NAME)\
  res &= BS_KERNEL.register_type (pd, NAME##_event::bs_type ()); BS_ASSERT (res);

  //bs stuff
  BLUE_SKY_TYPE_STD_CREATE (event_manager)
  BLUE_SKY_TYPE_STD_COPY (event_manager)
  BLUE_SKY_TYPE_IMPL (event_manager, objbase, "event_manager", "event_manager", "event_manager")

  bool
  well_events_register_type (const blue_sky::plugin_descriptor &pd)
  {
    bool res  = true;
    REG_EVENT_IN_KERNEL(WELSPECS)
    REG_EVENT_IN_KERNEL(WELLCON)
    REG_EVENT_IN_KERNEL(COMPDAT)
    REG_EVENT_IN_KERNEL(WCONPROD)
    REG_EVENT_IN_KERNEL(WCONHIST)
    REG_EVENT_IN_KERNEL(WCONINJE)
    REG_EVENT_IN_KERNEL(WECON)
    REG_EVENT_IN_KERNEL(WECONINJ)
    REG_EVENT_IN_KERNEL(WEFAC)
    REG_EVENT_IN_KERNEL(WELTARG)
    REG_EVENT_IN_KERNEL(WPIMULT)
    REG_EVENT_IN_KERNEL(COMPENSATION)
    REG_EVENT_IN_KERNEL(PERMFRAC)
    return res;
  }

  bool
  event_manager_register_types (const plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, event_manager::bs_type ()); BS_ASSERT (res);

    return res;
  }

}//ns bs

