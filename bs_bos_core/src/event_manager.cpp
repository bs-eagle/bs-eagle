/**
 * @file well_event.h
 * @brief event manager functions
 * @author Morozov Andrey
 * @date 2008-06-07
 */
#include "stdafx.h"

#include "event_manager.h"
#include "well_events.h"

#include <boost/regex.hpp>
#include "strategy_name.h"

//using namespace boost::spirit;
//using namespace boost::gregorian;

namespace blue_sky
  {

  template <typename strategy_t>
  event_manager<strategy_t>::~event_manager ()
  {

  }

  template <typename strategy_t>
  event_manager<strategy_t>::event_manager(bs_type_ctor_param /*param*/)
  {
  }

  template <typename strategy_t>
  event_manager<strategy_t>::event_manager(const event_manager& src)
  : bs_refcounter (src), objbase (src)
  {
    *this = src;
  }

  template <typename strategy_t>
  void
  event_manager <strategy_t>::process_event (const date_t &date, const std::string &event_name, const std::string &event_params)
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
  template <typename strategy_t>
  void
  event_manager <strategy_t>::end_event ()
  {
    current_event_ = 0;
  }

  template <typename strategy_t>
  typename event_manager<strategy_t>::sp_event_base
  event_manager<strategy_t>::create_event (const boost::posix_time::ptime &date, const std::string & event_name, const std::string & event_params)
  {
    boost::regex re_check_type (event_name + "(.*)_" + tools::strategy_name <strategy_t>::name ());

    sp_obj event_object;

    const std::vector <type_tuple> &types = BS_KERNEL.registered_types ();
    for (size_t i = 0, cnt = types.size (); i < cnt; ++i)
      {
        const type_descriptor &td = types[i].td_;
        
        if (boost::regex_match (td.stype_.c_str (), re_check_type))
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


#define EVENT_LIST_REGISTRATION(WHERE, strategy_t)	\
  REG_EVENT_IN_##WHERE(WELSPECS,			strategy_t)		\
  REG_EVENT_IN_##WHERE(WELLCON,				strategy_t)		\
  REG_EVENT_IN_##WHERE(COMPDAT,				strategy_t)		\
  REG_EVENT_IN_##WHERE(WCONPROD,			strategy_t)		\
  REG_EVENT_IN_##WHERE(WCONHIST,			strategy_t)		\
  REG_EVENT_IN_##WHERE(WCONINJE,			strategy_t)		\
  REG_EVENT_IN_##WHERE(WECON,					strategy_t)		\
  REG_EVENT_IN_##WHERE(WECONINJ,			strategy_t)		\
  REG_EVENT_IN_##WHERE(WEFAC,					strategy_t)		\
  REG_EVENT_IN_##WHERE(WELTARG,				strategy_t)		\
  REG_EVENT_IN_##WHERE(WPIMULT,				strategy_t)		\
  REG_EVENT_IN_##WHERE(COMPENSATION,	strategy_t)		\
  REG_EVENT_IN_##WHERE(FRACTURE,			strategy_t)		\
  REG_EVENT_IN_##WHERE(PERMFRAC,			strategy_t)

#define REG_EVENT_IN_FACTORY(NAME, strategy_t)\
  factory_register <NAME##_event <strategy_t> >     (#NAME);

#define REG_EVENT_IN_KERNEL(NAME, strategy_t)\
  res &= BS_KERNEL.register_type (pd, NAME##_event<strategy_t>::bs_type ()); BS_ASSERT (res);

  //bs stuff
  BLUE_SKY_TYPE_STD_CREATE_T_DEF (event_manager, (class))
  BLUE_SKY_TYPE_STD_COPY_T_DEF (event_manager, (class))
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (event_manager <base_strategy_fi>), 1, (objbase), "event_manager_seq_fi", "BOS_Core event_manager class", "BOS_Core event_manager class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (event_manager <base_strategy_di>), 1, (objbase), "event_manager_seq_di", "BOS_Core event_manager class", "BOS_Core event_manager class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (event_manager <base_strategy_mixi>), 1, (objbase), "event_manager_seq_mixi", "BOS_Core event_manager class", "BOS_Core event_manager class", false)

  bool
  well_events_register_type (const blue_sky::plugin_descriptor &pd)
  {
    bool res  = true;
    EVENT_LIST_REGISTRATION(KERNEL, base_strategy_fi);
    EVENT_LIST_REGISTRATION(KERNEL, base_strategy_di);
    EVENT_LIST_REGISTRATION(KERNEL, base_strategy_mixi);
    return res;
  }

  bool
  event_manager_register_types (const plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, event_manager <base_strategy_fi>::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, event_manager <base_strategy_di>::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, event_manager <base_strategy_mixi>::bs_type ()); BS_ASSERT (res);

    return true;
  }

}//ns bs

