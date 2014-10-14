/**
 *       \file  event_manager.cpp
 *      \brief  event_manager functions
 *     \author  Morozov Andrey
 *       \date  07.06.2008
 *  \copyright  This source code is released under the terms of
 *              the BSD License. See LICENSE for more details.
 * */
#include "event_manager.h"

#include <boost/regex.hpp>

//using namespace boost::spirit;
//using namespace boost::gregorian;

namespace pt = boost::posix_time;

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
  : bs_refcounter (src), event_manager_iface (src)
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
  event_manager::create_event (const double date, const std::string & event_name, const std::string & event_params)
  {
    BS_SP (event_base) event = BS_KERNEL.create_object (event_name);
    if (!event)
      {
        bs_throw_exception (boost::format ("Type (%s) not registered (params: %s)") % event_name % event_params);
      }

    event->init (event_name, event_params);
    event_map::iterator it = event_list.find (date);
    event_list_t &e = it->second;

    e.push_back (event);
    return event;
  }

  void
  event_manager::set_current_date (date_t date)
  {
    if (event_list.find (date) == event_list.end ())
      {
        event_list.insert (std::make_pair (date, event_list_t ()));
        current_date_ = date;
      }
  }

  event_manager::date_t
  event_manager::get_current_date () const
  {
    return current_date_;
  }

  void
  event_manager::finalize_events ()
  {
    if (event_list.empty ())
      {
        bs_throw_exception ("Events list is empty");
      }

    event_map::reverse_iterator it = event_list.rbegin ();
    if (it->second.size ())
      {
        event_list.insert (std::make_pair (it->first + 30.0, event_list_t ()));
      }
  }

  //bs stuff
  BLUE_SKY_TYPE_STD_CREATE (event_manager)
  BLUE_SKY_TYPE_STD_COPY (event_manager)
  BLUE_SKY_TYPE_IMPL (event_manager, event_manager_iface, "event_manager", "BOS_Core event_manager class", "BOS_Core event_manager class")

}//ns bs

