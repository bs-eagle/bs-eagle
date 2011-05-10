/**
 *       \file  event_manager.h
 *      \brief  declaration of event_manager class
 *     \author  Morozov Andrey
 *       \date  07.06.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef EVENT_MANAGER_H_
#define EVENT_MANAGER_H_

#include "event_base.h"
#include "event_manager_iface.hpp"

namespace blue_sky
  {
  class event_base;

  /**
   * \class event_manager
   * \brief class to store events and to manage of them
   * */
  class BS_API_PLUGIN event_manager : public event_manager_iface
    {
      //-----------------------------------------
      //  TYPES
      //-----------------------------------------
    public:
      typedef event_base                              event_base_t;       //!< event_base type
      typedef smart_ptr < event_base_t >              sp_event_base;      //!< smart_ptr to event_base type
      typedef event_manager                           self_t;             //!< shortname for this type
      typedef self_t                                  this_t;             //!< shortname for this type
      typedef std::list < sp_event_base >             event_list_t; //!< list of events
      typedef boost::posix_time::ptime                date_t;             //!< shortname for time
      typedef std::map <date_t, event_list_t>         event_map;          //!< events map by date

      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
      //blue-sky class declaration
    public:
      /**
       * \brief  blue-sky type declaration
       * */
      BLUE_SKY_TYPE_DECL(event_manager);

    public:
      /**
       * \brief  event_manager dtor
       * */
      virtual ~event_manager ();

      /**
       * \brief  creates new event or adds event to previous one
       *         and parses event params
       * \param  date date of the event
       * \param  event_name name of event
       * \param  event_params string with event params
       * \return may throw exception
       * */
      void 
      process_event (const date_t &date, const std::string &event_name, const std::string &event_params);

      /**
       * \brief  ends processing of event
       * \details if processing of event is not ended
       *          all folowing events will be treated as
       *          subevents for current event
       * */
      void
      end_event ();

      void
      set_current_date (date_t const &date);

      date_t const &
      get_current_date () const;

      /**
       * \brief  creats new event or adds event to previous one
       *         and parses event params
       * \param  date date of the event
       * \param  event_name name of event
       * \param  event_params string with event params
       * \return may throw exception
       * */
      sp_event_base 
      create_event(const boost::posix_time::ptime &date, const std::string & event_name, const std::string & event_params);


      /**
       * \brief adds empty event list at the end of date list if last date has events
       * \return throws exception if date list is empty
       * */
      virtual void
      finalize_events ();
      //-----------------------------------------
      //  VARIABLES
      //-----------------------------------------
    public:
      event_map         event_list;     //!< list of scheduled events

    private:
      sp_event_base     current_event_; //!< current processed event
      date_t            current_date_;
    };

}//ns bs

#endif  // #ifndef EVENT_MANAGER_H_
