/**
 *       \file  event_manager.h
 *      \brief  declaration of event_manager class
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  07.06.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef EVENT_MANAGER_H_
#define EVENT_MANAGER_H_

#include "event_base.h"

namespace blue_sky
  {

  template <typename strategy_t>
  class event_base;

  template <typename strategy_t>
  class reservoir;

  /**
   * \class event_manager
   * \brief class to store events and to manage of them
   * */
  template <typename strategy_t>
  class BS_API_PLUGIN event_manager : public objbase
    {
      //-----------------------------------------
      //  TYPES
      //-----------------------------------------
    public:
      typedef event_base <strategy_t>                 event_base_t;       //!< event_base type
      typedef smart_ptr < event_base_t >              sp_event_base;      //!< smart_ptr to event_base type
      typedef reservoir <strategy_t>                  sp_top;             //!< reservoir type
      typedef event_manager<strategy_t>               self_t;             //!< shortname for this type
      typedef self_t                                  this_t;             //!< shortname for this type
      typedef std::list < sp_event_base >             sp_event_base_list; //!< list of events
      typedef boost::posix_time::ptime                date_t;             //!< shortname for time
      typedef std::map <date_t, sp_event_base_list >  event_map;          //!< events map by date

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


      //-----------------------------------------
      //  VARIABLES
      //-----------------------------------------
    public:
      event_map         event_list;     //!< list of scheduled events

    private:
      sp_event_base     current_event_; //!< current processed event
    };

  /**
   * \brief  registers types of events in blue-sky kernel
   * \param  pd plugin_descriptor
   * \return true if all types registered successfully
   * */
  bool
  well_events_register_type (const blue_sky::plugin_descriptor &pd);

  /**
   * \brief  registers types of event_manager in blue-sky kernel
   * \param  pd plugin_descriptor
   * \return true if all types registered successfully
   * */
  bool
  event_manager_register_types (const plugin_descriptor &pd);

}//ns bs

#endif  // #ifndef EVENT_MANAGER_H_
