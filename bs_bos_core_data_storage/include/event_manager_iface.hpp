#ifndef EVENT_MANAGER_IFACE_HPP_5f2c21e8_7711_11e0_97e3_638dd1a9b36a
#define EVENT_MANAGER_IFACE_HPP_5f2c21e8_7711_11e0_97e3_638dd1a9b36a
/**
 *       \file  event_manager_iface.hpp
 *      \brief  interface class for events manager
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  05.05.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

namespace blue_sky 
{
  class BS_API_PLUGIN event_manager_iface : public objbase
  {
  public:
    typedef boost::posix_time::ptime date_t;             //!< shortname for time

    /**
     * \brief  ends processing of event
     * \details if processing of event is not ended
     *          all folowing events will be treated as
     *          subevents for current event
     * */
    virtual void
    end_event () = 0;

    /**
     * \brief  creates new event or adds event to previous one
     *         and parses event params
     * \param  date date of the event
     * \param  event_name name of event
     * \param  event_params string with event params
     * \return may throw exception
     * */
    virtual void 
    process_event (const date_t &date, const std::string &event_name, const std::string &event_params) = 0;

    /**
     * \brief adds empty event list at the end of date list if last date has events
     * \return throws exception if date list is empty
     * */
    virtual void
    finalize_events () = 0;

    /**
     * \brief set current date (usualy parsed by DATE, TSTEP keywords)
     * \param date 
     * */
    virtual void
    set_current_date (date_t const &date) = 0;

    /**
     * \brief returns current date
     * */
    virtual date_t const &
    get_current_date () const = 0;
  };

}


#endif

