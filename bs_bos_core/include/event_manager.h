/**
* @file well_event.h
* @brief declaration of event manager class
* @author Morozov Andrey
* @date 2008-06-07
*/
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
  * @brief even_manager class to store events and to manage of them.
  */
  template <typename strategy_t>
  class BS_API_PLUGIN event_manager : public objbase
    {
      //-----------------------------------------
      //  TYPES
      //-----------------------------------------
    public:
      typedef event_base <strategy_t>                 event_base_t;
      typedef smart_ptr < event_base_t >              sp_event_base;
      typedef reservoir <strategy_t>                  sp_top;
      typedef event_manager<strategy_t>               self_t;
      typedef self_t                                  this_t;
      typedef std::list < sp_event_base >             sp_event_base_list;
      typedef boost::posix_time::ptime                date_t;
      typedef std::map <date_t, sp_event_base_list >  event_map;

      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
      //blue-sky class declaration
    public:
      BLUE_SKY_TYPE_DECL(event_manager);

    public:
      //! destructor
      virtual ~event_manager ();

      void 
      process_event (const date_t &date, const std::string &event_name, const std::string &event_params);

      void
      end_event ();

      sp_event_base 
      create_event(const boost::posix_time::ptime &date, const std::string & event_name, const std::string & event_params);


      //void save_keyword(iterator_t first, iterator_t const& last);
      //void free_keyword(iterator_t first, iterator_t const& last);
      //void read_event(iterator_t first, iterator_t const& last);

      //-----------------------------------------
      //  VARIABLES
      //-----------------------------------------
    public:
      event_map         event_list;

    private:
      sp_event_base     current_event_;
    };

  bool
  well_events_register_type (const blue_sky::plugin_descriptor &pd);

  bool
  event_manager_register_types (const plugin_descriptor &pd);

}//ns bs

#endif  // #ifndef EVENT_MANAGER_H_
