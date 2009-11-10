/**
 * \file event_filter.h
 * \brief filter events
 * \author Sergey Miryanov
 * \date 15.12.2008
 * */
#ifndef BS_EVENT_FILTER_H_
#define BS_EVENT_FILTER_H_

namespace blue_sky
  {

  class BS_API_PLUGIN event_filter : public objbase
    {
    public:

      virtual ~event_filter () {}
      virtual bool accept_well (const std::string &name) const;

      virtual void add_filter_well (const std::string &well_name);
      virtual void set_reject_all (bool reject_all);

      BLUE_SKY_TYPE_DECL (event_filter);

    protected:

      typedef std::list <std::string> well_name_list_t;

      well_name_list_t            well_name_list_;
      auto_value <bool, false>    reject_all_;

    };


} // namespace blue_sky


#endif  // BS_EVENT_FILTER_H_

