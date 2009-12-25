/**
 *       \file  event_filter.h
 *      \brief  Filters events
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  15.12.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_EVENT_FILTER_H_
#define BS_EVENT_FILTER_H_

namespace blue_sky
  {

    /**
     * \class event_filter
     * \brief Filters events
     * */
  class BS_API_PLUGIN event_filter : public objbase
    {
    public:

      /**
       * \brief  event_filter dtor
       * */
      virtual ~event_filter () {}


      /**
       * \brief  Determines is well accepted
       * \param  name Name of well to test
       * \return True is well with name name accepted
       * */
      virtual bool 
      accept_well (const std::string &name) const;

      /**
       * \brief  Adds new well name that should be filtered
       * \param  well_name Name of well
       * */
      virtual void 
      add_filter_well (const std::string &well_name);

      /**
       * \brief  Sets reject_all flag
       * \param  reject_all New value of flag
       * */
      virtual void 
      set_reject_all (bool reject_all);

      //! blue-sky type declaration
      BLUE_SKY_TYPE_DECL (event_filter);

    protected:

      typedef std::list <std::string> well_name_list_t;   //!< Type for list of well names

      well_name_list_t            well_name_list_;        //!< List of well names
      auto_value <bool, false>    reject_all_;            //!< Flag to reject all wells

    };


} // namespace blue_sky


#endif  // BS_EVENT_FILTER_H_

