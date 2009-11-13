/**
 *       \file  facility_manager.h
 *      \brief  Facilities manager
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  16.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_FACILITY_MANAGER_H_
#define BS_FACILITY_MANAGER_H_

#include "facility_base.h"

namespace blue_sky
  {

  template <typename strategy_t>
  class well;

  class data_storage_interface;

  /**
   * \class facility_manager
   * \brief Facilities manager
   * */
  template <typename strategy_t>
  class BS_API_PLUGIN facility_manager : public objbase
    {
    public:

      typedef facility_base <strategy_t>                    facility_base_t;
      typedef well <strategy_t>                             well_t;

      typedef smart_ptr <facility_base_t>                   sp_facility_t;
      typedef bos_val_table <std::string, sp_facility_t>    facility_map_t;
      typedef smart_ptr <facility_map_t>                    sp_facility_map_t;
      typedef smart_ptr <well_t, true>                      sp_well_t;

      typedef smart_ptr <data_storage_interface, true>      sp_storage_t;

      typedef typename facility_map_t::iterator             facility_iterator_t;
      typedef typename facility_map_t::const_iterator       facility_const_iterator_t;

      //typedef well_iterator <locked_facility_map_t, facility_iterator_t, locked_well_t>  well_iterator_t;
      typedef facility_iterator_t                           well_iterator_t;
      typedef facility_const_iterator_t                     well_const_iterator_t;

    public:

      /**
       * \brief  Returns well with name well_name which belongs to group group_name
       * \param  group_name Name of group
       * \param  well_name Name of well
       * \return Instance of well or null pointer if not exists
       * */
      sp_well_t 
      get_well (const std::string &group_name, const std::string &well_name) const;

      /**
       * \brief  Returns first well with name well_name
       * \param  well_name Name of well
       * \return Instance of well or null pointer if not exists
       * */
      sp_well_t 
      get_well (const std::string &well_name) const;

      /**
       * \brief  Adds well
       * \param  well Instance of well to be added
       * */
      void 
      add_well (const sp_well_t &well);

      /**
       * \brief  Saves facilities data to storage
       * \param  storage
       * */
      void 
      save_data (const sp_storage_t &storage) const;

      /**
       * \brief  Returns begin iterator
       * \return Begin iterator
       * */
      well_const_iterator_t 
      wells_begin () const;

      /**
       * \brief  Returns end iterator
       * \return End iterator
       * */
      well_const_iterator_t 
      wells_end () const;

      //! blue-sky type declaration
      BLUE_SKY_TYPE_DECL_T (facility_manager);

    private:

      sp_facility_map_t         facility_list_;   //!< List of facilities

    };

  /**
   * \brief  Registers facility_manager types in blue-sky kernel
   * \param  pd plugin_descriptor
   * \return True if all types registered successfully
   * */
  bool
  facility_manager_register_type (const blue_sky::plugin_descriptor &pd);


} // namespace blue_sky


#endif  // #ifndef BS_FACILITY_MANAGER_H_

