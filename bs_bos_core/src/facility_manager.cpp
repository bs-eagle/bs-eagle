/**
 *       \file  facility_manager.cpp
 *      \brief  Implementation of facility manager
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  17.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"

#include "facility_manager.h"
#include "data_storage_interface.h"
#include "calc_well.h"
#include "well_connection.h"
#include "reservoir.h"
#include "facility_manager.h"

namespace blue_sky
  {

  BLUE_SKY_TYPE_IMPL_T_EXT (2, (bos_val_table <std::string, facility_manager<base_strategy_fi>::sp_facility_t>), 1, (objbase), "bos_map<str, facility>_seq_fi", "map of str and facility", "map of str and facility", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (2, (bos_val_table <std::string, facility_manager<base_strategy_di>::sp_facility_t>), 1, (objbase), "bos_map<str, facility>_seq_di", "map of str and facility", "map of str and facility", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (2, (bos_val_table <std::string, facility_manager<base_strategy_mixi>::sp_facility_t>), 1, (objbase), "bos_map<str, facility>_seq_mixi", "map of str and facility", "map of str and facility", false);

  /**
   * \brief  'default' ctor for facility_manager
   * \param  param Additional params for ctor
   * */
  template <typename strategy_t>
  facility_manager<strategy_t>::facility_manager (bs_type_ctor_param /* param = NULL */)
      : facility_list_ (BS_KERNEL.create_object (facility_map_t::bs_type ()))
  {

  }

  /**
   * \brief  copy-ctor for facility_manager
   * \param  fs Instance of facility_manager to be copied
   * */
  template <typename strategy_t>
  facility_manager<strategy_t>::facility_manager (const facility_manager &fs)
  : bs_refcounter (fs), objbase (fs), facility_list_ (fs.facility_list_)
  {

  }

  template <typename strategy_t>
  typename facility_manager<strategy_t>::sp_well_t
  facility_manager<strategy_t>::get_well (const std::string &/*group_name*/, const std::string &well_name) const
    {
      typename facility_map_t::const_iterator it = facility_list_->find (well_name);
      if (it == facility_list_->end ())
        {
          return 0;
        }

      sp_well_t w (it->second, bs_dynamic_cast ());
      BS_ASSERT (w) (well_name);

      return w;
    }

  template <typename strategy_t>
  typename facility_manager<strategy_t>::sp_well_t
  facility_manager<strategy_t>::get_well (const std::string &well_name) const
    {
      typename facility_map_t::const_iterator it = facility_list_->find (well_name);
      if (it == facility_list_->end ())
        {
          BS_ASSERT (false && "Can't find well") (well_name);
          return 0;
        }

      sp_well_t w (it->second, bs_dynamic_cast ());
      return w;
    }

  template <typename strategy_t>
  void
  facility_manager<strategy_t>::add_well (const sp_well_t &w)
  {
    BS_ASSERT (w);
    BS_ASSERT (w->get_name ().size ());

    facility_list_->insert (std::make_pair (w->get_name (), w));
    BOSOUT (section::wells, level::low) << "Inserted " << w->get_name () << bs_end;
  }

  template <typename strategy_t>
  void
  facility_manager<strategy_t>::save_data (const sp_storage_t &storage) const
    {
      typename facility_map_t::const_iterator it = facility_list_->begin (), e = facility_list_->end ();
      for (; it != e; ++it)
        {
          storage->save (it->second);
        }
    }

  template <typename strategy_t>
  typename facility_manager<strategy_t>::well_const_iterator_t
  facility_manager<strategy_t>::wells_begin () const
    {
      return facility_list_->begin ();
      //const sp_facility_map_t &locked_list (facility_list_);
      //return well_iterator_t (locked_list, locked_list->begin ());
    }
  template <typename strategy_t>
  typename facility_manager<strategy_t>::well_const_iterator_t
  facility_manager<strategy_t>::wells_end () const
    {
      return facility_list_->end ();
      //const sp_facility_map_t &locked_list (facility_list_);
      //return well_iterator_t (locked_list, locked_list->end ());
    }
  //template <typename strategy_t>
  //typename facility_manager<strategy_t>::well_iterator_t
  //facility_manager<strategy_t>::wells_begin ()
  //  {
  //    return facility_list_->begin ();
  //    //const sp_facility_map_t &locked_list (facility_list_);
  //    //return well_iterator_t (locked_list, locked_list->begin ());
  //  }
  //template <typename strategy_t>
  //typename facility_manager<strategy_t>::well_iterator_t
  //facility_manager<strategy_t>::wells_end ()
  //  {
  //    return facility_list_->end ();
  //    //const sp_facility_map_t &locked_list (facility_list_);
  //    //return well_iterator_t (locked_list, locked_list->end ());
  //  }

  /*template <typename strategy_t>
  unsigned int
  facility_manager<strategy_t>::wells_size () const {
    return facility_list_.size ();
  }*/

  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE_T_DEF (facility_manager, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF (facility_manager, (class));
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (facility_manager<base_strategy_fi>), 1, (objbase), "facility_manager_seq_fi", "facility_manager_seq_fi", "facility_manager_seq_fi", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (facility_manager<base_strategy_di>), 1, (objbase), "facility_manager_seq_di", "facility_manager_seq_di", "facility_manager_seq_di", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (facility_manager<base_strategy_mixi>), 1, (objbase), "facility_manager_seq_mixi", "facility_manager_seq_mixi", "facility_manager_seq_mixi", false);
  //////////////////////////////////////////////////////////////////////////

  bool
  facility_manager_register_type (const blue_sky::plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, facility_manager<base_strategy_fi>::bs_type ());
    BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, facility_manager<base_strategy_di>::bs_type ());
    BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, facility_manager<base_strategy_mixi>::bs_type ());
    BS_ASSERT (res);

    return res;
  }

} // namespace blue_sky
