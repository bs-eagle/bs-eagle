/**
 *       \file  facility_manager.cpp
 *      \brief  Implementation of facility manager
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  17.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "facility_manager.h"
#include "data_storage_interface.h"
#include "calc_well.h"
#include "well_connection.h"
#include "reservoir.h"
#include "facility_manager.h"

namespace blue_sky
  {

  BLUE_SKY_TYPE_IMPL_T_EXT (2, (bos_val_table <std::string, facility_manager::sp_facility_t>), 1, (objbase), "bos_map<str, facility>_seq", "map of str and facility", "map of str and facility", false);

  /**
   * \brief  'default' ctor for facility_manager
   * \param  param Additional params for ctor
   * */
  facility_manager::facility_manager (bs_type_ctor_param /* param = NULL */)
      : facility_list_ (BS_KERNEL.create_object (facility_map_t::bs_type ()))
  {

  }

  /**
   * \brief  copy-ctor for facility_manager
   * \param  fs Instance of facility_manager to be copied
   * */
  facility_manager::facility_manager (const facility_manager &fs)
  : bs_refcounter (fs), objbase (fs), facility_list_ (fs.facility_list_)
  {

  }

  facility_manager::sp_well_t
  facility_manager::get_well (const std::string &/*group_name*/, const std::string &well_name) const
    {
      facility_map_t::const_iterator it = facility_list_->find (well_name);
      if (it == facility_list_->end ())
        {
          return 0;
        }

      sp_well_t w (it->second, bs_dynamic_cast ());
      BS_ASSERT (w) (well_name);

      return w;
    }

  facility_manager::sp_well_t
  facility_manager::get_well (const std::string &well_name) const
    {
      facility_map_t::const_iterator it = facility_list_->find (well_name);
      if (it == facility_list_->end ())
        {
          BS_ASSERT (false && "Can't find well") (well_name);
          return 0;
        }

      sp_well_t w (it->second, bs_dynamic_cast ());
      return w;
    }

  void
  facility_manager::add_well (const sp_well_t &w)
  {
    BS_ASSERT (w);
    BS_ASSERT (w->get_name ().size ());

    facility_list_->insert (std::make_pair (w->get_name (), w));
    BOSOUT (section::wells, level::low) << "Inserted " << w->get_name () << bs_end;
  }

  void
  facility_manager::save_data (const sp_storage_t &storage) const
    {
      facility_map_t::const_iterator it = facility_list_->begin (), e = facility_list_->end ();
      for (; it != e; ++it)
        {
          storage->save (it->second);
        }
    }

  facility_manager::well_const_iterator_t
  facility_manager::wells_begin () const
    {
      return facility_list_->begin ();
      //const sp_facility_map_t &locked_list (facility_list_);
      //return well_iterator_t (locked_list, locked_list->begin ());
    }
  facility_manager::well_const_iterator_t
  facility_manager::wells_end () const
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
  BLUE_SKY_TYPE_STD_CREATE (facility_manager);
  BLUE_SKY_TYPE_STD_COPY (facility_manager);
  BLUE_SKY_TYPE_IMPL (facility_manager, objbase, "facility_manager_seq", "facility_manager_seq", "facility_manager_seq");
  //////////////////////////////////////////////////////////////////////////

  bool
  facility_manager_register_type (const blue_sky::plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, facility_manager::bs_type ()); BS_ASSERT (res);

    return res;
  }

} // namespace blue_sky
