/**
 *       \file  well_serializer.cpp
 *      \brief  Implementation of well_serializer
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  21.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Obsolete
 * */
#include "stdafx.h"

#include "calc_well.h"
#include "well_serializer.h"
#include "well_connection.h"
#include "reservoir.h"
#include "facility_manager.h"
namespace blue_sky
  {
  namespace wells
    {

    template <typename strategy_t>
    well_serializer<strategy_t>::well_serializer (bs_type_ctor_param param /* = NULL */)
    {
      handled_type_ = well_t::bs_type ();
    }
    template <typename strategy_t>
    well_serializer<strategy_t>::well_serializer (const well_serializer& w)
    : bs_refcounter (w), data_serializer (w)
    {
      handled_type_ = well_t::bs_type ();
    }

    template <typename strategy_t>
    void
    well_serializer<strategy_t>::save (const sp_storage_t &storage, const sp_obj &obj) const
      {
        sp_well_t w (obj, bs_dynamic_cast ());

        BS_ASSERT (w) (obj->bs_resolve_type ().stype_);

        const sp_storage_t &locked_storage (storage);
        data_storage &st = *locked_storage;

        st.save ("i_coord", w->i_coord_)
        .save ("j_coord", w->j_coord_)
        .save ("name",    w->name_)
        ;
      }

    //////////////////////////////////////////////////////////////////////////
    BLUE_SKY_TYPE_STD_CREATE_T_DEF (well_serializer, (class));
    BLUE_SKY_TYPE_STD_COPY_T_DEF (well_serializer, (class));
    BLUE_SKY_TYPE_IMPL_T_EXT (1, (well_serializer <base_strategy_fi>), 1, (data_serializer), "well_serializer_seq_fi", "well_serializer", "well_serializer", false);
    BLUE_SKY_TYPE_IMPL_T_EXT (1, (well_serializer <base_strategy_di>), 1, (data_serializer), "well_serializer_seq_di", "well_serializer", "well_serializer", false);

    //////////////////////////////////////////////////////////////////////////
    bool
    well_serializer_register_type (const plugin_descriptor &pd)
    {
      bool res = true;

      res &= BS_KERNEL.register_type (pd, well_serializer<base_strategy_fi>::bs_type ());
      BS_ASSERT (res);
      res &= BS_KERNEL.register_type (pd, well_serializer<base_strategy_di>::bs_type ());
      BS_ASSERT (res);

      return res;
    }
  } // namespace wells
} // namespace blue_sky
