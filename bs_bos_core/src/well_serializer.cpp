/**
 *       \file  well_serializer.cpp
 *      \brief  Implementation of well_serializer
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  21.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Obsolete
 * */

#include "calc_well.h"
#include "well_serializer.h"
#include "well_connection.h"
#include "reservoir.h"
#include "facility_manager.h"
namespace blue_sky
  {
  namespace wells
    {

    well_serializer::well_serializer (bs_type_ctor_param /*param = NULL */)
    {
      handled_type_ = well::bs_type ();
    }
    well_serializer::well_serializer (const well_serializer& w)
    : bs_refcounter (w), data_serializer (w)
    {
      handled_type_ = well::bs_type ();
    }

    void
    well_serializer::save (const sp_storage_t &storage, const sp_obj &obj) const
      {
        BS_SP (well) w (obj, bs_dynamic_cast ());

        BS_ASSERT (w) (obj->bs_resolve_type ().stype_);

        const sp_storage_t &locked_storage (storage);
        data_storage &st = *locked_storage;

        st.save ("i_coord", w->i_coord_)
        .save ("j_coord", w->j_coord_)
        .save ("name",    w->name_)
        ;
      }

    //////////////////////////////////////////////////////////////////////////
    BLUE_SKY_TYPE_STD_CREATE (well_serializer);
    BLUE_SKY_TYPE_STD_COPY (well_serializer);
    BLUE_SKY_TYPE_IMPL (well_serializer, data_serializer, "well_serializer_seq", "well_serializer", "well_serializer");

    //////////////////////////////////////////////////////////////////////////
    bool
    well_serializer_register_type (const plugin_descriptor &pd)
    {
      bool res = true;

      res &= BS_KERNEL.register_type (pd, well_serializer::bs_type ()); BS_ASSERT (res);

      return res;
    }
  } // namespace wells
} // namespace blue_sky
