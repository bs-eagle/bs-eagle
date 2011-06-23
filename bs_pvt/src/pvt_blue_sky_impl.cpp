/**
 * \file pvt_blue_sky_impl.cpp
 * \brief impl of necessary methods
 * \author Miryanov Sergey
 * \date 08.05.2008
 */
#include "bs_pvt_stdafx.h"

#include "pvt_base.h"
#include "pvt_dead_oil.h"
#include "pvt_gas.h"
#include "pvt_oil.h"
#include "pvt_water.h"

namespace blue_sky
  {

  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE (pvt_dead_oil);
  BLUE_SKY_TYPE_STD_COPY (pvt_dead_oil);
  BLUE_SKY_TYPE_IMPL (pvt_dead_oil, pvt_base, "pvt_dead_oil", "pvt_dead_oil", "pvt_dead_oil");

  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE (pvt_gas);
  BLUE_SKY_TYPE_STD_COPY (pvt_gas);
  BLUE_SKY_TYPE_IMPL (pvt_gas, pvt_base, "pvt_gas", "pvt_gas", "pvt_gas");

  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE (pvt_oil);
  BLUE_SKY_TYPE_STD_COPY (pvt_oil);
  BLUE_SKY_TYPE_IMPL (pvt_oil, pvt_base, "pvt_oil", "pvt_oil", "pvt_oil");

  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE (pvt_water);
  BLUE_SKY_TYPE_STD_COPY (pvt_water);
  BLUE_SKY_TYPE_IMPL (pvt_water, pvt_base, "pvt_water", "pvt_water", "pvt_water");

  //////////////////////////////////////////////////////////////////////////
  bool
  pvt_register_types (const plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, pvt_dead_oil::bs_type ());  BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, pvt_gas::bs_type ());       BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, pvt_oil::bs_type ());       BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, pvt_water::bs_type ());     BS_ASSERT (res);

    return res;
  }


} // namespace blue_sky
