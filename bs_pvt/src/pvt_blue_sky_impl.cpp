/**
 * \file pvt_blue_sky_impl.cpp
 * \brief impl of necessary methods
 * \author Miryanov Sergey
 * \date 08.05.2008
 */

#include "pvt_base.h"
#include "pvt_dead_oil.h"
#include "pvt_gas.h"
#include "pvt_oil.h"
#include "pvt_water.h"
#include "pvt_dummy.h"
#include "pvt_3p.h"

namespace blue_sky
  {
#if 0
  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE_T_DEF (pvt_dead_oil,(class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF (pvt_dead_oil,(class));
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (pvt_dead_oil<base_strategy_fi>), 1, (pvt_base<base_strategy_fi>), "pvt_dead_oil_fi", "pvt_dead_oil_fi", "pvt_dead_oil_fi", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (pvt_dead_oil<base_strategy_di>), 1, (pvt_base<base_strategy_di>), "pvt_dead_oil_di", "pvt_dead_oil_di", "pvt_dead_oil_di", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (pvt_dead_oil<base_strategy_mixi>), 1, (pvt_base<base_strategy_mixi>), "pvt_dead_oil_mixi", "pvt_dead_oil_mixi", "pvt_dead_oil_mixi", false);

  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE_T_DEF (pvt_gas,(class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF (pvt_gas,(class));
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (pvt_gas<base_strategy_fi>), 1, (pvt_base<base_strategy_fi>), "pvt_gas_fi", "pvt_gas_fi", "pvt_gas_fi", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (pvt_gas<base_strategy_di>), 1, (pvt_base<base_strategy_di>), "pvt_gas_di", "pvt_gas_di", "pvt_gas_di", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (pvt_gas<base_strategy_mixi>), 1, (pvt_base<base_strategy_mixi>), "pvt_gas_mixi", "pvt_gas_mixi", "pvt_gas_mixi", false);

  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE_T_DEF (pvt_oil,(class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF (pvt_oil,(class));
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (pvt_oil<base_strategy_fi>), 1, (pvt_base<base_strategy_fi>), "pvt_oil_fi", "pvt_oil_fi", "pvt_oil_fi", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (pvt_oil<base_strategy_di>), 1, (pvt_base<base_strategy_di>), "pvt_oil_di", "pvt_oil_di", "pvt_oil_di", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (pvt_oil<base_strategy_mixi>), 1, (pvt_base<base_strategy_mixi>), "pvt_oil_mixi", "pvt_oil_mixi", "pvt_oil_mixi", false);

  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE_T_DEF (pvt_water,(class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF (pvt_water,(class));
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (pvt_water<base_strategy_fi>), 1, (pvt_base<base_strategy_fi>), "pvt_water_fi", "pvt_water_fi", "pvt_water_fi", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (pvt_water<base_strategy_di>), 1, (pvt_base<base_strategy_di>), "pvt_water_di", "pvt_water_di", "pvt_water_di", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (pvt_water<base_strategy_mixi>), 1, (pvt_base<base_strategy_mixi>), "pvt_water_mixi", "pvt_water_mixi", "pvt_water_mixi", false);
#endif 

  //////////////////////////////////////////////////////////////////////////
  bool
  pvt_register_types (const plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, pvt_dead_oil::bs_type ());
    BS_ASSERT (res);

    res &= BS_KERNEL.register_type (pd, pvt_gas::bs_type ());
    BS_ASSERT (res);

    res &= BS_KERNEL.register_type (pd, pvt_oil::bs_type ());
    BS_ASSERT (res);

    res &= BS_KERNEL.register_type (pd, pvt_water::bs_type ());
    BS_ASSERT (res);
    
    res &= BS_KERNEL.register_type (pd, pvt_dummy::bs_type ());
    BS_ASSERT (res);

    res &= BS_KERNEL.register_type (pd, pvt_3p::bs_type ());
    BS_ASSERT (res);
    
    return res;
  }


} // namespace blue_sky
