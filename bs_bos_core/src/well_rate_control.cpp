/**
 * \file well_rate_contorl.cpp
 * \brief well_rate_control impl
 * \author Sergey Miryanov
 * \date 21.11.2008
 * */
#include "stdafx.h"
#include "well_rate_control.h"

namespace blue_sky {
namespace wells {

  template <typename strategy_t>
  well_rate_control <strategy_t>::well_rate_control (bs_type_ctor_param /*param = NULL */)
  {

  }
  template <typename strategy_t>
  well_rate_control <strategy_t>::well_rate_control (const well_rate_control &x)
  : bs_refcounter (x), objbase (x)
  {

  }

  BLUE_SKY_TYPE_STD_CREATE_T_DEF (well_rate_control, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF (well_rate_control, (class));
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (well_rate_control <base_strategy_fi>), 1, (objbase), "well_rate_control_seq_fi", "Base class for well_rate_control", "Base class for well_rate_control", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (well_rate_control <base_strategy_di>), 1, (objbase), "well_rate_control_seq_di", "Base class for well_rate_control", "Base class for well_rate_control", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (well_rate_control <base_strategy_mixi>), 1, (objbase), "well_rate_control_seq_mixi", "Base class for well_rate_control", "Base class for well_rate_control", false);


} // namespace wells
} // namespace blue_sky

