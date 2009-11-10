/**
 * \file
 * \brief
 * \author
 * \date
 * */
#include "stdafx.h"

#include "calc_well.h"
#include "well_connection.h"
#include "well_controller.h"
#include "well_rate_control.h"

#include "calc_model_data_accessors.h"

#include "rate_control_type.h"

#include "well_rate_control_prod_mobility.h"
#include "well_rate_control_inj_mobility.h"

#include "well_rate_control_bhp_deriv.h"
#include "well_rate_control_rate_deriv.h"
#include "well_rate_control_rate.h"

#include "well_rate_control_impl.h"
#include "well_rate_control_impl_type.h"

#include "reservoir.h"
#include "facility_manager.h"

namespace blue_sky
  {


  template <typename impl_type_t>
  well_rate_control_impl <impl_type_t>::well_rate_control_impl (bs_type_ctor_param param /* = NULL */)
  {
  }

  template <typename impl_type_t>
  well_rate_control_impl <impl_type_t>::well_rate_control_impl (const well_rate_control_impl &i)
  : bs_refcounter (i)
  {
  }

  BLUE_SKY_TYPE_STD_CREATE_T_DEF (well_rate_control_impl, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF (well_rate_control_impl, (class));

#define DEFINE_WELL_RATE_CONTROL_IMPL(TYPE) \
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (well_rate_control_impl <TYPE <base_strategy_fi> >), 	1, (objbase), BOOST_PP_STRINGIZE(well_rate_control_impl <TYPE <base_strategy_fi> >), BOOST_PP_STRINGIZE(well_rate_control_impl <TYPE <base_strategy_fi> >), BOOST_PP_STRINGIZE(well_rate_control_impl <TYPE <base_strategy_fi> >), false);  \
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (well_rate_control_impl <TYPE <base_strategy_di> >), 	1, (objbase), BOOST_PP_STRINGIZE(well_rate_control_impl <TYPE <base_strategy_di> >), BOOST_PP_STRINGIZE(well_rate_control_impl <TYPE <base_strategy_di> >), BOOST_PP_STRINGIZE(well_rate_control_impl <TYPE <base_strategy_di> >), false);\
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (well_rate_control_impl <TYPE <base_strategy_mixi> >), 	1, (objbase), BOOST_PP_STRINGIZE(well_rate_control_impl <TYPE <base_strategy_di> >), BOOST_PP_STRINGIZE(well_rate_control_impl <TYPE <base_strategy_mixi> >), BOOST_PP_STRINGIZE(well_rate_control_impl <TYPE <base_strategy_mixi> >), false);

  using namespace wells;
  DEFINE_WELL_RATE_CONTROL_IMPL (compute_bhp_3p_type);
  DEFINE_WELL_RATE_CONTROL_IMPL (compute_rate_3p_type);

  DEFINE_WELL_RATE_CONTROL_IMPL (compute_bhp_2p_ow_type);
  DEFINE_WELL_RATE_CONTROL_IMPL (compute_bhp_2p_og_type);
  DEFINE_WELL_RATE_CONTROL_IMPL (compute_rate_2p_ow_type);
  DEFINE_WELL_RATE_CONTROL_IMPL (compute_rate_2p_og_type);

  DEFINE_WELL_RATE_CONTROL_IMPL (compute_bhp_1p_o_type);
  DEFINE_WELL_RATE_CONTROL_IMPL (compute_bhp_1p_w_type);
  DEFINE_WELL_RATE_CONTROL_IMPL (compute_bhp_1p_g_type);
  DEFINE_WELL_RATE_CONTROL_IMPL (compute_rate_1p_o_type);
  DEFINE_WELL_RATE_CONTROL_IMPL (compute_rate_1p_w_type);
  DEFINE_WELL_RATE_CONTROL_IMPL (compute_rate_1p_g_type);

  DEFINE_WELL_RATE_CONTROL_IMPL (compute_dummy_type);

} // namespace blue_sky

