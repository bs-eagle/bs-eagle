/**
 *       \file  default_well_rate_control_factory.cpp
 *      \brief  Impementation of default_well_rate_control_factory
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  21.05.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Obsolete, should be removed
 * */
#include "stdafx.h"
#include "default_well_rate_control_factory.h"

#include "calc_well.h"
#include "well_controller.h"
#include "well_connection.h"

#include "calc_model_data_accessors.h"

#include "well_rate_control.h"
#include "well_rate_control_impl.h"

#include "well_rate_control_prod_mobility.h"
#include "well_rate_control_inj_mobility.h"

#include "well_rate_control_rate.h"
#include "well_rate_control_bhp_deriv.h"
#include "well_rate_control_rate_deriv.h"

#include "well_rate_control_impl_type.h"
#include "reservoir.h"
#include "facility_manager.h"

namespace blue_sky {
namespace wells {

  template <typename strategy_t>
  default_well_rate_control_factory <strategy_t>::default_well_rate_control_factory (bs_type_ctor_param param)
  : base_t (param)
  {
  }
  template <typename strategy_t>
  default_well_rate_control_factory <strategy_t>::default_well_rate_control_factory (const default_well_rate_control_factory &f)
        : bs_refcounter (), base_t (f)
  {
  }

  template <typename strategy_t>
  typename default_well_rate_control_factory <strategy_t>::sp_well_rate_control_t
  default_well_rate_control_factory <strategy_t>::create_control (rate_control_type control_type, bool is_bhp, 
                                                                  bool /*is_production*/, 
                                                                  const sp_calc_model_t &calc_model)
  {
      typedef typename strategy_t::index_t index_t;

      if (control_type == null_control)
        {
          return BS_KERNEL.create_object (well_rate_control_impl <compute_dummy_type <strategy_t> >::bs_type (), true);
        }


      index_t n_phases = calc_model->n_phases;
      bool is_o = calc_model->is_oil ();
      bool is_w = calc_model->is_water ();
      bool is_g = calc_model->is_gas ();

      if (n_phases == 3)
        {
          return is_bhp
                 ? BS_KERNEL.create_object (well_rate_control_impl <compute_bhp_3p_type <strategy_t> >::bs_type (), true)
                 : BS_KERNEL.create_object (well_rate_control_impl <compute_rate_3p_type <strategy_t> >::bs_type (), true)
                 ;
        }
      else if (is_o && is_w)
        {
          return is_bhp
                 ? BS_KERNEL.create_object (well_rate_control_impl <compute_bhp_2p_ow_type <strategy_t> >::bs_type (), true)
                 : BS_KERNEL.create_object (well_rate_control_impl <compute_rate_2p_ow_type <strategy_t> >::bs_type (), true)
                 ;
        }
      else if (is_o && is_g)
        {
          return is_bhp
                 ? BS_KERNEL.create_object (well_rate_control_impl <compute_bhp_2p_og_type <strategy_t> >::bs_type (), true)
                 : BS_KERNEL.create_object (well_rate_control_impl <compute_rate_2p_og_type <strategy_t> >::bs_type (), true)
                 ;
        }
      else if (is_o)
        {
          return is_bhp
                 ? BS_KERNEL.create_object (well_rate_control_impl <compute_bhp_1p_o_type <strategy_t> >::bs_type (), true)
                 : BS_KERNEL.create_object (well_rate_control_impl <compute_rate_1p_o_type <strategy_t> >::bs_type (), true)
                 ;
        }
      else if (is_w)
        {
          return is_bhp
                 ? BS_KERNEL.create_object (well_rate_control_impl <compute_bhp_1p_w_type <strategy_t> >::bs_type (), true)
                 : BS_KERNEL.create_object (well_rate_control_impl <compute_rate_1p_w_type <strategy_t> >::bs_type (), true)
                 ;
        }
      else if (is_g)
        {
          return is_bhp
                 ? BS_KERNEL.create_object (well_rate_control_impl <compute_bhp_1p_g_type <strategy_t> >::bs_type (), true)
                 : BS_KERNEL.create_object (well_rate_control_impl <compute_rate_1p_g_type <strategy_t> >::bs_type (), true)
                 ;
        }
      else
        {
          bs_throw_exception ("can't create well_rate_control_impl");
        }
  }

  BLUE_SKY_TYPE_STD_CREATE_T_DEF (default_well_rate_control_factory, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF (default_well_rate_control_factory, (class));
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (default_well_rate_control_factory <base_strategy_fi>), 1, (objbase), "default_well_rate_control_fi", "Default factory for well_rate_control", "Default factory for well_rate_control", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (default_well_rate_control_factory <base_strategy_di>), 1, (objbase), "default_well_rate_control_di", "Default factory for well_rate_control", "Default factory for well_rate_control", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (default_well_rate_control_factory <base_strategy_mixi>), 1, (objbase), "default_well_rate_control_mixi", "Default factory for well_rate_control", "Default factory for well_rate_control", false);

} // namespace wells
} // namespace blue_sky

