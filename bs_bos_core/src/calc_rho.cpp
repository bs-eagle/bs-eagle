/**
 *
 * */
#include "stdafx.h"

#include "calc_rho.h"
#include "calc_model.h"
#include "calc_well.h"
#include "well_connection.h"
#include "calc_model_data_accessors.h"
#include "reservoir.h"
#include "facility_manager.h"

namespace blue_sky
  {

  template <typename strategy_t>
  calc_total_average_rho <strategy_t>::calc_total_average_rho (bs_type_ctor_param /* param = NULL */)
  {
  }

  template <typename strategy_t>
  calc_total_average_rho <strategy_t>::calc_total_average_rho (const calc_total_average_rho &rhs)
        : bs_refcounter ()
  {
    *this = rhs;
  }

  template <typename strategy_t>
  void
  calc_total_average_rho<strategy_t>::calculate (const sp_well_t &well, const sp_calc_model_t &calc_model, 
                                                 const sp_mesh_iface_t & /*mesh*/) const
    {
      BS_ASSERT (!well->is_shut ()) (well->name ());
      BS_ASSERT (well->get_connections_count ()) (well->name ());

      typedef typename well_t::connection_t               connection_t;
      typedef typename well_t::sp_connection_t            sp_connection_t;
      typedef typename calc_model_t::data_t               calc_model_data_t;
      typedef typename calc_model_t::phase_d_t            phase_d_t;

      const phase_d_t &phase_d = calc_model->phase_d;
      const item_array_t &saturation = calc_model->saturation_3p;

      bool is_w = calc_model->is_water ();
      bool is_g = calc_model->is_gas ();
      bool is_o = calc_model->is_oil ();

      item_t rhop_satp = 0;
      item_t sat = 0;

      index_t i_w = phase_d[FI_PHASE_WATER];
      index_t i_g = phase_d[FI_PHASE_GAS];
      index_t i_o = phase_d[FI_PHASE_OIL];

      for (size_t i = 0, cnt = well->get_connections_count (); i < cnt; ++i)
        {
          const sp_connection_t &c = well->get_connection (i);
          if (c->is_shut ())
            continue;

          index_t n_block = c->n_block ();
          const calc_model_data_t &data = calc_model->get_data (n_block);

          if (is_w)
            {
              rhop_satp += DENSITY (data, phase_d, FI_PHASE_WATER) * saturation[i_w];
              sat += saturation[i_w];
            }
          if (is_g)
            {
              rhop_satp += DENSITY (data, phase_d, FI_PHASE_GAS) * saturation[i_g];
              sat += saturation[i_g];
            }
          if (is_o)
            {
              rhop_satp += DENSITY (data, phase_d, FI_PHASE_OIL) * saturation[i_o];
              sat += saturation[i_o];
            }
        }

      BS_ASSERT (sat != 0);
      item_t rho = rhop_satp / sat;

      for (size_t i = 0, cnt = well->get_connections_count (); i < cnt; ++i)
        {
          const sp_connection_t &c (well->get_connection (i));
          if (!c->is_shut ())
            c->density = rho;
        }
    }

  BLUE_SKY_TYPE_STD_CREATE_T_DEF (calc_total_average_rho, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF (calc_total_average_rho, (class));
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (calc_total_average_rho<base_strategy_fi>), 1, (objbase), "calc_total_average_fi", "calc_total_average_fi", "calc_total_average_fi", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (calc_total_average_rho<base_strategy_di>), 1, (objbase), "calc_total_average_di", "calc_total_average_di", "calc_total_average_di", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (calc_total_average_rho<base_strategy_mixi>), 1, (objbase), "calc_total_average_mixi", "calc_total_average_mixi", "calc_total_average_mixi", false);

  bool
  calc_rho_register_types (const blue_sky::plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, calc_total_average_rho <base_strategy_fi>::bs_type ());
    BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, calc_total_average_rho <base_strategy_di>::bs_type ());
    BS_ASSERT (res);

    return res;
  }

} // namespace blue_sky

