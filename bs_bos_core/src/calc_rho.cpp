/**
 *       \file  calc_rho.cpp
 *      \brief  Implementation of calc_total_average_rho
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  26.09.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"

#include "calc_rho.h"
#include "calc_model.h"
#include "calc_well.h"
#include "well_connection.h"
#include "calc_model_data_accessors.h"
// FIXME:
//#include "reservoir.h"
//#include "facility_manager.h"

namespace blue_sky
  {

    /**
     * \brief  'default' ctor for calc_total_average_rho 
     * \param  param additional ctor params
     * */
  calc_total_average_rho::calc_total_average_rho (bs_type_ctor_param /* param = NULL */)
  {
  }

  /**
   * \brief  copy-ctor for calc_total_average_rho
   * \param  rhs calc_total_average_rho instance to be copied
   * */
  calc_total_average_rho::calc_total_average_rho (const calc_total_average_rho &rhs)
        : bs_refcounter ()
  {
    *this = rhs;
  }

  void
  calc_total_average_rho::calculate (const sp_well_t &well, const sp_calc_model_t &calc_model, 
                                                 const sp_mesh_iface_t & /*mesh*/) const
    {
      BS_ASSERT (!well->is_shut ()) (well->name ());
      BS_ASSERT (!well->is_no_connections ()) (well->name ());

      typedef well::connection_t               connection_t;
      typedef well::sp_connection_t            sp_connection_t;
      typedef calc_model::phase_d_t            phase_d_t;

      const phase_d_t &phase_d = calc_model->phase_d;
      const t_double *saturation = &(*calc_model->saturation_3p)[0];

      bool is_w = calc_model->is_water ();
      bool is_g = calc_model->is_gas ();
      bool is_o = calc_model->is_oil ();

      item_t rhop_satp = 0;
      item_t sat = 0;

      index_t i_w = phase_d[FI_PHASE_WATER];
      index_t i_g = phase_d[FI_PHASE_GAS];
      index_t i_o = phase_d[FI_PHASE_OIL];

      typename well_t::connection_iterator_t it = well->connections_begin (), 
               e = well->connections_end ();
      for (; it != e; ++it)
        {
          const sp_connection_t &c = *it;
          if (c->is_shut ())
            continue;

          index_t n_block = c->n_block ();
          const calc_model_data &data = calc_model->get_data (n_block);

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

      it = well->connections_begin ();
      for (; it != e; ++it)
        {
          const sp_connection_t &c (*it);
          if (!c->is_shut ())
            c->density = rho;
        }
    }

  BLUE_SKY_TYPE_STD_CREATE (calc_total_average_rho);
  BLUE_SKY_TYPE_STD_COPY (calc_total_average_rho);
  BLUE_SKY_TYPE_IMPL (calc_total_average_rho, objbase, "calc_total_average", "calc_total_average", "calc_total_average");

  bool
  calc_rho_register_types (const blue_sky::plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, calc_total_average_rho::bs_type ());    BS_ASSERT (res);

    return res;
  }

} // namespace blue_sky

