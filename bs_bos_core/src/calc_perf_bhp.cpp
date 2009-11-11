/**
 *       \file  calc_perf_bhp.cpp
 *      \brief  Implementation of calc_perf_bhp
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  18.11.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"

#include "calc_perf_bhp.h"
#include "calc_model.h"
#include "calc_well.h"
#include "well_connection.h"
#include "reservoir.h"
#include "facility_manager.h"

namespace blue_sky
  {

    /**
     * \brief  'default' calc_perf_bhp ctor
     * \param  param additional ctor params
     * */
  template <typename strategy_t>
  calc_perf_bhp <strategy_t>::calc_perf_bhp(bs_type_ctor_param /*param = NULL */)
  {

  }
  /**
   * \brief  copy-ctor for calc_perf_bhp
   * \param  src calc_perf_bhp instance to be copied
   * */
  template <typename strategy_t>
  calc_perf_bhp<strategy_t>::calc_perf_bhp(const calc_perf_bhp<strategy_t> & /*x*/)
        : bs_refcounter ()
  {

  }

  template <typename strategy_t>
  void
  calc_perf_bhp <strategy_t>::calculate (sp_well_t &well, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh) const
    {
      BS_ASSERT (!well->is_shut ()) (well->name ());

      typedef typename base_t::well_t::connection_t           connection_t;
      typedef typename base_t::well_t::sp_connection_t        sp_connection_t;
      typedef typename base_t::calc_model_t::sat_d_t          sat_d_t;
      typedef typename base_t::calc_model_t::phase_d_t        phase_d_t;
      typedef typename base_t::calc_model_t::data_t           calc_model_data_t;

      if (well->get_connections_count () == 0)
        {
          BOSOUT (section::wells, level::debug) 
            << "[" << well->name () << "] calc_perf_bhp: connection list is empty"
            << bs_end;
            
          return ;
        }

      item_t gravity = calc_model->internal_constants.gravity_constant;
      item_t ptop = well->bhp ();
      item_t dtop = well->get_reference_depth (mesh);

      BOSOUT (section::wells, level::low)
        << "[" << well->name () << "] calc_perf_bhp: ptop: " 
        << ptop 
        << " dtop: " << dtop << bs_end;

      sp_connection_t prev_connection;
      const typename base_t::calc_model_t::item_array_t &pressure = calc_model->pressure;
      for (size_t i = 0, cnt = well->get_connections_count (); i < cnt; ++i)
        {
          const sp_connection_t &c (well->get_connection (i));

          if (!c->is_shut ())
            {
              if (prev_connection)
                {
                  c->set_head_term (prev_connection->head_term + 0.5 * (prev_connection->density + c->density) * gravity * (c->connection_depth - prev_connection->connection_depth));
                }
              else
                {
                  c->set_head_term (c->density * gravity * (c->connection_depth - dtop));
                }
              c->set_cur_bhp (ptop + c->head_term);
              c->set_bulkp (pressure[c->n_block ()]);
#ifdef _DEBUG
              BOSOUT (section::wells, level::low)
                << "[" << well->name () << " : " << c->n_block () << "] calc_perf_bhp: density: " 
                << c->density 
                << " connection_depth: " << c->connection_depth
                << " head_term: " << c->head_term
                << " cur_bhp: " << c->cur_bhp
                << bs_end;
#endif
            }

          prev_connection = c;
        }
    }

  BLUE_SKY_TYPE_STD_CREATE_T_DEF (calc_perf_bhp, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF (calc_perf_bhp, (class));
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (calc_perf_bhp <base_strategy_fi>), 1, (objbase), "calc_perf_bhp_fi", "calc_perf_bhp_fi", "calc_perf_bhp_fi", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (calc_perf_bhp <base_strategy_di>), 1, (objbase), "calc_perf_bhp_di", "calc_perf_bhp_di", "calc_perf_bhp_di", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (calc_perf_bhp <base_strategy_mixi>), 1, (objbase), "calc_perf_bhp_mixi", "calc_perf_bhp_mixi", "calc_perf_bhp_mixi", false);

} // namespace blue_sky

