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
// FIXME:
//#include "reservoir.h"
//#include "facility_manager.h"

namespace blue_sky
  {

    /**
     * \brief  'default' calc_perf_bhp ctor
     * \param  param additional ctor params
     * */
  calc_perf_bhp::calc_perf_bhp(bs_type_ctor_param /*param = NULL */)
  {

  }
  /**
   * \brief  copy-ctor for calc_perf_bhp
   * \param  src calc_perf_bhp instance to be copied
   * */
  calc_perf_bhp::calc_perf_bhp(const calc_perf_bhp& /*x*/)
        : bs_refcounter ()
  {

  }

  void
  calc_perf_bhp::calculate (sp_well_t &well, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh) const
    {
      BS_ASSERT (!well->is_shut ()) (well->name ());

      typedef base_t::well_t::connection_t           connection_t;
      typedef base_t::well_t::sp_connection_t        sp_connection_t;
      typedef base_t::calc_model_t::sat_d_t          sat_d_t;
      typedef base_t::calc_model_t::phase_d_t        phase_d_t;
      typedef base_t::calc_model_t::data_t           calc_model_data_t;

      if (well->is_no_connections ())
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
      const t_double *pressure = &(*calc_model->pressure)[0];
      base_t::well_t::connection_iterator_t it = well->connections_begin (), e = well->connections_end ();
      for (; it != e; ++it)
        {
          const sp_connection_t &c (*it);

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

  BLUE_SKY_TYPE_STD_CREATE (calc_perf_bhp);
  BLUE_SKY_TYPE_STD_COPY (calc_perf_bhp);
  BLUE_SKY_TYPE_IMPL (calc_perf_bhp, objbase, "calc_perf_bhp", "calc_perf_bhp", "calc_perf_bhp");

} // namespace blue_sky

