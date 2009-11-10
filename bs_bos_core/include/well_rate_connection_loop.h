/**
 * \file well_rate_connection_loop.h
 * \brief loop through connection for well derivs computations
 * \author Sergey Miryanov
 * \date 21.11.2008
 * */
#ifndef BS_WELLS_WELL_RATE_CONTROL_CONNECTION_LOOP_H_
#define BS_WELLS_WELL_RATE_CONTROL_CONNECTION_LOOP_H_

namespace blue_sky
  {

  template <typename params_t, typename inj_inner_loop_t, typename prod_inner_loop_t>
  inline void
  connection_loop (params_t &params, const inj_inner_loop_t &inj_inner_loop, const prod_inner_loop_t &prod_inner_loop)
  {
    BS_ASSERT (!params.well_->is_shut ()) (params.well_->name ());
    if (!params.well_->get_connections_count ())
      {
        BOSOUT (section::wells, level::debug)
          << "[" << params.well_->name () << "] connection loop: connection list is empty"
          << bs_end;

        return ;
      }

    typedef typename params_t::sp_connection_t      sp_connection_t;
    typedef typename params_t::data_t               data_t;

    params.depth = params.well_->get_connection (0)->connection_depth;
    for (size_t i = 0, cnt = params.well_->get_connections_count (); i < cnt; ++i)
      {
        const sp_connection_t &c (params.well_->get_connection (i));
        if (c->is_shut ())
          {
            params.depth    = c->connection_depth;
            continue;
          }

        params.n_block      = c->n_block ();
        params.perf_bhp     = c->cur_bhp;
        params.diff_depth   = c->connection_depth - params.depth;
        params.depth        = c->connection_depth;
        params.rho          = c->density;
        params.gw           = c->get_fact ();
        params.main_var     = params.main_vars[params.n_block];
        const data_t &data  = params.data_array[params.n_block];

        //
        compute_potentials (data, params);

        // TODO: check crossflow
        //if (params.is_prod)
        if (is_prod_potential (params))
          {
            prod_inner_loop (c, data, params);
          }
        else
          {
            params.compute_perf_vars (data, params.inj_params_);
            inj_inner_loop (c, data, params);
          }
      }
  }

  template <typename params_t, typename inner_loop_t>
  inline void
  update_wr_connection_loop (params_t &params, const inner_loop_t &inner_loop)
  {
    BS_ASSERT (!params.well_->is_shut ()) (params.well_->name ());
    if (!params.well_->get_connections_count ())
      {
        BOSOUT (section::wells, level::debug)
          << "[" << params.well_->name () << "] connection loop: connection list is empty"
          << bs_end;

        return ;
      }

    typedef typename params_t::sp_connection_t      sp_connection_t;
    typedef typename params_t::item_t               item_t;

    //item_t ww = 1.0 / params.ww_value[0];
    item_t ww = params.ww_value[0];
    if (fabs (ww) >= 10e-16)
      {
        ww = 1.0 / ww;
      }

    for (size_t i = 0, cnt = params.well_->get_connections_count (); i < cnt; ++i)
      {
        const sp_connection_t &c (params.well_->get_connection (i));
        if (c->is_shut ())
          continue;

        inner_loop (c, ww);
      }
  }

  template <typename params_t, typename inner_loop_t>
  inline void 
  apply_wefac_connection_loop (params_t &params, const inner_loop_t &inner_loop)
  {
    BS_ASSERT (!params.well_->is_shut ()) (params.well_->name ());
    if (!params.well_->get_connections_count ())
      {
        BOSOUT (section::wells, level::debug)
          << "[" << params.well_->name () << "] connection loop: connection list is empty"
          << bs_end;

        return ;
      }

    typedef typename params_t::sp_connection_t      sp_connection_t;
    typedef typename params_t::item_t               item_t;

    for (size_t i = 0, cnt = params.well_->get_connections_count (); i < cnt; ++i)
      {
        const sp_connection_t &c (params.well_->get_connection (i));
        if (c->is_shut ())
          continue;

        inner_loop (c, params);
      }

  }

} // namespace blue_sky


#endif // #ifndef BS_WELLS_WELL_RATE_CONTROL_CONNECTION_LOOP_H_

