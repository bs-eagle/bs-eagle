/**
 *       \file  calc_well_pressure.cpp
 *      \brief  Implementation of calc_well_pressure
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  26.09.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"
#include "calc_well_pressure.h"
#include "calc_model.h"
#include "calc_well.h"
#include "well_connection.h"
#include "well_rate_control.h"
#include "calc_model_data_accessors.h"
#include "reservoir.h"
#include "facility_manager.h"

namespace blue_sky
  {

  /**
   * \brief  'default' ctor for calc_well_pressure
   * \param  param Additional params for ctor
   * */
  template <typename strategy_t>
  calc_well_pressure <strategy_t>::calc_well_pressure (bs_type_ctor_param /*param = NULL */)
  {
  }

  /**
   * \brief  copy-ctor for calc_well_pressure
   * \param  rhs Instance of calc_well_pressure to be copied
   * */
  template <typename strategy_t>
  calc_well_pressure <strategy_t>::calc_well_pressure (const calc_well_pressure &rhs)
  : bs_refcounter ()
  {
    *this = rhs;
  }

  template <typename strategy_t>
  bool
  calc_well_pressure<strategy_t>::calculate (sp_well_t &well, const sp_calc_model_t &calc_model) const
    {
      BS_ASSERT (!well->is_shut ()) (well->name ());

      if (well->is_rate ())
        {
          return calculate_for_rate (well, calc_model);
        }
      else
        {
          return false;
        }
    }

  /**
   * \class calc_for_rate
   * \brief Calculates BHP for well if well controlled by rate
   * */
  template <typename calc_well_pressure_t>
  struct calc_for_rate
    {
      typedef typename calc_well_pressure_t::item_t       item_t;
      typedef typename calc_well_pressure_t::index_t      index_t;
      typedef typename calc_well_pressure_t::item_array_t item_array_t;
      typedef typename calc_well_pressure_t::calc_model_t calc_model_t;
      typedef typename calc_model_t::data_t               calc_model_data_t;
      typedef typename calc_model_t::sat_d_t              sat_d_t;
      typedef typename calc_model_t::phase_d_t            phase_d_t;

      typedef typename calc_well_pressure_t::sp_well_t    sp_well_t;


      /**
       * \brief  ctor for calc_for_rate, initializes references and flags
       * \param  well
       * \param  phase_d
       * \param  sat_d
       * \param  n_phases
       * \param  gas_oil_ratio
       * */
      calc_for_rate (const sp_well_t &well, const phase_d_t &phase_d, const sat_d_t &sat_d, index_t n_phases, const item_array_t &gas_oil_ratio)
          : is_prod (well->get_well_controller ()->is_production ())
          , phase_d (phase_d)
          , sat_d (sat_d)
          , n_phases (n_phases)
          , n_block (0)
          , gas_oil_ratio (gas_oil_ratio)
      {
        // TODO: BAD DESIGN
        if (is_prod)
          {
            using namespace wells;
            rate_control_type control_type = well->get_well_controller ()->get_control_type ();
            is_w = control_type == water_rate_control || control_type == liquid_rate_control;
            is_g = control_type == gas_rate_control;
            is_o = control_type == oil_rate_control || control_type == liquid_rate_control;
          }
        else
          {
            using namespace wells;
            injection_type injection = well->get_well_controller ()->injection ();

            is_w = injection == injection_water;
            is_g = injection == injection_gas;
            is_o = injection == injection_oil;
          }

        // TODO: HACK
        is_w = phase_d[FI_PHASE_WATER] != -1 && is_w;
        is_g = phase_d[FI_PHASE_GAS] != -1 && is_g;
        is_o = phase_d[FI_PHASE_OIL] != -1 && is_o;
      }

      /**
       * \brief      Calculates denom and numer for production well
       * \param[in]  data
       * \param[out] denom
       * \param[out] numer
       * */
      void
      calc_prod (const calc_model_data_t &data, item_t &denom, item_t &numer)
      {
        if (is_w)
          {
            const item_t &mw    = MOBILITY (data, phase_d, FI_PHASE_WATER);
            item_t pcwo         = n_phases == 1 ? 0 : CAP_PRESSURE (data, phase_d, FI_PHASE_WATER);
            item_t y            = gw * mw;
            item_t z            = y * (H  - po - pcwo);

            denom += y;
            numer += z;
          }
        if (is_g)
          {
            const item_t &mg    = MOBILITY (data, phase_d, FI_PHASE_GAS);
            item_t mo           = phase_d[FI_PHASE_OIL]   == -1 ? 0 : MOBILITY (data, phase_d, FI_PHASE_OIL);
            item_t gor          = phase_d[FI_PHASE_OIL]   == -1 ? 0 : gas_oil_ratio [n_block];
            item_t pcgo         = n_phases == 1 ? 0 : CAP_PRESSURE (data, phase_d, FI_PHASE_GAS);
            item_t Po           = (H - po);
            item_t Pg           = (Po - pcgo);
            item_t y            = gw * (mg + gor * mo);
            item_t z            = gw * (mg * Pg + gor * mo * Po);

            denom += y;
            numer += z;
          }
        if (is_o)
          {
            item_t mo   = MOBILITY (data, phase_d, FI_PHASE_OIL);
            item_t y    = gw * mo;
            item_t z    = y * (H - po);

            denom += y;
            numer += z;
          }
      }

      /**
       * \brief      Calculates denom and numer for injection well
       * \param[in]  data
       * \param[out] denom
       * \param[out] numer
       * */
      void
      calc_inj (const calc_model_data_t &data, item_t &denom, item_t &numer)
      {
        item_t mw   = 0;
        item_t mg   = 0;
        item_t mo   = 0;
        item_t pcwo = 0;
        item_t pcgo = 0;

        if (n_phases == 1)
          {
            mw = phase_d[FI_PHASE_WATER] == -1 ? 0 : INVERS_VISCOSITY (data, phase_d, FI_PHASE_WATER);
            mg = phase_d[FI_PHASE_GAS]   == -1 ? 0 : INVERS_VISCOSITY (data, phase_d, FI_PHASE_GAS);
            mo = phase_d[FI_PHASE_OIL]   == -1 ? 0 : INVERS_VISCOSITY (data, phase_d, FI_PHASE_OIL);
          }
        else
          {
            mw = phase_d[FI_PHASE_WATER] == -1 ? 0 : (RELATIVE_PERM (data, phase_d, FI_PHASE_WATER) * INVERS_VISCOSITY (data, phase_d, FI_PHASE_WATER));
            mg = phase_d[FI_PHASE_GAS]   == -1 ? 0 : (RELATIVE_PERM (data, phase_d, FI_PHASE_GAS)   * INVERS_VISCOSITY (data, phase_d, FI_PHASE_GAS));
            mo = phase_d[FI_PHASE_OIL]   == -1 ? 0 : (RELATIVE_PERM (data, phase_d, FI_PHASE_OIL)   * INVERS_VISCOSITY (data, phase_d, FI_PHASE_OIL));
            pcwo = phase_d[FI_PHASE_WATER] == -1 ? 0 : CAP_PRESSURE (data, phase_d, FI_PHASE_WATER);
            pcgo = phase_d[FI_PHASE_GAS]   == -1 ? 0 : CAP_PRESSURE (data, phase_d, FI_PHASE_GAS);
          }

        if (is_w)
          {
            item_t bw   = INVERS_FVF (data, phase_d, FI_PHASE_WATER);
            mw          = mw * bw;
            mg          = mg * bw;
            mo          = mo * bw;

            item_t y    = gw * (mw + mg + mo);
            item_t z    = y * (H - po - pcwo);

            denom += y;
            numer += z;
          }
        else if (is_g)
          {
            item_t bg   = INVERS_FVF (data, phase_d, FI_PHASE_GAS);
            mw          = mw * bg;
            mg          = mg * bg;
            mo          = mo * bg;

            item_t y    = gw * (mw + mg + mo);
            item_t z    = y * (H - po - pcgo);

            denom += y;
            numer += z;
          }
        else if (is_o)
          {
            item_t bo   = INVERS_FVF (data, phase_d, FI_PHASE_OIL);
            mw          = mw * bo;
            mg          = mg * bo;
            mo          = mo * bo;

            item_t y    = gw * (mw + mg + mo);
            item_t z    = y * (H - po);

            denom += y;
            numer += z;
          }
      }

public:

      bool is_prod;
      bool is_w, is_g, is_o;

      auto_value <item_t> gw;
      auto_value <item_t> po;
      auto_value <item_t> H;

      const phase_d_t &phase_d;
      const sat_d_t &sat_d;
      
      index_t n_phases;
      index_t n_block;

      const item_array_t &gas_oil_ratio;
    };

  template <typename strategy_t>
  bool
  calc_well_pressure <strategy_t>::calculate_for_rate (sp_well_t &well, const sp_calc_model_t &calc_model) const
    {
      BS_ASSERT (!well->is_shut ()) (well->name ());

      typedef typename base_t::well_t::connection_t         connection_t;
      typedef typename base_t::well_t::sp_connection_t      sp_connection_t;

      typedef typename base_t::calc_model_t::sat_d_t        sat_d_t;
      typedef typename base_t::calc_model_t::phase_d_t      phase_d_t;

      typedef typename base_t::calc_model_t::data_t         calc_model_data_t;

      if (well->get_connections_count () == 0)
        {
          return false;
        }

      item_t gravity                = calc_model->internal_constants.gravity_constant;
      item_t numer                  = 0;
      item_t denom                  = 0;
      item_t prev_depth             = well->get_connection (0)->connection_depth;

      calc_for_rate <this_t> calc (well, calc_model->phase_d, calc_model->sat_d, calc_model->n_phases, calc_model->gas_oil_ratio);
      for (size_t i = 0, cnt = well->get_connections_count (); i < cnt; ++i)
        {
          const sp_connection_t &c = well->get_connection (i);
          if (c->is_shut ())
            {
              prev_depth                  = c->connection_depth;
              continue;
            }

          calc.n_block                    = c->n_block ();
          const calc_model_data_t &data   = calc_model->get_data (calc.n_block);
          main_var_type main_var          = calc_model->main_variable[calc.n_block];

          item_t rho                      = c->density;
          item_t diff_h                   = c->connection_depth - prev_depth;
          prev_depth                      = c->connection_depth;
          calc.gw                         = c->get_fact ();
          calc.po                         = calc_model->pressure[calc.n_block];
          calc.H                          = rho * gravity * diff_h;

          BS_ASSERT (main_var != FI_NULL) (main_var);

          if (calc.is_prod)
            {
              calc.calc_prod (data, denom, numer);
            }
          else
            {
              calc.calc_inj (data, denom, numer);
            }
        }

      item_t q_rate = well->get_input_rate ();
      item_t bhp = (q_rate - numer) / denom;
      bool switched_to_bhp = false;

      BOSOUT (section::wells, level::low) << "[" << well->name () << "] calc_well_pressure: q_rate: " << q_rate << " numer: " << numer << " denom: " << denom << " bhp: " << bhp << bs_end;
      if ((fabs (denom) < 1e-10) || (bhp < 1.0 || bhp > 1000.0))
        {
          BOSOUT (section::wells, level::low) << "[" << well->name () << "] calc_well_pressure: switched to bhp" << bs_end;

          well->get_well_controller ()->switch_to_bhp (well);
          switched_to_bhp = true;
        }
      else
        {
          BOSOUT (section::wells, level::low) << "[" << well->name () << "] calc_well_pressure: none" << bs_end;

          well->set_bhp (bhp);
        }

      return switched_to_bhp;
    }


  BLUE_SKY_TYPE_STD_CREATE_T_DEF (calc_well_pressure, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF (calc_well_pressure, (class));
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (calc_well_pressure<base_strategy_fi>), 1, (objbase), "calc_well_pressure_fi", "calc_well_pressure_fi", "calc_well_pressure_fi", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (calc_well_pressure<base_strategy_di>), 1, (objbase), "calc_well_pressure_di", "calc_well_pressure_di", "calc_well_pressure_di", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (calc_well_pressure<base_strategy_mixi>), 1, (objbase), "calc_well_pressure_mixi", "calc_well_pressure_mixi", "calc_well_pressure_mixi", false);

  bool
  calc_well_pressure_register_types (const blue_sky::plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, calc_well_pressure <base_strategy_fi>::bs_type ());
    BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, calc_well_pressure <base_strategy_di>::bs_type ());
    BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, calc_well_pressure <base_strategy_mixi>::bs_type ());
    BS_ASSERT (res);

    return res;
  }

} // namespace blue_sky

