/**
 *       \file  wellbore_density_calc.cpp
 *      \brief  Implementation of wellbore_density_calc
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  18.11.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"

#include "wellbore_density_calc.h"

#include "calc_model.h"
#include "calc_well.h"
#include "well_connection.h"

#include "reservoir.h"
#include "facility_manager.h"
namespace blue_sky
  {

    /**
     * \brief  'default' ctor for wellbore_density_calc
     * \param  param Additional parameters for ctor
     * */
  template <typename strategy_t>
  wellbore_density_calc<strategy_t>::wellbore_density_calc (bs_type_ctor_param param /* = NULL */)
  {

  }
  /**
   * \brief  copy-ctor for wellbore_density_calc
   * \param  x Instance of wellbore_density_calc to be copied
   * */
  template <typename strategy_t>
  wellbore_density_calc<strategy_t>::wellbore_density_calc (const wellbore_density_calc &x)
  : bs_refcounter (x)
  {

  }

  template <typename strategy_t>
  void
  wellbore_density_calc <strategy_t>::calculate (sp_well_t &well, const sp_calc_model_t &calc_model) const
    {
      BS_ASSERT (!well->is_shut ()) (well->name ());

      typedef typename base_t::well_t                     well_t;
      typedef typename base_t::calc_model_t               calc_model_t;
      typedef typename well_t::connection_t               connection_t;
      typedef typename well_t::sp_connection_t            sp_connection_t;
      typedef typename well_t::sp_connection_t            sp_connection_t;

      item_t atm_pressure = calc_model->internal_constants.atmospheric_pressure;

      for (size_t i = 0, cnt = well->get_connections_count (); i < cnt; ++i)
        {
          const sp_connection_t &c (well->get_connection (i));

          item_t pbhp = calc_model->pressure[c->n_block ()];
          if (c->cur_bhp >= atm_pressure)
            {
              pbhp = c->cur_bhp;
            }
          else if (well->bhp () >= atm_pressure)
            {
              pbhp = well->bhp ();
            }

          if (density_calc (well, calc_model, pbhp))
            return ;
        }

      // TODO: throw exception
      return ;
    }

  template <typename strategy_t>
  bool
  wellbore_density_calc <strategy_t>::density_calc (sp_well_t &well, const sp_calc_model_t &calc_model, item_t p_bhp) const
    {
      BS_ASSERT (!well->is_shut ()) (well->name ());
      BS_ASSERT (well->get_connections_count ()) (well->name ());

      typedef typename base_t::index_t                          index_t;
      typedef typename base_t::well_t::sp_connection_t          sp_connection_t;
      typedef typename base_t::calc_model_t::sp_pvt_water       sp_pvt_water_t;
      typedef typename base_t::calc_model_t::sp_pvt_gas         sp_pvt_gas_t;
      typedef typename base_t::calc_model_t::sp_pvt_dead_oil    sp_pvt_dead_oil_t;
      typedef typename base_t::well_t::sp_connection_t          sp_connection_t;

      item_t gravity_const  = calc_model->internal_constants.gravity_constant;
      index_t n_phases      = calc_model->n_phases;

      index_t d_w           = calc_model->phase_d[FI_PHASE_WATER];
      index_t d_g           = calc_model->phase_d[FI_PHASE_GAS];
      index_t d_o           = calc_model->phase_d[FI_PHASE_OIL];
      index_t d_s_w         = calc_model->sat_d[FI_PHASE_WATER];
      index_t d_s_g         = calc_model->sat_d[FI_PHASE_GAS];

      bool is_w             = calc_model->is_water ();
      bool is_g             = calc_model->is_gas ();
      bool is_o             = calc_model->is_oil ();

      boost::array <item_t, FI_PHASE_TOT> local_invers_fvf;
      local_invers_fvf.assign (0);

      BS_ERROR (well->get_connections_count (), "density_calc");

      const sp_connection_t &first_con = well->get_connection (0);
      BS_ASSERT (first_con);

      item_t ppo = 0;
      if (!p_bhp)
        {
          // Average pressure is equal BHP at first connection ...
          if (first_con->cur_bhp)
            ppo = first_con->cur_bhp;
          // or reservoir pressure at connection
          else if (well->bhp ())
            ppo = well->bhp ();
          else
            ppo = calc_model->pressure[first_con->n_block ()];
        }
      else
        ppo = p_bhp;

      item_t ppo1     = ppo;
      item_t p_water  = 0;
      item_t p_gas    = 0;
      item_t p_oil    = 0;
      item_t rrso     = 0.0f;

      item_t sat_w    = 0;
      item_t sat_g    = 0;
      item_t sat_o    = 0;

      typename base_t::well_t::sp_well_controller_t well_controller_ = well->get_well_controller ();
      wells::injection_type injection = well_controller_->injection ();

      const rate_data <strategy_t> &rate = well->rate ();

      sp_connection_t prev_con;
      for (size_t i = 0, cnt = well->get_connections_count (); i < cnt; ++i)
        {
          const sp_connection_t &con (well->get_connection (i));

          index_t n_block   = con->n_block ();
          index_t n_pvt_reg = calc_model->pvt_regions[n_block];
          index_t i_w       = n_block * n_phases + d_w;
          index_t i_g       = n_block * n_phases + d_g;
          index_t i_o       = n_block * n_phases + d_o;

          item_t con_denom  = 0;
          item_t con_upnom  = 0;

          if (prev_con && prev_con->cur_bhp)
            ppo1 = prev_con->cur_bhp + 0.5 * (prev_con->density + con->density) * gravity_const *
                   (con->connection_depth - prev_con->connection_depth);
          else
            ppo1 = ppo;

          sp_pvt_water_t pvt_water;
          sp_pvt_gas_t pvt_gas;
          sp_pvt_dead_oil_t pvt_oil;

          if (is_w)
            pvt_water = calc_model->pvt_water_array [n_pvt_reg];
          if (is_g)
            pvt_gas = calc_model->pvt_gas_array [n_pvt_reg];
          if (is_o)
            pvt_oil = calc_model->pvt_oil_array [n_pvt_reg];

          if (ppo1)
            p_oil = ppo1;
          else
            p_oil = ppo;
          p_water = p_gas = p_oil;

          if (n_phases > 1)
            {
              // two or three phases
              // water phase
              if (is_w)
                p_water += calc_model->data[n_block].cap_pressure[d_s_w];
              // gas phase
              if (is_g)
                p_gas += calc_model->data[n_block].cap_pressure[d_s_g];
            }

          // calculate 1 / (formation volume factor), 1 / (viscosity), 1 / (viscosity) / (fvf)
          // and derivates
          // water phase
          // calculate pvt properties
          if (is_w)
            {
              pvt_water->calc (p_water, &local_invers_fvf[d_w], 0, 0, 0, 0, 0);
            }
          if (is_g)
            {
              pvt_gas->calc (p_gas, &local_invers_fvf[d_g], 0, 0, 0, 0, 0);
            }
          if (is_o && is_g)
            {
              rrso = calc_model->gas_oil_ratio[n_block];
              if (!pvt_oil->calc (is_g, calc_model->main_variable[n_block],
                                  p_oil, rrso,
                                  &local_invers_fvf[d_o], 0, 0,
                                  0, 0, 0,
                                  0, 0, 0,
                                  0, 0))
                {
                  BOSERR (section::wells, level::error) << "pvt_oil (o && g) can't calculate (" << p_oil << ", " << rrso << ")" << bs_end;
                  BOSERR (section::wells, level::error) << "p_bhp: " << p_bhp << bs_end;
                  if (first_con->cur_bhp)
                    BOSERR (section::wells, level::error) << "first_con: " << first_con->cur_bhp << bs_end;
                  else if (well->bhp ())
                    BOSERR (section::wells, level::error) << "well: " << well->bhp () << bs_end;
                  else
                    BOSERR (section::wells, level::error) << "pressure: " << calc_model->pressure[first_con->n_block ()] << bs_end;
                  // TODO: BUG:
                  //throw bs_exception ("wellbore_density_calc", "pvt->oil (o && g) can't calc");
                  return false;
                }
            }
          else if (is_o)
            {
              if (!pvt_oil->calc (is_g, -1, p_oil, 0,
                             &local_invers_fvf[d_o], 0, 0,
                             0, 0, 0,
                             0, 0, 0,
                             0, 0))
                {
                  BOSERR (section::wells, level::error) << "pvt_oil (o) can't calculate (" << p_oil << ", " << rrso << ")" << bs_end;
                  BOSERR (section::wells, level::error) << "p_bhp: " << p_bhp << bs_end;
                  if (first_con->cur_bhp)
                    BOSERR (section::wells, level::error) << "first_con: " << first_con->cur_bhp << bs_end;
                  else if (well->bhp ())
                    BOSERR (section::wells, level::error) << "well: " << well->bhp () << bs_end;
                  else
                    BOSERR (section::wells, level::error) << "pressure: " << calc_model->pressure[first_con->n_block ()] << bs_end;
                  // TODO: BUG:
                  //throw bs_exception ("wellbore_density_calc", "pvt->oil (o) can't calc");
                  return false;
                }
            }

          // production well
          if (well_controller_->is_production ())
            {
              // water phase
              if (is_w)
                {
                  con_denom += rate.prod.water / local_invers_fvf[d_w];
                  con_upnom += pvt_water->get_surface_density () * rate.prod.water / local_invers_fvf[d_w];
                }

              // gas phase
              if (is_g)
                {
                  con_denom += (rate.prod.gas - rrso * rate.prod.oil) / local_invers_fvf[d_g];
                  con_upnom += pvt_gas->get_surface_density () * (rate.prod.gas - rrso * rate.prod.oil) / local_invers_fvf[d_g];
                }

              // oil phase
              if (is_o)
                {
                  con_denom += rate.prod.oil / local_invers_fvf[d_o];
                  con_upnom += pvt_oil->get_surface_density () * rate.prod.oil / local_invers_fvf[d_o];
                }

              // calculate density for current connection
              if (fabs (con_denom) > EPS_DIV)
                {
                  con->density = con_upnom / con_denom;
                }
              else
                {
                  con->density = 0.0;
                  sat_w = 0.0;
                  sat_g = 0.0;
                  sat_o = 0.0;

                  if (n_phases > 1)
                    {
                      if (is_w)
                        {
                          sat_w = calc_model->saturation_3p[i_w];
                          con->density += sat_w * pvt_water->get_surface_density () * local_invers_fvf[d_w];
                        }

                      if (is_g && is_o)
                        {
                          sat_g = calc_model->saturation_3p[i_g];
                          sat_o = calc_model->saturation_3p[i_o];
                          con->density += sat_g * pvt_gas->get_surface_density () * local_invers_fvf[d_g] +
                                          sat_o * (pvt_oil->get_surface_density () +
                                                   calc_model->gas_oil_ratio[n_block] * pvt_gas->get_surface_density ()) * local_invers_fvf[d_o];
                        }
                      else if (is_o)
                        {
                          sat_o = calc_model->saturation_3p[i_o];
                          con->density += sat_o * pvt_oil->get_surface_density () * local_invers_fvf[d_o];
                        }
                      else if (is_g)
                        {
                          sat_g = calc_model->saturation_3p[i_g];
                          con->density += sat_g * pvt_gas->get_surface_density () * local_invers_fvf[d_g];
                        }
                    }
                  else
                    {
                      // one phase only
                      if (is_w)
                        {
                          con->density = pvt_water->get_surface_density () * local_invers_fvf[d_w];
                        }
                      if (is_g)
                        {
                          con->density = pvt_gas->get_surface_density () * local_invers_fvf[d_g];
                        }
                      if (is_o)
                        {
                          con->density = pvt_oil->get_surface_density () * local_invers_fvf[d_o];
                        }
                    }
                }
            }
          // injection well
          else //if (well->current_status > 0)
            {
              if (injection == wells::injection_water)
                {
                  if (is_w)
                    con->density = pvt_water->get_surface_density () * local_invers_fvf[d_w];
                  else
                    return false;
                }
              else if (injection == wells::injection_gas)
                {
                  if (is_g)
                    con->density = pvt_gas->get_surface_density () * local_invers_fvf[d_g];
                  else
                    return false;
                }
              else if (injection == wells::injection_oil)
                {
                  if (is_o)
                    con->density = pvt_oil->get_surface_density () * local_invers_fvf[d_o];
                  else
                    return false;
                }
            }
        }

      return true;
    }

  BLUE_SKY_TYPE_STD_CREATE_T_DEF (wellbore_density_calc, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF (wellbore_density_calc, (class));
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (wellbore_density_calc <base_strategy_fi>), 1, (objbase), "wellbore_density_calc_fi", "wellbore_density_calc_fi", "wellbore_density_calc_fi", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (wellbore_density_calc<base_strategy_di>), 1, (objbase), "wellbore_density_calc_di", "wellbore_density_calc_di", "wellbore_density_calc_di", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (wellbore_density_calc<base_strategy_mixi>), 1, (objbase), "wellbore_density_calc_mixi", "wellbore_density_calc_mixi", "wellbore_density_calc_mixi", false);


} // namespace blue_sky

