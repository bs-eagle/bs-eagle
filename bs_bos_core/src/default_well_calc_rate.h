/**
 * \file default_well_calc_rate.h
 * \brief calculate rates for default_well
 * \author Sergey Miryanov
 * \date 13.10.2009
 * */
#ifndef BS_BOS_CORE_DEFAULT_WELL_CALC_RATE_H_
#define BS_BOS_CORE_DEFAULT_WELL_CALC_RATE_H_

#include "matrix_vector_op.h"
#include "pp_index.h"

#define RELATIVE_PERM_W data.relative_perm[d_w]
#define RELATIVE_PERM_G data.relative_perm[d_g]
#define RELATIVE_PERM_O data.relative_perm[d_o]
#define P_DERIV_INVERS_VISC_W data.p_deriv_invers_viscosity[d_w]
#define P_DERIV_INVERS_VISC_G data.p_deriv_invers_viscosity[d_g]
#define P_DERIV_INVERS_VISC_O data.p_deriv_invers_viscosity[d_o]
#define S_DERIV_RELATIVE_PERM_WW data.s_deriv_relative_perm[d_w * n_phases + d_w]
#define S_DERIV_RELATIVE_PERM_GG data.s_deriv_relative_perm[d_g * n_phases + d_g]
#define S_DERIV_RELATIVE_PERM_GO data.s_deriv_relative_perm[d_g * n_phases + d_o]
#define S_DERIV_RELATIVE_PERM_OW data.s_deriv_relative_perm[d_o * n_phases + d_w]
#define S_DERIV_RELATIVE_PERM_OG data.s_deriv_relative_perm[d_o * n_phases + d_g]
#define S_DERIV_RELATIVE_PERM_OO data.s_deriv_relative_perm[d_o * n_phases + d_o]
#define S_DERIV_CAP_PRESSURE_W data.s_deriv_cap_pressure[d_w]
#define S_DERIV_CAP_PRESSURE_G data.s_deriv_cap_pressure[d_g]
#define INVERS_VISC_W data.invers_viscosity[d_w]
#define INVERS_VISC_G data.invers_viscosity[d_g]
#define INVERS_VISC_O data.invers_viscosity[d_o]
#define GOR_DERIV_INVERS_VISC data.gor_deriv_invers_viscosity
#define P_DERIV_INVERS_FVF_W data.p_deriv_invers_fvf[d_w]
#define P_DERIV_INVERS_FVF_G data.p_deriv_invers_fvf[d_g]
#define P_DERIV_INVERS_FVF_O data.p_deriv_invers_fvf[d_o]
#define INVERS_FVF_W data.invers_fvf[d_w]
#define INVERS_FVF_G data.invers_fvf[d_g]
#define INVERS_FVF_O data.invers_fvf[d_o]
#define S_DERIV_INVERS_FVF_W data.p_deriv_invers_fvf[d_w] * data.s_deriv_cap_pressure[d_w] 
#define S_DERIV_INVERS_FVF_G data.p_deriv_invers_fvf[d_g] * data.s_deriv_cap_pressure[d_g] 
#define S_DERIV_INVERS_VISC_W data.p_deriv_invers_viscosity[d_w] * data.s_deriv_cap_pressure[d_w]
#define S_DERIV_INVERS_VISC_G data.p_deriv_invers_viscosity[d_g] * data.s_deriv_cap_pressure[d_g]
#define CAP_PRESSURE_W data.cap_pressure[d_w]
#define CAP_PRESSURE_G data.cap_pressure[d_g]

namespace blue_sky {
namespace wells {

  namespace detail {
    template <bool is_prod>
    struct rate_data
    {
    };

    template <>
    struct rate_data <true>
    {
      template <typename T>
      static typename T::rate_data_inner &
      get (T &t)
      {
        return t.prod;
      }
    };
    template <>
    struct rate_data <false>
    {
      template <typename T>
      static typename T::rate_data_inner &
      get (T &t)
      {
        return t.inj;
      }
    };

    template <bool is_prod, main_var_type main_var>
    struct mobility
    {
    };

    template <main_var_type main_var>
    struct mobility <true, main_var>
    {
      template <typename data_t, typename index_t, typename T>
      static typename data_t::item_t
      water (const data_t &data, index_t d_w, T *)
      {
        return data.mobility[d_w];
      }
      template <typename data_t, typename index_t, typename T>
      static typename data_t::item_t
      oil (const data_t &data, index_t d_o, T *)
      {
        return data.mobility[d_o];
      }
      template <typename data_t, typename index_t, typename T>
      static typename data_t::item_t
      gas (const data_t &, index_t, T *)
      {
        return 0;
      }
    };

    template <main_var_type main_var>
    struct mobility <false, main_var>
    {
      template <typename data_t, typename index_t, typename T>
      static typename data_t::item_t
      get_mobility (const data_t &data, index_t d, T *t)
      {
        if (main_var == FI_SG_VAR)
          {
            return data.invers_fvf[d] * t->krp_tetap;
          }
        else
          {
            return data.invers_fvf[d] * t->krow_tetaow;
          }
      }

      template <typename data_t, typename index_t, typename T>
      static typename data_t::item_t
      water (const data_t &data, index_t d_w, T *t)
      {
        return get_mobility (data, d_w, t);
      }
      template <typename data_t, typename index_t, typename T>
      static typename data_t::item_t
      oil (const data_t &data, index_t d_o, T *t)
      {
        return get_mobility (data, d_o, t);
      }
      template <typename data_t, typename index_t, typename T>
      static typename data_t::item_t
      gas (const data_t &data, index_t d_g, T *t)
      {
        return get_mobility (data, d_g, t);
      }
    };
  }

  template <typename strategy_t, bool is_w, bool is_g, bool is_o, bool is_production_well>
  struct calc_rate_and_derivs_t
  {
    typedef typename strategy_t::item_t                 item_t;
    typedef typename strategy_t::index_t                index_t;
    typedef typename strategy_t::item_array_t           item_array_t;
    typedef typename strategy_t::rhs_item_array_t       rhs_item_array_t;

    typedef default_well <strategy_t>                   well_t;
    typedef default_connection <strategy_t>             connection_t;
    typedef typename well_t::connection_list_t          connection_list_t;
    typedef typename well_t::base_t::rate_data_t        rate_data_t;
    typedef typename rate_data_t::rate_data_inner       rate_data_inner_t;
    typedef calc_model <strategy_t>                     calc_model_t;
    typedef calc_model_data <strategy_t>                calc_model_data_t;
    typedef typename calc_model_t::data_array_t         cell_data_t;
    typedef typename calc_model_t::main_var_array_t     main_vars_t;
    typedef typename main_vars_t::value_type            main_var_t;
    typedef jacobian_matrix <strategy_t>                jmatrix_t;

    typedef smart_ptr <calc_model_t, true>              sp_calc_model_t;
    typedef smart_ptr <jmatrix_t, true>                 sp_jmatrix_t;

    enum
      {
        n_phases = is_w + is_g + is_o,
        b_sqr = n_phases * n_phases,
      };

    enum
      {
        is_1p = n_phases == 1,
        is_2p = n_phases == 2,
        is_3p = n_phases == 3,
      };

      enum
        {
          gas_sg = blue_sky::detail::pp_index <n_phases, is_w, is_g>::gas_sg,
          gas_so = blue_sky::detail::pp_index <n_phases, is_w, is_g>::gas_so,
          gas_po = blue_sky::detail::pp_index <n_phases, is_w, is_g>::gas_po,
          oil_sg = blue_sky::detail::pp_index <n_phases, is_w, is_g>::oil_sg,
          oil_so = blue_sky::detail::pp_index <n_phases, is_w, is_g>::oil_so,
          oil_po = blue_sky::detail::pp_index <n_phases, is_w, is_g>::oil_po,
          wat_sg = blue_sky::detail::pp_index <n_phases, is_w, is_g>::wat_sg,
          wat_so = blue_sky::detail::pp_index <n_phases, is_w, is_g>::wat_so,
          wat_po = blue_sky::detail::pp_index <n_phases, is_w, is_g>::wat_po,
        };

      enum
        {
          gas_idx = 0,
          oil_idx = is_g,
          wat_idx = is_g + is_o,

          gas_sw = !is_1p ? b_sqr + gas_idx : -1,
          oil_sw = !is_1p ? b_sqr + oil_idx : -1,
          wat_sw = !is_1p ? b_sqr + wat_idx : -1,
        };

    calc_rate_and_derivs_t (const sp_calc_model_t &calc_model, const sp_jmatrix_t &jmatrix, well_t *well, connection_list_t &list_)
    : cell_data_ (calc_model->data)
    , main_vars_ (calc_model->main_variable)
    , pressure_ (calc_model->pressure)
    , gas_oil_ratio_ (calc_model->gas_oil_ratio)
    , invers_fvf_average_ (calc_model->invers_fvf_average)
    , well_rate_ (well->rate_)
    , well_rate_rc_ (well->rate_rc_)
    , gor_ (well->gor_)
    , d_w (0)
    , d_g (0)
    , d_o (0)
    , list_ (list_)
    , injection_type_ (well->get_well_controller ()->injection ())
    , control_type_ (well->get_well_controller ()->get_control_type ())
    , sp_diagonal_ (jmatrix->get_sp_diagonal ())
    , sec_rhs_ (jmatrix->get_sec_rhs ())
    , n_sec_vars_ (calc_model->n_sec_vars)
    , ww (well->ww_value)
    , wefac (well->exploitation_factor_)
    {
      if (is_w) d_w = calc_model->phase_d[FI_PHASE_WATER];
      if (is_g) d_g = calc_model->phase_d[FI_PHASE_GAS];
      if (is_o) d_o = calc_model->phase_d[FI_PHASE_OIL];
    }

    const cell_data_t         &cell_data_;
    const main_vars_t         main_vars_;
    const item_array_t        &pressure_;
    const item_array_t        &gas_oil_ratio_;
    const typename calc_model_t::invers_fvf_avgerage_t &invers_fvf_average_;
    rate_data_t               &well_rate_;
    rate_data_t               &well_rate_rc_;
    item_t                    &gor_;
    index_t                   d_w;
    index_t                   d_g;
    index_t                   d_o;
    bool                      is_production;
    connection_list_t         &list_;
    wells::injection_type     injection_type_;
    wells::rate_control_type  control_type_;

    const rhs_item_array_t    &sp_diagonal_;
    const rhs_item_array_t    &sec_rhs_;
    int                       n_sec_vars_;

    item_t                    krow_tetaow;
    item_t                    krp_tetap;
    item_t                    sw_part;
    item_t                    po_part;
    item_t                    sg_part;

    item_t                    &ww;
    item_t                    wefac;

    void
    calc_rate (well_t *well)
    {
      connection_loop <true, true> (well);
      if (is_production_well)
        {
          compute_bw (well->bw_value, well->rate_.prod, well->get_well_controller ()->rate ().prod);
        }
      else
        {
          compute_bw (well->bw_value, well->rate_.inj, well->get_well_controller ()->rate ().inj);
        }
    }

    void
    calc_derivs (well_t *well, bool is_rate)
    {
      if (is_rate)
        {
          connection_loop <false, true> (well);
          update_loop (well);
        }
      else
        {
          connection_loop <false, false> (well);
        }
    }

  private:

    template <bool is_rate_loop, bool is_rate>
    void 
    connection_loop (well_t *well)
    {
      if (list_.empty ())
        {
          BOSOUT (section::wells, level::debug)
            << "[" << well->name () << "] connection loop: connection list is empty"
            << bs_end;

          return ;
        }

      for (size_t i = 0, cnt = list_.size (); i < cnt; ++i)
        {
          typename well_t::sp_default_connection_t &c (list_[i]);
          if (c->is_shut ())
            {
              continue;
            }

          index_t n_block = c->n_block ();
          item_t Po  = c->cur_bhp - pressure_[n_block];
          if (Po < item_t (0.0))
            {
              connection_loop_main_var <is_rate_loop, is_rate, true> (c, cell_data_[n_block], main_vars_[n_block], Po);
            }
          else
            {
              connection_loop_main_var <is_rate_loop, is_rate, false> (c, cell_data_[n_block], main_vars_[n_block], Po);
            }
        }
    }

    template <bool is_rate_loop, bool is_rate, bool is_prod>
    void
    connection_loop_main_var (connection_t *c, const calc_model_data_t &data, main_var_t main_var, item_t Po)
    {
      if (main_var == FI_SG_VAR)
        connection_loop_injection <is_rate_loop, is_rate, is_prod, FI_SG_VAR> (c, data, Po);
      else if (main_var == FI_RO_VAR)
        connection_loop_injection <is_rate_loop, is_rate, is_prod, FI_RO_VAR> (c, data, Po);
      else if (main_var == FI_MOMG_VAR)
        connection_loop_injection <is_rate_loop, is_rate, is_prod, FI_MOMG_VAR> (c, data, Po);
    }

    template <bool is_rate_loop, bool is_rate, bool is_prod, main_var_t main_var>
    void
    connection_loop_injection (connection_t *c, const calc_model_data_t &data, item_t Po)
    {
      if (is_water_injection <is_prod> () && is_gas_injection <is_prod> () && is_oil_injection <is_prod> ())
        {
          connection_loop_concrete <is_rate_loop, is_rate, is_prod, main_var, true, true, true> (c, data, Po);
        }
      else if (is_water_injection <is_prod> () && is_oil_injection <is_prod> ())
        {
          connection_loop_concrete <is_rate_loop, is_rate, is_prod, main_var, true, false, true> (c, data, Po);
        }
      else if (is_gas_injection <is_prod> () && is_oil_injection <is_prod> ())
        {
          connection_loop_concrete <is_rate_loop, is_rate, is_prod, main_var, false, true, true> (c, data, Po);
        }
      else if (is_water_injection <is_prod> ())
        {
          connection_loop_concrete <is_rate_loop, is_rate, is_prod, main_var, true, false, false> (c, data, Po);
        }
      else if (is_gas_injection <is_prod> ())
        {
          connection_loop_concrete <is_rate_loop, is_rate, is_prod, main_var, false, true, false> (c, data, Po);
        }
      else if (is_oil_injection <is_prod> ())
        {
          connection_loop_concrete <is_rate_loop, is_rate, is_prod, main_var, false, false, true> (c, data, Po);
        }
    }

    template <bool is_rate_loop, bool is_rate, bool is_prod, main_var_t main_var, bool is_w_inj, bool is_g_inj, bool is_o_inj>
    void
    connection_loop_concrete (connection_t *c, const calc_model_data_t &data, item_t Po)
    {
      if (is_rate_loop)
        {
          calc_rate <is_prod, main_var, is_w_inj, is_g_inj, is_o_inj> (c, data, Po);
        }
      else
        {
          calc_derivs <is_prod, is_rate, main_var, is_w_inj, is_g_inj, is_o_inj> (c, data, Po);
        }
    }


    void
    update_loop (well_t *well)
    {
      if (list_.empty ())
        {
          BOSOUT (section::wells, level::debug)
            << "[" << well->name () << "] connection loop: connection list is empty"
            << bs_end;

          return ;
        }

      item_t inv_ww = 1.0;
      if (fabs (ww) >= 10e-16)
        {
          inv_ww = 1.0 / ww;
        }
      item_t ww_bw = well->bw_value * inv_ww;

      for (size_t i = 0, cnt = list_.size (); i < cnt; ++i)
        {
          typename well_t::sp_default_connection_t &c (list_[i]);
          if (c->is_shut ())
            {
              continue;
            }

          update_wr (c, inv_ww); 
          update_rate (c, ww_bw);
        }
    }

    template <bool is_prod, main_var_t main_var, bool is_w_inj, bool is_g_inj, bool is_o_inj>
    void
    calc_rate (connection_t *c, const calc_model_data_t &data, item_t Po)
    {
      index_t n_block = c->n_block ();
      item_t gw  = c->get_fact ();
      item_t Pw  = Po;
      item_t Pg  = Po;
      item_t wat = 0, gas = 0, oil = 0, free_gas = 0, solution_gas = 0;

      if (is_w_inj)
        Pw -= data.cap_pressure [d_w];
      if (is_g_inj)
        Pg -= data.cap_pressure[d_g];

      if (!is_prod)
        {
          krow_tetaow = 0;
          krp_tetap   = 0;

          if (is_w)         krow_tetaow += RELATIVE_PERM_W * INVERS_VISC_W;
          if (is_o)         krow_tetaow += RELATIVE_PERM_O * INVERS_VISC_O;
                            krp_tetap   += krow_tetaow;
          if (is_g)         krp_tetap   += RELATIVE_PERM_G * INVERS_VISC_G;
        }

      if (is_g_inj) 
        {
          if (is_prod)
            {
              free_gas                  = data.mobility[d_g];
              solution_gas              = data.mobility[d_o] * gas_oil_ratio_[n_block];

              if (main_var == FI_SG_VAR)
                {
                  c->mobility_value[gas_idx] = free_gas + solution_gas;
                }
              else
                {
                  c->mobility_value[gas_idx] = solution_gas;
                }

              free_gas                  = gw * Pg * free_gas;
              solution_gas              = gw * Po * solution_gas;
              gas                       = free_gas + solution_gas;
            }
          else
            {
              gas                       = detail::mobility <is_prod, main_var>::gas (data, d_g, this);
              c->mobility_value[gas_idx]= gas;
              gas                       = gw * Pg * gas;
            }
        }
      if (is_o_inj) 
        {
          oil                           = detail::mobility <is_prod, main_var>::oil (data, d_o, this);
          c->mobility_value[oil_idx]    = oil;
          oil                           = gw * Po * oil;
        }
      if (is_w_inj) 
        {
          wat                           = detail::mobility <is_prod, main_var>::water (data, d_w, this);
          c->mobility_value[wat_idx]    = wat;
          wat                           = gw * Pw * wat;
        }

      if (is_g_inj) c->rate_value[gas_idx]      = gas;
      if (is_o_inj) c->rate_value[oil_idx]      = oil;
      if (is_w_inj) c->rate_value[wat_idx]      = wat;

      rate_data_inner_t &rate                   = detail::rate_data <is_prod>::get (c->rate_);
      if (is_g_inj) rate.gas                    = gas;
      if (is_o_inj) rate.oil                    = oil;
      if (is_w_inj) rate.water                  = wat;
      if (is_w_inj || is_o_inj) rate.liquid     = wat + oil;

      if (is_prod && is_g_inj)
        {
          c->rate_.free_gas                     = free_gas;
          c->rate_.solution_gas                 = solution_gas;
        }

      rate_data_inner_t &rate_rc                = detail::rate_data <is_prod>::get (c->rate_rc_);
      if (is_prod && is_g_inj) 
        {
          rate_rc.gas                           = c->rate_.free_gas + c->rate_.solution_gas / data.invers_fvf[d_g];
          c->rate_rc_.free_gas                  = c->rate_.free_gas / invers_fvf_average_ [d_g];
          c->rate_rc_.solution_gas              = c->rate_.solution_gas;
        }
      else if (!is_prod && is_g_inj)
        {
          rate_rc.gas                           = c->rate_.free_gas / invers_fvf_average_ [d_g];
        }
      if (is_w_inj) rate_rc.water               = wat / invers_fvf_average_ [d_w];
      if (is_o_inj) rate_rc.oil                 = oil / invers_fvf_average_ [d_o];
      if (is_w_inj || is_o_inj) rate_rc.liquid  = rate_rc.water + rate_rc.oil;

      rate_data_inner_t &well_rate              = detail::rate_data <is_prod>::get (well_rate_);
      if (is_g_inj) well_rate.gas               += gas;
      if (is_o_inj) well_rate.oil               += oil;
      if (is_w_inj) well_rate.water             += wat;
      if (is_w_inj || is_o_inj) 
        well_rate.liquid                        += wat + oil;

      if (is_prod)
        {
          well_rate_.free_gas                   += c->rate_.free_gas;
          well_rate_.solution_gas               += c->rate_.solution_gas;
        }

      rate_data_inner_t &well_rate_rc           = detail::rate_data <is_prod>::get (well_rate_rc_);
      if (is_g_inj) well_rate_rc.gas            += rate_rc.gas;
      if (is_o_inj) well_rate_rc.oil            += rate_rc.oil;
      if (is_w_inj) well_rate_rc.water          += rate_rc.water;
      if (is_w_inj || is_o_inj)
        well_rate_rc.liquid                     += rate_rc.water + rate_rc.oil;

      if (is_prod)
        {
          well_rate_rc_.free_gas                += c->rate_rc_.free_gas;
          gor_                                  = well_rate_.solution_gas / rate.oil;
        }
    }

    template <bool is_prod, main_var_t main_var>
    void
    wat_deriv (connection_t *c, const calc_model_data_t &data, const boost::array <item_t, FI_PHASE_TOT> &mobility, item_t gw, item_t Pw)
    {
      if (is_w)
        {
          item_t sw_deriv = 0, po_deriv = 0, so_deriv = 0, sg_deriv = 0;

          if (is_prod)
            {
              if (is_g)                           sg_deriv = data.s_deriv_mobility[d_w * n_phases + d_g];
              if (is_o)                           so_deriv = data.s_deriv_mobility[d_w * n_phases + d_o];
                                                  po_deriv = data.p_deriv_mobility[d_w];
              if (!is_1p)                         sw_deriv = data.s_deriv_mobility[d_w * n_phases + d_w];
            }
          else
            {
              if (is_g && main_var == FI_SG_VAR)  sg_deriv = INVERS_FVF_W * sg_part;
              if (is_g && main_var != FI_SG_VAR)  sg_deriv = INVERS_FVF_W * RELATIVE_PERM_O * GOR_DERIV_INVERS_VISC;
              if (is_o)                           so_deriv = S_DERIV_RELATIVE_PERM_OO * INVERS_VISC_O;
              item_t                              common   = INVERS_FVF_W * krp_tetap * P_DERIV_INVERS_FVF_W;
                                                  po_deriv = (RELATIVE_PERM_O * P_DERIV_INVERS_VISC_O 
                                                              + RELATIVE_PERM_W * P_DERIV_INVERS_VISC_W 
                                                              - common);
              if (!is_1p)                         sw_deriv = (S_DERIV_RELATIVE_PERM_WW * INVERS_VISC_W 
                                                              + S_DERIV_RELATIVE_PERM_OW * INVERS_VISC_O 
                                                              + RELATIVE_PERM_W * P_DERIV_INVERS_VISC_W * S_DERIV_CAP_PRESSURE_W 
                                                              - common);
            }


          if (is_g)   c->rr_value[wat_sg] = -gw * (Pw * sg_deriv);
          if (is_o)   c->rr_value[wat_so] = -gw * (Pw * so_deriv);
          if (is_w)   c->rr_value[wat_po] = -gw * (Pw * po_deriv - mobility[wat_idx]);
          if (!is_1p) c->rr_value[wat_sw] = -gw * (Pw * sw_deriv - mobility[wat_idx] * S_DERIV_CAP_PRESSURE_W);
        }
    }

    template <bool is_prod, main_var_t main_var>
    void
    gas_deriv (connection_t *c, const calc_model_data_t &data, const boost::array <item_t, FI_PHASE_TOT> &mobility, item_t gw, item_t Pg, item_t Po, index_t n_block)
    {
      if (is_g)
        {
          item_t sw_deriv = 0, po_deriv = 0, so_deriv = 0, sg_deriv = 0;
          if (is_prod)
            {
              if (!is_1p && main_var == FI_SG_VAR)  sg_deriv = Pg * data.s_deriv_mobility[d_g * n_phases + d_g] 
                                                                - data.mobility[d_g] * data.s_deriv_cap_pressure[d_g] 
                                                                + gas_oil_ratio_[n_block] * data.s_deriv_mobility[d_o * n_phases + d_g] * Po;
              if (!is_1p && main_var != FI_SG_VAR)  sg_deriv = Po * (data.mobility[d_o] + gas_oil_ratio_[n_block] * data.s_deriv_mobility[d_o * n_phases + d_g]);
              if (is_o)                             so_deriv = Pg * data.s_deriv_mobility[d_g * n_phases + d_o] + Po * data.s_deriv_mobility[d_o * n_phases + d_o] * gas_oil_ratio_[n_block];
                                                    po_deriv = gas_oil_ratio_[n_block] * (data.p_deriv_mobility[d_o] * Po - data.mobility[d_o]);
              if (main_var == FI_SG_VAR)            po_deriv += data.p_deriv_mobility[d_g] * Pg - data.mobility[d_g] + data.p_deriv_gas_oil_ratio * data.mobility[d_o] * Po;
              if (is_w)                             sw_deriv = Pg * data.s_deriv_mobility[d_g * n_phases + d_w] + Po * data.s_deriv_mobility[d_o * n_phases + d_w] * gas_oil_ratio_[n_block];
            }
          else
            {
              if (!is_1p && main_var == FI_SG_VAR)  sg_deriv = Pg * (S_DERIV_INVERS_FVF_G * krp_tetap + INVERS_FVF_G * sg_part) - mobility[gas_idx] * S_DERIV_CAP_PRESSURE_G;
              if (!is_1p && main_var != FI_SG_VAR)  sg_deriv = Pg * INVERS_FVF_G * RELATIVE_PERM_O * data.gor_deriv_invers_viscosity;
              if (is_o)                             so_deriv = Pg * INVERS_FVF_G * S_DERIV_RELATIVE_PERM_GO;
              if (main_var == FI_SG_VAR)            po_deriv = Pg * (P_DERIV_INVERS_FVF_G * krp_tetap + INVERS_FVF_G * po_part) - INVERS_FVF_G * krp_tetap;
              if (main_var != FI_SG_VAR)            po_deriv = Pg * (P_DERIV_INVERS_FVF_G * krp_tetap + INVERS_FVF_G * po_part) - INVERS_FVF_G * krow_tetaow;
              if (is_w)                             sw_deriv = Pg * INVERS_FVF_G * sw_part;
            }

          if (!is_1p) c->rr_value[gas_sg] = -gw * sg_deriv;
          if (is_o)   c->rr_value[gas_so] = -gw * so_deriv;
          if (is_g)   c->rr_value[gas_po] = -gw * po_deriv;
          if (is_w)   c->rr_value[gas_sw] = -gw * sw_deriv;
        }
    }

    template <bool is_prod, main_var_t main_var>
    void
    oil_deriv (connection_t *c, const calc_model_data_t &data, const boost::array <item_t, FI_PHASE_TOT> &mobility, item_t gw, item_t Po)
    {
      if (is_o)
        {
          item_t sw_deriv = 0, po_deriv = 0, so_deriv = 0, sg_deriv = 0;

          if (is_prod)
            {
              if (is_g && main_var == FI_SG_VAR)  sg_deriv = data.s_deriv_mobility[d_o * n_phases + d_g];
              if (is_g && main_var != FI_SG_VAR)  sg_deriv = data.relative_perm[d_o] * data.gor_deriv_invers_visc_fvf;
              if (!is_1p)                         so_deriv = data.s_deriv_mobility[d_o * n_phases + d_o];
                                                  po_deriv = data.p_deriv_mobility[d_o];
              if (is_w)                           sw_deriv = data.s_deriv_mobility[d_o * n_phases + d_w];
            }
          else
            {
              if (is_g && main_var == FI_SG_VAR)  sg_deriv = INVERS_FVF_O * sg_part;
              if (is_g && main_var != FI_SG_VAR)  sg_deriv = RELATIVE_PERM_O * (INVERS_FVF_O * data.gor_deriv_invers_viscosity + INVERS_VISC_O * data.gor_deriv_invers_fvf);
              if (!is_1p)                         so_deriv = INVERS_FVF_O * S_DERIV_RELATIVE_PERM_OO;
                                                  po_deriv = P_DERIV_INVERS_FVF_O * krp_tetap + INVERS_FVF_O * po_part;
              if (is_w)                           sw_deriv = INVERS_FVF_O * sw_part;
            }

          if (is_g)   c->rr_value[oil_sg] = -gw * (Po * sg_deriv);
          if (!is_1p) c->rr_value[oil_so] = -gw * (Po * so_deriv);
          if (is_o)   c->rr_value[oil_po] = -gw * (Po * po_deriv - mobility[oil_idx]);
          if (is_w)   c->rr_value[oil_sw] = -gw * (Po * sw_deriv);
        }
    }

    template <bool is_prod, bool is_rate, main_var_t main_var, bool is_w_inj, bool is_g_inj, bool is_o_inj>
    void
    calc_derivs (connection_t *c, const calc_model_data_t &data, item_t Po)
    {
      index_t n_block = c->n_block ();
      item_t gw  = c->get_fact ();
      item_t Pw  = Po;
      item_t Pg  = Po;

      if (is_w)
        Pw -= data.cap_pressure [d_w];
      if (is_g)
        Pg -= data.cap_pressure[d_g];

      if (!is_prod)
        {
          sw_part = 0;
          po_part = 0;
          sg_part = 0;

          if (is_g)         sg_part += S_DERIV_RELATIVE_PERM_GG * INVERS_VISC_G;
          if (is_g)         sg_part += RELATIVE_PERM_G * S_DERIV_INVERS_VISC_G * CAP_PRESSURE_G;
          if (is_g && is_o) sg_part += S_DERIV_RELATIVE_PERM_OG * INVERS_VISC_O;

          if (is_g)         po_part += RELATIVE_PERM_G * P_DERIV_INVERS_VISC_G;
          if (is_o)         po_part += RELATIVE_PERM_O * P_DERIV_INVERS_VISC_O;
          if (is_w)         po_part += RELATIVE_PERM_W * P_DERIV_INVERS_VISC_W;

          if (is_w && is_o) sw_part += S_DERIV_RELATIVE_PERM_OW * INVERS_VISC_O;
          if (is_w)         sw_part += S_DERIV_RELATIVE_PERM_WW * INVERS_VISC_W;
          if (is_w)         sw_part += RELATIVE_PERM_W * S_DERIV_INVERS_VISC_W * CAP_PRESSURE_W;
          if (is_w)         sw_part -= P_DERIV_INVERS_FVF_W * CAP_PRESSURE_W;


        }

      if (is_g && is_g_inj) gas_deriv <is_prod, main_var> (c, data, c->mobility_value, gw, Pg, Po, n_block);
      if (is_o && is_o_inj) oil_deriv <is_prod, main_var> (c, data, c->mobility_value, gw, Po);
      if (is_w && is_w_inj) wat_deriv <is_prod, main_var> (c, data, c->mobility_value, gw, Pw);

      m_minus_vv_prod <n_phases>::eliminate (&c->rr_value[b_sqr], &sp_diagonal_[n_block * n_phases], c->rr_value);

      if (is_rate)
        {
          if (is_w && is_water_control <is_prod> ())
            {
              if (n_phases > 0) c->wr_value[0]          += c->rr_value [wat_idx * n_phases + 0];
              if (n_phases > 1) c->wr_value[1]          += c->rr_value [wat_idx * n_phases + 1];
              if (n_phases > 2) c->wr_value[2]          += c->rr_value [wat_idx * n_phases + 2];
            }
          if (is_g && is_gas_control <is_prod> ())
            {
              if (n_phases > 0) c->wr_value[0]          += c->rr_value [gas_idx * n_phases + 0];
              if (n_phases > 1) c->wr_value[1]          += c->rr_value [gas_idx * n_phases + 1];
              if (n_phases > 2) c->wr_value[2]          += c->rr_value [gas_idx * n_phases + 2];
            }
          if (is_o && is_oil_control <is_prod> ())
            {
              if (n_phases > 0) c->wr_value[0]          += c->rr_value [oil_idx * n_phases + 0];
              if (n_phases > 1) c->wr_value[1]          += c->rr_value [oil_idx * n_phases + 1];
              if (n_phases > 2) c->wr_value[2]          += c->rr_value [oil_idx * n_phases + 2];
            }

          item_t rw_wat = 0, rw_gas = 0, rw_oil = 0;

          if (is_w && is_water_injection <is_prod> ())  rw_wat = c->rw_value[wat_idx]    = -gw * c->mobility_value[wat_idx];
          if (is_g && is_gas_injection <is_prod> ())    rw_gas = c->rw_value[gas_idx]    = -gw * c->mobility_value[gas_idx];
          if (is_o && is_oil_injection <is_prod> ())    rw_oil = c->rw_value[oil_idx]    = -gw * c->mobility_value[oil_idx];

          if (is_w && is_water_control <is_prod> ())    ww += rw_wat;
          if (is_g && is_gas_control <is_prod> ())      ww += rw_gas;
          if (is_o && is_oil_control <is_prod> ())      ww += rw_oil;
        }

      if (wefac > 0.0)
        {
          if (n_phases > 0) c->rr_value[0] = wefac * c->rr_value[0];
          if (n_phases > 1) c->rr_value[1] = wefac * c->rr_value[1];
          if (n_phases > 1) c->rr_value[2] = wefac * c->rr_value[2];
          if (n_phases > 1) c->rr_value[3] = wefac * c->rr_value[3];
          if (n_phases > 2) c->rr_value[4] = wefac * c->rr_value[4];
          if (n_phases > 2) c->rr_value[5] = wefac * c->rr_value[5];
          if (n_phases > 2) c->rr_value[6] = wefac * c->rr_value[6];
          if (n_phases > 2) c->rr_value[7] = wefac * c->rr_value[7];
          if (n_phases > 2) c->rr_value[8] = wefac * c->rr_value[8];

          if (is_rate)
            {
              if (n_phases > 0) c->rw_value[0] = wefac * c->rw_value[0];
              if (n_phases > 1) c->rw_value[1] = wefac * c->rw_value[1];
              if (n_phases > 2) c->rw_value[2] = wefac * c->rw_value[2];
            }
        }

      v_minus_vs_prod <n_phases>::eliminate (&c->rr_value[b_sqr], sec_rhs_[n_sec_vars_ * n_block], c->rate_value);
      //c->ps_value[0] = c->rr_value[b_sqr + 0];
      //c->ps_value[1] = c->rr_value[b_sqr + 1];
      //c->ps_value[2] = c->rr_value[b_sqr + 2];
      //c->rr_value[b_sqr + 0] = 0;
      //c->rr_value[b_sqr + 1] = 0;
      //c->rr_value[b_sqr + 2] = 0;
    }

    void
    compute_bw (item_t &bw, const rate_data_inner_t &rate, const rate_data_inner_t &limit)
    {
      bw = 0;
      item_t mult = is_production_well ? -1.0 : 1.0;
      if (is_o && is_w && ((is_oil_control <is_production_well> () && is_water_control <is_production_well> ()) || is_liquid_control <is_production_well> ()))
        {
          bw += (rate.oil + rate.water) - mult * limit.liquid;
        }
      else if (is_o && (is_oil_control <is_production_well> () || is_liquid_control <is_production_well> ()))
        {
          bw += (rate.oil) - mult * limit.oil;
        }
      else if (is_w && (is_water_control <is_production_well> () || is_liquid_control <is_production_well> ()))
        {
          bw += (rate.water) - mult * limit.water;
        }
      else if (is_g && is_gas_control <is_production_well> ())
        {
          bw += (rate.gas) - mult * limit.gas;
        }
    }

    void
    update_wr (connection_t *c, item_t inv_ww)
    {
      if (n_phases > 0) c->wr_value[0] *= inv_ww;
      if (n_phases > 1) c->wr_value[1] *= inv_ww;
      if (n_phases > 2) c->wr_value[2] *= inv_ww;
    }

    void
    update_rate (connection_t *c, item_t ww_bw)
    {
      if (n_phases > 0) c->rate_value[0] -= c->rw_value[0] * ww_bw;
      if (n_phases > 1) c->rate_value[1] -= c->rw_value[1] * ww_bw;
      if (n_phases > 2) c->rate_value[2] -= c->rw_value[2] * ww_bw;
    }

    template <bool is_prod>
    bool 
    is_liquid_control () const
    {
      if (is_prod)
        {
          return control_type_ == wells::liquid_rate_control;
        }
      else
        {
          return injection_type_ == wells::injection_oil || injection_type_ == wells::injection_water;
        }
    }
    template <bool is_prod>
    bool
    is_oil_control () const
    {
      if (is_prod)
        {
          return is_o && (control_type_ == wells::oil_rate_control || control_type_ == wells::liquid_rate_control);
        }
      else
        {
          return is_o && injection_type_ == wells::injection_oil;
        }
    }
    template <bool is_prod>
    bool
    is_water_control () const
    {
      if (is_prod)
        {
          return is_w && (control_type_ == wells::water_rate_control || control_type_ == wells::liquid_rate_control);
        }
      else
        {
          return is_w && injection_type_ == wells::injection_water;
        }
    }
    template <bool is_prod>
    bool
    is_gas_control () const
    {
      if (is_prod)
        {
          return is_g && control_type_ == wells::gas_rate_control;
        }
      else
        {
          return is_g && injection_type_ == wells::injection_gas;
        }
    }
    template <bool is_prod>
    bool
    is_water_injection () const
    {
      if (is_prod)
        {
          return is_w;
        }
      else
        {
          return is_w && injection_type_ == wells::injection_water;
        }
    }
    template <bool is_prod>
    bool
    is_gas_injection () const
    {
      if (is_prod)
        {
          return is_g;
        }
      else
        {
          return is_g && injection_type_ == wells::injection_gas;
        }
    }
    template <bool is_prod>
    bool
    is_oil_injection () const
    {
      if (is_prod)
        {
          return is_o;
        }
      else
        {
          return is_o && injection_type_ == wells::injection_oil;
        }
    }
  };

} // namespace wells
} // namespace blue_sky


#endif // #ifndef BS_BOS_CORE_DEFAULT_WELL_CALC_RATE_H_

