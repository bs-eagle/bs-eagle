/**
 *       \file  fi_operator_block_connections_mpfa.h
 *      \brief  Old (two point) mpfa implementation
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  06.04.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Obsolete, should be removed
 * */
#ifndef BS_FI_OPERATOR_BLOCK_CONNECTIONS_MPFA_H_
#define BS_FI_OPERATOR_BLOCK_CONNECTIONS_MPFA_H_

namespace blue_sky {

  namespace tpfa 
  {
    template <typename strategy_t>
    struct mpfa_impl 
    {
    public:
      typedef typename strategy_t::item_t             item_t;
      typedef typename strategy_t::index_t            index_t;
      typedef typename strategy_t::item_array_t       item_array_t;
      typedef typename strategy_t::index_array_t      index_array_t;
      typedef boost::array <index_t, FI_PHASE_TOT>    up_cell_array_t;
      typedef calc_model <strategy_t>                 calc_model_t;
      typedef typename calc_model_t::data_t           data_t;
      typedef typename calc_model_t::data_array_t     data_array_t;
      typedef typename calc_model_t::main_var_array_t main_var_array_t;

    public:
      mpfa_impl (const fi_operator_impl <strategy_t> &fi_operator_, item_array_t &rhs_, item_array_t &reg_values_)
      : fi_operator_ (fi_operator_),
      n_phases (fi_operator_.n_phases),
      is_w (fi_operator_.is_w),
      is_g (fi_operator_.is_g),
      is_o (fi_operator_.is_o),
      d_w (fi_operator_.d_w),
      d_g (fi_operator_.d_g),
      d_o (fi_operator_.d_o),
      saturation_3p_ (fi_operator_.saturation_3p_),
      pressure_ (fi_operator_.pressure_),
      gas_oil_ratio_ (fi_operator_.gas_oil_ratio_),
      main_vars_ (fi_operator_.main_vars_),
      data_ (fi_operator_.data_),
      gravity_ (fi_operator_.calc_model_->internal_constants.gravity_constant),
      n_sec_vars (fi_operator_.n_sec_vars),
      rhs_ (rhs_),
      reg_values_ (reg_values_)
      {
      }

      /*BS_FORCE_INLINE*/ void
      mpfa_calc_avarage_density(item_t *ave_density,
          const index_t *cell_ind_block,
          const int &n_cells_in_conn) const
      {
        item_t counter[FI_PHASE_TOT] = {0};
        index_t flag[FI_PHASE_TOT] = {0};

        memset (ave_density, 0, sizeof (item_t) * FI_PHASE_TOT);
        //ave_density[0] = ave_density[1] = ave_density[2] = 0;

        for (int j = 0; j < n_cells_in_conn; ++j)
          {
            index_t ic = cell_ind_block[j];
            index_t idx = ic * n_phases;
            const data_t &data_ic = data_[ic];

            if (n_phases > 1)
              {
                bool is_w_ = is_w && saturation_3p_[idx + d_w] > EPS_DIV;
                bool is_g_ = is_g && saturation_3p_[idx + d_g] > EPS_DIV;
                bool is_o_ = is_o && saturation_3p_[idx + d_o] > EPS_DIV;

                if (is_w_)
                  {
                    ave_density[d_w] += data_ic.density[d_w];
                    ++counter[d_w];
                    ++flag[d_w];
                  }
                if (is_g_)
                  {
                    ave_density[d_g] += data_ic.density[d_g];
                    ++counter[d_g];
                    ++flag[d_g];
                  }
                if (is_o_)
                  {
                    ave_density[d_o] += data_ic.density[d_o];
                    ++counter[d_o];
                    ++flag[d_o];
                  }
              }
            else
              {
                ave_density[0] += data_ic.density[0];
                ++counter[0];
                ++flag[0];
              }
          }
        for (int j = 0; j < n_phases; ++j)
          {
            if (flag[j])
              ave_density[j] /= counter[j];
          }
      }

      BS_FORCE_INLINE void
      mpfa_calc_potential (item_t *cell_pot,
          item_t *sum_cell_pot,
          item_t *ave_density,
          const item_t *truns_block,
          const item_array_t &depth,
          const index_t * cell_ind_block,
          index_t cell_m,
          item_t dt,
          int n_cells_in_conn) const
      {
        item_t g_h      = gravity_;
        item_t base_h   = depth[cell_m]; 
        sum_cell_pot[0] = sum_cell_pot[1] = sum_cell_pot[2] = 0;

        for (index_t j = 0; j < n_cells_in_conn; ++j)
          {
            index_t ic            = cell_ind_block[j];
            item_t dh             = depth[ic] - base_h; 
            const data_t &data_ic = data_[ic];
            item_t truns_block_dt = truns_block[j] * dt;

            if (is_w)
              {
                index_t k = j * n_phases + d_w;
                cell_pot[k] = pressure_[ic] - ave_density[d_w] * g_h * dh;
                if (this->n_phases > 1)
                  cell_pot[k] += data_ic.cap_pressure[d_w];
                sum_cell_pot[d_w] += cell_pot[k] * truns_block_dt;
              }
            if (is_g)
              {
                index_t k = j * this->n_phases + d_g;
                cell_pot[k] = pressure_[ic] - ave_density[d_g] * g_h * dh;
                if (this->n_phases > 1)
                  cell_pot[k] += data_ic.cap_pressure[d_g];
                sum_cell_pot[d_g] += cell_pot[k] * truns_block_dt;
              }
            if (is_o)
              {
                index_t k = j * this->n_phases + d_o;
                cell_pot[k] = pressure_[ic] - ave_density[d_o] * g_h * dh;
                sum_cell_pot[d_o] += cell_pot[k] * truns_block_dt;
              }
          }
      }

      /*BS_FORCE_INLINE */void
      mpfa_calc_upstream(up_cell_array_t &up_cell,
          item_t *sum_cell_pot,
          const index_t &cell_m,
          const index_t &cell_p) const
      {
        index_t cell_[] = {cell_m, cell_p};
        for (int j = 0; j < n_phases; ++j)
          {
            up_cell[j] = cell_[sum_cell_pot[j] > 0 ? 0 : 1];
          }
      }

      /*BS_FORCE_INLINE */void
      mpfa_fill_rhs(item_t * & rhs_m,
                                            item_t * & rhs_p,
                                            const up_cell_array_t &up_cell,
                                            const item_t *truns_block,
                                            item_t *cell_pot,
                                            item_array_t &flux_rhs,
                                            const item_t &dt,
                                            const int &n_cells_in_conn,
                                            const index_t &cell_m,
                                            const index_t &cell_p,
                                            const index_t &equ_w, const index_t equ_g, const index_t &equ_o) const
      {
        rhs_m = &flux_rhs[0] + cell_m * this->n_phases;
        rhs_p = &flux_rhs[0] + cell_p * this->n_phases;
        item_t rhs_block_m[3] = {0}, rhs_block_p[3] = {0};
        //rhs_block_m[0] = rhs_block_m[1] = rhs_block_m[2] = 0;
        //rhs_block_p[0] = rhs_block_p[1] = rhs_block_p[2] = 0;

        static data_t dummy_data;
        const data_t &data_m = data_[cell_m];
        const data_t &data_w = is_w ? data_[up_cell[d_w]] : dummy_data;
        const data_t &data_g = is_g ? data_[up_cell[d_g]] : dummy_data;
        const data_t &data_o = is_o ? data_[up_cell[d_o]] : dummy_data;

        item_t flow_m = data_m.truns_mult * dt;
        item_t flow_w = is_w ? flow_m * data_w.mobility[d_w] : 0;
        item_t flow_g = is_g ? flow_m * data_g.mobility[d_g] : 0;
        item_t flow_o = is_o ? flow_m * data_o.mobility[d_o] : 0;

        index_t k_j = 0;
        for (int j = 0; j < n_cells_in_conn; ++j, k_j += n_phases)
          {
            // water
            if (is_w)
              {
                index_t k = k_j + d_w;
                item_t flow = flow_w * truns_block[j] * cell_pot[k];
                rhs_block_m[equ_w] -= flow;
                rhs_block_p[equ_w] += flow;
              }
            // oil
            if (is_o)
              {
                index_t k = k_j + d_o;
                item_t flow = flow_o * truns_block[j] * cell_pot[k];
                rhs_block_m[equ_o] -= flow;
                rhs_block_p[equ_o] += flow;
                // disolved gas
                if (is_g)
                  {
                    flow *= gas_oil_ratio_[up_cell[d_o]];
                    rhs_block_m[equ_g] -= flow;
                    rhs_block_p[equ_g] += flow;
                  }
              }
            // gas
            if (is_g)
              {
                index_t k = k_j + d_g;
                item_t flow = flow_g * truns_block[j] * cell_pot[k];
                rhs_block_m[equ_g] -= flow;
                rhs_block_p[equ_g] += flow;
              }
          }
        for (int j = 0; j < n_phases; ++j)
          {
            rhs_m[j] += rhs_block_m[j];
            rhs_p[j] += rhs_block_p[j];
          }
      }

      /*BS_FORCE_INLINE*/ void
      mpfa_fill_jacobian( boost::array <item_t, 18>  &m_mem,
          boost::array <item_t, 18>  &p_mem,
          const up_cell_array_t &up_cell,
          const item_t *cell_pot,
          const item_t *truns_block,
          item_t * &rhs_m,
          item_t * &rhs_p,
          const index_t *cell_ind_block,
          item_array_t &sp_diag,
          item_array_t &s_rhs,
          const item_t &dt,
          const int &n_cells_in_conn,
          const index_t &cell_m,
          const index_t &cell_p,
          const index_t &equ_w, const index_t equ_g, const index_t &equ_o) const
      {
        //m_mem.assign(0);
        //p_mem.assign(0);
        assign (m_mem, item_t (0));
        assign (p_mem, item_t (0));

        static data_t dummy_data;

        item_t *sp_diag_val = &sp_diag[0];
        item_t *s_rhs_val = &s_rhs[0];

        const data_t *data_m = &data_[cell_m];
        const data_t *data_w = &dummy_data;
        const data_t *data_o = &dummy_data;
        const data_t *data_g = &dummy_data;

        index_t up_cell_dw  = 0;
        index_t up_cell_do  = 0;
        index_t up_cell_dg  = 0;

        index_t sp_size_w   = 0;
        index_t sp_size_o   = 0;
        index_t sp_size_g   = 0;

        index_t rhs_size_w  = 0;
        index_t rhs_size_o  = 0;
        index_t rhs_size_g  = 0;

        item_t truns_mob_w  = 0;
        item_t truns_mob_o  = 0;
        item_t truns_mob_g  = 0;

        item_t truns_p_mob_w  = 0;
        item_t truns_p_mob_o  = 0;
        item_t truns_p_mob_g  = 0;

        item_t truns_s_mob_ww  = 0;
        item_t truns_s_mob_ow  = 0;
        item_t truns_s_mob_og  = 0;
        item_t truns_s_mob_oo  = 0;
        item_t truns_s_mob_gg  = 0;

        // calculate mobility indexies
        index_t mob_w_dw = d_w * n_phases + d_w;
        index_t mob_o_dw = d_o * n_phases + d_w;
        index_t mob_o_dg = d_o * n_phases + d_g;
        index_t mob_o_do = d_o * n_phases + d_o;
        index_t mob_g_dg = d_g * n_phases + d_g;

        if (is_w)
          {
            data_w      = &data_[up_cell[d_w]];
            up_cell_dw  = up_cell[d_w];
            sp_size_w   = n_sec_vars * this->n_phases * up_cell_dw;
            rhs_size_w  = n_sec_vars * up_cell_dw;
            truns_mob_w = data_m->truns_mult * dt * data_w->mobility[d_w];
            truns_p_mob_w = data_m->truns_mult *dt * data_w->p_deriv_mobility[d_w];
            truns_s_mob_ww  = data_m->truns_mult * dt * data_w->s_deriv_mobility[mob_w_dw];
         }

        if (is_o)
          {
            data_o      = &data_[up_cell[d_o]];
            up_cell_do  = up_cell[d_o];
            rhs_size_o  = n_sec_vars * up_cell_do;
            sp_size_o   = n_sec_vars * this->n_phases * up_cell_do;
            truns_mob_o = data_m->truns_mult * dt * data_o->mobility[d_o];
            truns_p_mob_o = data_m->truns_mult *dt * data_o->p_deriv_mobility[d_o];
            truns_s_mob_ow  = data_m->truns_mult * dt * data_o->s_deriv_mobility[mob_o_dw];
            truns_s_mob_og  = data_m->truns_mult * dt * data_o->s_deriv_mobility[mob_o_dg];
            truns_s_mob_oo  = data_m->truns_mult * dt * data_o->s_deriv_mobility[mob_o_do];
          }
        if (is_g)
          {
            data_g      = &data_[up_cell[d_g]];
            up_cell_dg  = up_cell[d_g];
            sp_size_g   = n_sec_vars * this->n_phases * up_cell_dg;
            rhs_size_g  = n_sec_vars * up_cell_dg;
            truns_mob_g = data_m->truns_mult * dt * data_g->mobility[d_g];
            truns_p_mob_g = data_m->truns_mult *dt * data_g->p_deriv_mobility[d_g];
            truns_s_mob_gg  = data_m->truns_mult * data_g->s_deriv_mobility[mob_g_dg];
         }

        index_t k_j = 0;
        for (int j = 0; j < n_cells_in_conn; ++j, k_j += n_phases)
          {
            index_t ic = cell_ind_block[j];

            index_t k;
            item_t dsw, dsg, dso, dpo;
            item_t dsw_up, dsg_up, dpo_up, dso_up;

            item_t *pp_block_m = &m_mem [j * 9];
            item_t *pp_block_p = &p_mem [j * 9];
            item_t *pp_block_m_up;
            item_t *pp_block_p_up;
            item_t ps_block_p_up[3];
            item_t ps_block_m_up[3];
            item_t *sp_block_up;
            item_t *s_rhs_block_up;

            item_t ps_block_m[FI_PHASE_TOT] = {0};   // temporary block for PS
            item_t ps_block_p[FI_PHASE_TOT] = {0};   // temporary block for PS

            item_t *sp_block = sp_diag_val + n_sec_vars * this->n_phases * ic;
            item_t *s_rhs_block = s_rhs_val + n_sec_vars * ic;

            item_t truns_block_j = truns_block[j];

            // water
            if (is_w)
              {
                ps_block_p_up[0]  = ps_block_p_up[1] = ps_block_p_up[2] = 0;
                ps_block_m_up[0]  = ps_block_m_up[1] = ps_block_m_up[2] = 0;
                sp_block_up       = sp_diag_val + sp_size_w;
                s_rhs_block_up    = s_rhs_val + rhs_size_w;

                if (up_cell_dw == cell_ind_block[0])
                  {
                    pp_block_m_up = &m_mem[0];
                    pp_block_p_up = &p_mem[0];
                  }
                else
                  {
                    pp_block_m_up = &m_mem[9];
                    pp_block_p_up = &p_mem[9];
                  }
                k = k_j + d_w;
                item_t cell_pot_k = cell_pot[k];
                dpo = truns_mob_w * truns_block[j];
                dpo_up = 0;
                if (ic == up_cell_dw)
                  {
                    dpo += truns_p_mob_w * truns_block_j * cell_pot_k;
                  }
                else
                  {
                    dpo_up = truns_p_mob_w * truns_block_j * cell_pot_k;
                  }

                if (n_phases > 1)
                  {
                    dsg = 0;
                    dso = 0;
                    dsw = truns_mob_w * truns_block_j * data_[ic].s_deriv_cap_pressure[d_w];
                    dsg_up = dso_up = dsw_up = 0;
                    if (ic == up_cell_dw)
                      {
                        dsw += truns_s_mob_ww * truns_block_j * cell_pot_k;
                      }
                    else
                      {
                        dsw_up = truns_s_mob_ww * truns_block_j * cell_pot_k;
                      }
                    if (n_phases == 3)
                      {
                        pp_block_m[p3_wat_po]     += dpo;
                        pp_block_m_up[p3_wat_po]  += dpo_up;
                        pp_block_p[p3_wat_po]     -= dpo;
                        pp_block_p_up[p3_wat_po]  -= dpo_up;

                        ps_block_m[p3_wat]        += dsw;
                        ps_block_m_up[p3_wat]     += dsw_up;
                        ps_block_p[p3_wat]        += -dsw;
                        ps_block_p_up[p3_wat]     += -dsw_up;
                      }
                    else
                      {
                        pp_block_m[p2ow_wat_po]     += dpo;
                        pp_block_m_up[p2ow_wat_po]  += dpo_up;
                        pp_block_p[p2ow_wat_po]     -= dpo;
                        pp_block_p_up[p2ow_wat_po]  -= dpo_up;

                        ps_block_m[p2ow_wat]        += dsw;
                        ps_block_m_up[p2ow_wat]     += dsw_up;
                        ps_block_p[p2ow_wat]        += -dsw;
                        ps_block_p_up[p2ow_wat]     += -dsw_up;
                      }

                    M_MINUS_VV_PROD (n_phases, ps_block_m_up, sp_block_up, pp_block_m_up);
                    V_MINUS_VS_PROD (n_phases, ps_block_m_up, s_rhs_block_up, rhs_m);

                    M_MINUS_VV_PROD (n_phases, ps_block_p_up, sp_block_up, pp_block_p_up);
                    V_MINUS_VS_PROD (n_phases, ps_block_p_up, s_rhs_block_up, rhs_p);
                  }
                else
                  {
                    pp_block_m[0]     += dpo;
                    pp_block_m_up[0]  += dpo_up;
                    pp_block_p[0]     -= dpo;
                    pp_block_p_up[0]  -= dpo_up;
                  }
              }
            // oil
            if (is_o)
              {
                item_t g_dpo, g_dso, g_dsw, g_dsg;
                item_t g_dpo_up, g_dso_up, g_dsw_up, g_dsg_up;
                k = k_j + d_o;
                item_t cell_pot_k = cell_pot[k];

                ps_block_p_up[0] = ps_block_p_up[1] = ps_block_p_up[2] = 0;
                ps_block_m_up[0] = ps_block_m_up[1] = ps_block_m_up[2] = 0;
                sp_block_up = sp_diag_val + sp_size_o;
                s_rhs_block_up = s_rhs_val + rhs_size_o;

                if (up_cell_do == cell_ind_block[0])
                  {
                    pp_block_m_up = &m_mem[0];
                    pp_block_p_up = &p_mem[0];
                  }
                else
                  {
                    pp_block_m_up = &m_mem[9];
                    pp_block_p_up = &p_mem[9];
                  }

                // calculate pressure derivatives
                dpo = truns_mob_o * truns_block_j;
                dpo_up = 0;
                if (ic == up_cell_do)
                  {
                    dpo += truns_p_mob_o * truns_block_j * cell_pot_k;
                  }
                else
                  {
                    dpo_up = truns_p_mob_o * truns_block_j * cell_pot_k;
                  }
                if (is_g)
                  {
                    g_dpo = gas_oil_ratio_[up_cell_do] * dpo;
                    g_dpo_up = gas_oil_ratio_[up_cell_do] * dpo_up;
                    if (ic == up_cell_do && FI_CHK_SG (main_vars_, ic))
                      {
                        g_dpo += truns_mob_o * truns_block_j * data_o->p_deriv_gas_oil_ratio * cell_pot_k;
                      }
                    else if (FI_CHK_SG (main_vars_, up_cell_do))
                      {
                        g_dpo_up += truns_mob_o * truns_block_j * data_o->p_deriv_gas_oil_ratio * cell_pot_k;
                      }
                  }

                dsw = 0;
                dso = 0;
                dsg = 0;
                dsw_up = 0;
                dso_up = 0;
                dsg_up = 0;

                // calculate saturation deriv
                if (ic == up_cell[d_o])
                  {
                    if (is_w)
                      dsw = truns_s_mob_ow * truns_block_j * cell_pot_k;
                    if (is_g)
                      dsg = truns_s_mob_og * truns_block_j * cell_pot_k;
                    if (n_phases > 1)
                      dso = truns_s_mob_oo * truns_block_j * cell_pot_k;
                  }
                else
                  {
                    if (is_w)
                      dsw_up = truns_s_mob_ow * truns_block_j * cell_pot_k;
                    if (is_g)
                      dsg_up = truns_s_mob_og * truns_block_j * cell_pot_k;
                    if (n_phases > 1)
                      dso_up = truns_s_mob_oo * truns_block_j * cell_pot_k;
                  }
                if (is_g)
                  {
                    g_dso = gas_oil_ratio_[up_cell_do] * dso;
                    g_dso_up = gas_oil_ratio_[up_cell_do] * dso_up;
                    if (is_w)
                      {
                        g_dsw = gas_oil_ratio_[up_cell_do] * dsw;
                        g_dsw_up = gas_oil_ratio_[up_cell_do] * dsw_up;
                      }

                    g_dsg = gas_oil_ratio_[up_cell_do] * dsg;
                    g_dsg_up = gas_oil_ratio_[up_cell_do] * dsg_up;

                    if (ic == up_cell[d_o] && FI_CHK_RO (main_vars_, ic))
                      {
                        g_dsg += truns_s_mob_og * truns_block_j * cell_pot_k;
                      }
                    else if (FI_CHK_RO (main_vars_, ic))
                      {
                        g_dsg_up += truns_s_mob_og * truns_block_j * cell_pot_k;
                      }
                  }

                if (n_phases == 3)
                  {
                    //oil
                    pp_block_m[p3_oil_po] += dpo;
                    pp_block_m[p3_oil_so] += dso;
                    pp_block_m[p3_oil_sg] += dsg;
                    ps_block_m[p3_oil] += dsw;

                    pp_block_p[p3_oil_po] -= dpo;
                    pp_block_p[p3_oil_so] -= dso;
                    pp_block_p[p3_oil_sg] -= dsg;
                    ps_block_p[p3_oil] += -dsw;

                    pp_block_m_up[p3_oil_po] += dpo_up;
                    pp_block_m_up[p3_oil_so] += dso_up;
                    pp_block_m_up[p3_oil_sg] += dsg_up;
                    ps_block_m_up[p3_oil] += dsw_up;

                    pp_block_p_up[p3_oil_po] -= dpo_up;
                    pp_block_p_up[p3_oil_so] -= dso_up;
                    pp_block_p_up[p3_oil_sg] -= dsg_up;
                    ps_block_p_up[p3_oil] += -dsw_up;

                    // gas
                    pp_block_m[p3_gas_po] += g_dpo;
                    pp_block_m[p3_gas_so] += g_dso;
                    pp_block_m[p3_gas_sg] += g_dsg;
                    ps_block_m[p3_gas] += g_dsw;

                    pp_block_p[p3_gas_po] -= g_dpo;
                    pp_block_p[p3_gas_so] -= g_dso;
                    pp_block_p[p3_gas_sg] -= g_dsg;
                    ps_block_p[p3_gas] += -g_dsw;

                    pp_block_m_up[p3_gas_po] += g_dpo_up;
                    pp_block_m_up[p3_gas_so] += g_dso_up;
                    pp_block_m_up[p3_gas_sg] += g_dsg_up;
                    ps_block_m_up[p3_gas] += g_dsw_up;

                    pp_block_p_up[p3_gas_po] -= g_dpo_up;
                    pp_block_p_up[p3_gas_so] -= g_dso_up;
                    pp_block_p_up[p3_gas_sg] -= g_dsg_up;
                    ps_block_p_up[p3_gas] += -g_dsw_up;
                  }
                else if (is_w)
                  {
                    pp_block_m[p2ow_oil_po] += dpo;
                    pp_block_m[p2ow_oil_so] += dso;
                    ps_block_m[p2ow_oil] += dsw;

                    pp_block_p[p2ow_oil_po] -= dpo;
                    pp_block_p[p2ow_oil_so] -= dso;
                    ps_block_p[p2ow_oil] += -dsw;

                    pp_block_m_up[p2ow_oil_po] += dpo_up;
                    pp_block_m_up[p2ow_oil_so] += dso_up;
                    ps_block_m_up[p2ow_oil] += dsw_up;

                    pp_block_p_up[p2ow_oil_po] -= dpo_up;
                    pp_block_p_up[p2ow_oil_so] -= dso_up;
                    ps_block_p_up[p2ow_oil] += -dsw_up;
                  }
                else if (is_g)
                  {
                    // oil
                    pp_block_m[p2og_oil_po] += dpo;
                    pp_block_m[p2og_oil_sg] += dsg;
                    ps_block_m[p2og_oil] += dso;

                    pp_block_p[p2og_oil_po] -= dpo;
                    pp_block_p[p2og_oil_sg] -= dsg;
                    ps_block_p[p2og_oil] += -dso;

                    pp_block_m_up[p2og_oil_po] += dpo_up;
                    pp_block_m_up[p2og_oil_sg] += dsg_up;
                    ps_block_m_up[p2og_oil] += dso_up;

                    pp_block_p_up[p2og_oil_po] -= dpo_up;
                    pp_block_p_up[p2og_oil_sg] -= dsg_up;
                    ps_block_p_up[p2og_oil] += -dso_up;

                    // gas
                    pp_block_m[p2og_gas_po] += g_dpo;
                    pp_block_m[p2og_gas_sg] += g_dsg;
                    ps_block_m[p2og_gas] += g_dso;

                    pp_block_p[p2og_gas_po] -= g_dpo;
                    pp_block_p[p2og_gas_sg] -= g_dsg;
                    ps_block_p[p2og_gas] += -g_dso;

                    pp_block_m_up[p2og_gas_po] += g_dpo_up;
                    pp_block_m_up[p2og_gas_sg] += g_dsg_up;
                    ps_block_m_up[p2og_gas] += g_dso_up;

                    pp_block_p_up[p2og_gas_po] -= g_dpo_up;
                    pp_block_p_up[p2og_gas_sg] -= g_dsg_up;
                    ps_block_p_up[p2og_gas] += -g_dso_up;
                  }
                else
                  {
                    pp_block_m[0] += dpo;
                    pp_block_p[0] -= dpo;

                    pp_block_m_up[0] += dpo_up;
                    pp_block_p_up[0] -= dpo_up;
                  }
                if (n_phases > 1)
                  {
                    M_MINUS_VV_PROD (n_phases, ps_block_m_up, sp_block_up, pp_block_m_up);
                    V_MINUS_VS_PROD (n_phases, ps_block_m_up, s_rhs_block_up, rhs_m);

                    M_MINUS_VV_PROD (n_phases, ps_block_p_up, sp_block_up, pp_block_p_up);
                    V_MINUS_VS_PROD (n_phases, ps_block_p_up, s_rhs_block_up, rhs_p);
                  }
              }
            // gas
            if (is_g)
              {
                //index_t mobg = 0;
                k = k_j + d_g;

                ps_block_p_up[0] = ps_block_p_up[1] = ps_block_p_up[2] = 0;
                ps_block_m_up[0] = ps_block_m_up[1] = ps_block_m_up[2] = 0;
                sp_block_up = sp_diag_val + sp_size_g;
                s_rhs_block_up = s_rhs_val + rhs_size_g;

                if (up_cell_dg == cell_ind_block[0])
                  {
                    pp_block_m_up = &m_mem[0];
                    pp_block_p_up = &p_mem[0];
                  }
                else
                  {
                    pp_block_m_up = &m_mem[9];
                    pp_block_p_up = &p_mem[9];
                  }

                // calculate pressure derivaties
                dpo = truns_mob_g * truns_block[j];
                dpo_up = 0;
                if (ic == up_cell_dg)
                  dpo += truns_p_mob_g * truns_block[j] * cell_pot[k];
                else
                  dpo_up = truns_p_mob_g * truns_block[j] * cell_pot[k];

                dsw = 0;
                dso = 0;
                dsg = 0;
                dsw_up = dso_up = dsg_up = 0;
                // calculate mobility indexies
                //mobg  = d_g * n_phases + d_g;

                if (is_o)
                  {
                    dsg += truns_mob_g * truns_block[j] * data_[ic].s_deriv_cap_pressure[d_g];
                    dsg_up = 0;
                    if (ic == up_cell_dg)
                      {
                        dsg += truns_s_mob_gg * truns_block[j] * cell_pot[k];
                      }
                    else
                      {
                        dsg_up = truns_s_mob_gg * truns_block[j] * cell_pot[k];
                      }
                  }
                if (n_phases > 1 && !(FI_CHK_SG (main_vars_, up_cell_dg)))
                  {
                    dpo = dpo_up = 0;
                    dso = dso_up = 0;
                    dsg = dsg_up = 0;
                    dsw = dsw_up = 0;
                  }

                if (n_phases == 3)
                  {
                    pp_block_m[p3_gas_po] += dpo;
                    pp_block_m[p3_gas_so] += dso;
                    pp_block_m[p3_gas_sg] += dsg;
                    ps_block_m[p3_gas] += dsw;

                    pp_block_p[p3_gas_po] -= dpo;
                    pp_block_p[p3_gas_so] -= dso;
                    pp_block_p[p3_gas_sg] -= dsg;
                    ps_block_p[p3_gas] += -dsw;

                    // up block
                    pp_block_m_up[p3_gas_po] += dpo_up;
                    pp_block_m_up[p3_gas_so] += dso_up;
                    pp_block_m_up[p3_gas_sg] += dsg_up;
                    ps_block_m_up[p3_gas] += dsw_up;

                    pp_block_p_up[p3_gas_po] -= dpo_up;
                    pp_block_p_up[p3_gas_so] -= dso_up;
                    pp_block_p_up[p3_gas_sg] -= dsg_up;
                    ps_block_p_up[p3_gas] += -dsw_up;
                  }
                else if (n_phases == 2)
                  {
                    pp_block_m[p2og_gas_po] += dpo;
                    pp_block_m[p2og_gas_sg] += dsg;
                    ps_block_m[p2og_gas] += dso;

                    pp_block_p[p2og_gas_po] -= dpo;
                    pp_block_p[p2og_gas_sg] -= dsg;
                    ps_block_p[p2og_gas] += -dso;

                    // up block
                    pp_block_m_up[p2og_gas_po] += dpo_up;
                    pp_block_m_up[p2og_gas_sg] += dsg_up;
                    ps_block_m_up[p2og_gas] += dso_up;

                    pp_block_p_up[p2og_gas_po] -= dpo_up;
                    pp_block_p_up[p2og_gas_sg] -= dsg_up;
                    ps_block_p_up[p2og_gas] += -dso_up;
                  }
                else
                  {
                    pp_block_m[0] += dpo;
                    pp_block_p[0] -= dpo;

                    // up block
                    pp_block_m_up[0] += dpo_up;
                    pp_block_p_up[0] -= dpo_up;
                  }
                if (n_phases > 1)
                  {
                    M_MINUS_VV_PROD (n_phases, ps_block_m_up, sp_block_up, pp_block_m_up);
                    M_MINUS_VV_PROD (n_phases, ps_block_p_up, sp_block_up, pp_block_p_up);
                    V_MINUS_VS_PROD (n_phases, ps_block_m_up, s_rhs_block_up, rhs_m);
                    V_MINUS_VS_PROD (n_phases, ps_block_p_up, s_rhs_block_up, rhs_p);
                  }
              }
            M_MINUS_VV_PROD (n_phases, ps_block_m, sp_block, pp_block_m);
            M_MINUS_VV_PROD (n_phases, ps_block_p, sp_block, pp_block_p);

            V_MINUS_VS_PROD (n_phases, ps_block_m, s_rhs_block, rhs_m);
            V_MINUS_VS_PROD (n_phases, ps_block_p, s_rhs_block, rhs_p);
          }

      }

    public:

      const fi_operator_impl <strategy_t> &fi_operator_;
      index_t                               n_phases;
      bool                                  is_w;
      bool                                  is_g;
      bool                                  is_o;
      index_t                               d_w;
      index_t                               d_g;
      index_t                               d_o;
      const item_array_t                    &saturation_3p_;
      const item_array_t                    &pressure_;
      const item_array_t                    &gas_oil_ratio_;
      const main_var_array_t                &main_vars_;
      const data_array_t                    &data_;
      item_t                                gravity_;
      index_t                               n_sec_vars;
      item_array_t                          &rhs_;
      item_array_t                          &reg_values_;
    };
  } // tpfa

  template <typename strategy_t>
  bool
  fi_operator_impl <strategy_t>::block_connections_mpfa (const item_t &dt)
  {
    //// TODO:
    //// Miryanov: we keep a new version for performance measurement on geostation on large models
    //if (n_phases == 3)
    //  {
    //    return mpfa_calc_3phase <strategy_t> (this).calc (dt);
    //  }
    //else if (n_phases == 2 && is_w)
    //  {
    //    return mpfa_calc_2phase_ow <strategy_t> (this).calc (dt);
    //  }
    //else if (n_phases == 2 && is_g)
    //  {
    //    return mpfa_calc_2phase_og <strategy_t> (this).calc (dt);
    //  }
    //else if (is_w)
    //  {
    //    return mpfa_calc_1phase_w <strategy_t> (this).calc (dt);
    //  }
    //else if (is_g)
    //  {
    //    return mpfa_calc_1phase_g <strategy_t> (this).calc (dt);
    //  }
    //else if (is_o)
    //  {
    //    return mpfa_calc_1phase_o <strategy_t> (this).calc (dt);
    //  }
    //else
    //  {
    //    throw bs_exception ("block_connections_mpfa", "NOT IMPL YET");
    //  }

    index_t b_sqr = n_phases * n_phases;
    index_t equ_w = 0, equ_g = 0, equ_o = 0;

    boost::array <item_t, FI_PHASE_TOT> ave_density;     // average density of phase over all cells involved in MPFA connection

    //potential
    item_t cell_pot[2 * FI_PHASE_TOT];    // potential for cell
    item_t sum_cell_pot[FI_PHASE_TOT];

    //up_stream
    typedef tpfa::mpfa_impl <strategy_t> mpfa_impl_t;
    typename mpfa_impl_t::up_cell_array_t up_cell;    // index of upstream cell

    //rhs
    item_t * rhs_m;
    item_t * rhs_p;

    //jacobian
    boost::array <item_t, 2 * 9> m_mem; 
    boost::array <item_t, 2 * 9> p_mem; 

    assign (flux_rhs_, 0);

    // TPFA: (SHOULD BE REMOVED AFTER FULL MPFA IMPL)
    // number of cells used to build MPFA for one connection
    index_t n_cells_in_conn = 2;          

    if (n_phases == 3)
      {
        equ_w = p3_wat;
        equ_g = p3_gas;
        equ_o = p3_oil;
      }
    else if (n_phases == 2 && is_w)
      {
        equ_w = p2ow_wat;
        equ_o = p2ow_oil;
      }
    else if (n_phases == 2 && is_g)
      {
        equ_g = p2og_gas;
        equ_o = p2og_oil;
      }
    else
      {
        equ_g = equ_w = equ_o = 0;
      }

    mpfa_impl_t mpfa (*this, flux_rhs_, reg_values_);

    // main loop through connections
    for (index_t i = 0; i < n_connections_; ++i)
      {
        index_t row_i = trns_rows_ptr_[i];
        const item_t *truns_block = &trns_values_[row_i];       // temporary array for TPFA connection (SHOULD BE REMOVED AFTER FULL MPFA IMPL)

        // temporary array for TPFA cell ind (SHOULD BE REMOVED AFTER FULL MPFA IMPL)
        const index_t *cell_ind_block = &trns_cols_ptr_[row_i];

        index_t cell_m = cell_ind_block[0];                           //pervii v cols index v stroke          // index of primary cell in MPFA connection (minus)
        index_t cell_p = cell_ind_block[1];                           //vtoroi          // index of primary cell in MPFA connection (plus)
        const index_t *m_mem_ptr = &m_array_[row_i];      // ptr to the block in matrix row cell_m for columns with indexies #cell_ind_block
        const index_t *p_mem_ptr = &p_array_[row_i];      // ptr to the block in matrix row cell_p for columns with indexies #cell_ind_block

        //fills ave_density and ++s counter
        mpfa.mpfa_calc_avarage_density(&ave_density[0], cell_ind_block, n_cells_in_conn);

        //fills cell_pot and sum_cell_pot
        mpfa.mpfa_calc_potential(cell_pot, sum_cell_pot, &ave_density[0], truns_block, depths_,
                            cell_ind_block, cell_m, dt, n_cells_in_conn);

        //fills up_cell
        mpfa.mpfa_calc_upstream(up_cell, sum_cell_pot, cell_m, cell_p);

        //fills rhs_m and rhs_p
        mpfa.mpfa_fill_rhs(rhs_m, rhs_p, up_cell, truns_block, cell_pot, mpfa.rhs_, dt, n_cells_in_conn,
                      cell_m, cell_p, equ_w, equ_g, equ_o);


        //fills m_mem and p_mem
        mpfa.mpfa_fill_jacobian(m_mem, p_mem, up_cell, cell_pot, truns_block, rhs_m, rhs_p,
                           cell_ind_block, sp_diag_, s_rhs_, dt, n_cells_in_conn,
                           cell_m, cell_p, equ_w, equ_g, equ_o);

        item_t * block_m1 = &mpfa.reg_values_[m_mem_ptr[0] * b_sqr];
        item_t * block_p1 = &mpfa.reg_values_[p_mem_ptr[0] * b_sqr];
        item_t * block_m2 = &mpfa.reg_values_[m_mem_ptr[1] * b_sqr];
        item_t * block_p2 = &mpfa.reg_values_[p_mem_ptr[1] * b_sqr];

        for (index_t j = 0; j < b_sqr; ++j)
          {
            block_m1[j] += m_mem[j];
            block_p1[j] += p_mem[j];

            block_m2[j] += m_mem[9 + j];
            block_p2[j] += p_mem[9 + j];
          }

        // END OF JACOBIAN filling
      }

    return 0;
  }

} // namespace blue_sky


#endif  // BS_FI_OPERATOR_BLOCK_CONNECTIONS_MPFA_H_

