/**
 * \file fi_operator_block_connections_mpfa_2.h
 * \brief new (full mpfa) mpfa impl
 * \author Sergey Miryanov
 * \date 06.04.2009
 * */
#ifndef BS_FI_OPERATOR_BLOCK_CONNECTION_MPFA_2_H_
#define BS_FI_OPERATOR_BLOCK_CONNECTION_MPFA_2_H_

#include "calc_model.h"
#include "matrix_vector_op.h"
#include "pp_index.h"

#define IS_STABLE 0
#define UNKNOWN -1

namespace blue_sky {

  namespace mpfa
  {
    #define PSI_W psi_[d_w]
    #define PSI_G psi_[d_g]
    #define PSI_O psi_[d_o]
    #define RHO_W rho_[d_w]
    #define RHO_G rho_[d_g]
    #define RHO_O rho_[d_o]
    #define UP_W up_[d_w]
    #define UP_G up_[d_g]
    #define UP_O up_[d_o]
    #define UP_CELL_W up_cell_[d_w]
    #define UP_CELL_G up_cell_[d_g]
    #define UP_CELL_O up_cell_[d_o]
    #define FLOW_W flow_[d_w]
    #define FLOW_G flow_[d_g]
    #define FLOW_O flow_[d_o]

    #define CAP_PRESSURE_K_W data_[k_cell].cap_pressure [d_w]
    #define CAP_PRESSURE_K_G data_[k_cell].cap_pressure [d_g]
    #define PRESSURE_K pressure_[k_cell]
    #define S_DERIV_CAP_PRESSURE_K_W data_[k_cell].s_deriv_cap_pressure [ds_w]
    #define S_DERIV_CAP_PRESSURE_K_G data_[k_cell].s_deriv_cap_pressure [ds_g]

    #define MOBILITY_UP_W           mobility_up_w
    #define MOBILITY_UP_G           mobility_up_g
    #define MOBILITY_UP_O           mobility_up_o
    #define GOR_UP_O                gor_up_o
    #define SW_DERIV_MOBILITY_UP_W  sw_deriv_mobility_up_w
    #define SW_DERIV_MOBILITY_UP_G  sw_deriv_mobility_up_g
    #define SW_DERIV_MOBILITY_UP_O  sw_deriv_mobility_up_o
    #define SG_DERIV_MOBILITY_UP_W  sg_deriv_mobility_up_w
    #define SG_DERIV_MOBILITY_UP_G  sg_deriv_mobility_up_g
    #define SG_DERIV_MOBILITY_UP_O  sg_deriv_mobility_up_o
    #define SO_DERIV_MOBILITY_UP_W  so_deriv_mobility_up_w
    #define SO_DERIV_MOBILITY_UP_G  so_deriv_mobility_up_g
    #define SO_DERIV_MOBILITY_UP_O  so_deriv_mobility_up_o
    #define P_DERIV_MOBILITY_UP_W   p_deriv_mobility_up_w
    #define P_DERIV_MOBILITY_UP_G   p_deriv_mobility_up_g
    #define P_DERIV_MOBILITY_UP_O   p_deriv_mobility_up_o
    #define MAIN_VAR_UP_G           main_var_up_g
    #define MAIN_VAR_UP_O           main_var_up_o
    #define P_DERIV_GOR_UP_O        p_deriv_gor_up_o

    template <typename strategy_t, bool is_w, bool is_g, bool is_o>
    struct mpfa_base_impl
    {
    public:
      typedef typename strategy_t::item_t             item_t;
      typedef typename strategy_t::rhs_item_t         rhs_item_t;
      typedef typename strategy_t::index_t            index_t;
      typedef typename strategy_t::item_array_t       item_array_t;
      typedef typename strategy_t::rhs_item_array_t   rhs_item_array_t;
      typedef typename strategy_t::index_array_t      index_array_t;
      typedef boost::array <index_t, FI_PHASE_TOT>    up_cell_array_t;
      typedef calc_model <strategy_t>                 calc_model_t;
      typedef typename calc_model_t::data_t           data_t;
      typedef typename calc_model_t::data_array_t     data_array_t;
      typedef typename calc_model_t::main_var_array_t main_var_array_t;

      typedef typename calc_model_t::sp_pvt_dead_oil_array_t  sp_pvt_dead_oil_array_t;
      typedef typename calc_model_t::sp_pvt_gas_array_t       sp_pvt_gas_array_t;
      typedef typename calc_model_t::sp_pvt_water_array_t     sp_pvt_water_array_t;


      typedef rs_mesh_iface <strategy_t>               mesh_iface_t;
      //typedef nc_ptr <mesh_iface_t>                           nc_mesh_iface_t;
      typedef smart_ptr <mesh_iface_t, true>           sp_mesh_iface_t;
      typedef rs_smesh_iface <strategy_t>              rs_smesh_iface_t;


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
          gas_sg = detail::pp_index <n_phases, is_w, is_g>::gas_sg,
          gas_so = detail::pp_index <n_phases, is_w, is_g>::gas_so,
          gas_po = detail::pp_index <n_phases, is_w, is_g>::gas_po,
          oil_sg = detail::pp_index <n_phases, is_w, is_g>::oil_sg,
          oil_so = detail::pp_index <n_phases, is_w, is_g>::oil_so,
          oil_po = detail::pp_index <n_phases, is_w, is_g>::oil_po,
          wat_sg = detail::pp_index <n_phases, is_w, is_g>::wat_sg,
          wat_so = detail::pp_index <n_phases, is_w, is_g>::wat_so,
          wat_po = detail::pp_index <n_phases, is_w, is_g>::wat_po,
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

    public:
      mpfa_base_impl (const fi_operator_impl <strategy_t, is_w, is_g, is_o> &fi_operator_, rhs_item_array_t &rhs_, rhs_item_array_t &reg_values_)
      : d_w (fi_operator_.d_w),
      d_g (fi_operator_.d_g),
      d_o (fi_operator_.d_o),
      ds_w (fi_operator_.ds_w),
      ds_g (fi_operator_.ds_g),
      saturation_3p_ (fi_operator_.saturation_3p_),
      pressure_ (fi_operator_.pressure_),
      gas_oil_ratio_ (fi_operator_.gas_oil_ratio_),
      main_vars_ (fi_operator_.main_vars_),
      data_ (fi_operator_.data_),
      depths_ (fi_operator_.depths_),
      gravity_ (fi_operator_.calc_model_->internal_constants.gravity_constant),
      n_sec_vars (fi_operator_.n_sec_vars),
      sp_diag_ (fi_operator_.sp_diag_),
      s_rhs_ (fi_operator_.s_rhs_),
      reg_values_ (reg_values_),
      rhs_ (rhs_),
      m_mem (0),
      p_mem (0),
      depth_top (0),
      mesh_ (fi_operator_.mesh_)

      {
      }

      BS_FORCE_INLINE void
      fill_jacobian_k_derivs (index_t k) const
      {
        rhs_item_t *jacobian_ik = &reg_values_ [m_mem[k] * b_sqr];
        rhs_item_t *jacobian_jk = &reg_values_ [p_mem[k] * b_sqr];

        for (index_t i = 0; i < b_sqr; ++i)
          {
            jacobian_ik [i] += pp[i];
            jacobian_jk [i] -= pp[i];
          }
      }

      void
      compute_rho (const index_t *cells, index_t cells_count)
      {
        boost::array <item_t, n_phases> counter;
        boost::array <index_t, n_phases> flag;

        assign (rho_, 0);
        assign (counter, 0);
        assign (flag, 0);

        for (index_t j = 0; j < cells_count; ++j)
          {
            index_t ic = cells[j];
            index_t idx = ic * n_phases;
            const data_t &data_ic = data_[ic];

            if (n_phases > 1) // BUG: should be checked
              {
                bool is_g_ = is_g && saturation_3p_[idx + d_g] > EPS_DIV;
                bool is_o_ = is_o && saturation_3p_[idx + d_o] > EPS_DIV;
                bool is_w_ = is_w && saturation_3p_[idx + d_w] > EPS_DIV;

                if (is_g_)
                  {
                    rho_[d_g] += data_ic.density[d_g];
                    ++counter[d_g];
                    ++flag[d_g];
                  }
                if (is_o_)
                  {
                    rho_[d_o] += data_ic.density[d_o];
                    ++counter[d_o];
                    ++flag[d_o];
                  }
                if (is_w_)
                  {
                    rho_[d_w] += data_ic.density[d_w];
                    ++counter[d_w];
                    ++flag[d_w];
                  }
              }
            else
              {
                rho_[0] += data_ic.density[0];
                ++counter[0];
                ++flag[0];
              }
          }
        for (index_t j = 0; j < n_phases; ++j)
          {
            if (flag[j])
              rho_[j] /= counter[j];
          }
      }

      void
      compute_psi (const rhs_item_t *truns, const index_t *cells, index_t cells_count)
      {
        assign (psi_, 0);

        depth_top = depths_ [cells[0]];
        for (index_t k = 0; k < cells_count; ++k)
          {
            index_t k_cell  = cells [k];
            rhs_item_t truns_k  = truns[k];
            item_t g_h      = gravity_ * (depths_[k_cell] - depth_top);
            item_t p        = PRESSURE_K;

            if (is_g) PSI_G += truns_k * (p + CAP_PRESSURE_K_G - RHO_G * g_h);
            if (is_o) PSI_O += truns_k * (p - RHO_O * g_h);
            if (is_w) PSI_W += truns_k * (p + CAP_PRESSURE_K_W - RHO_W * g_h);
        }
      }

      void
      compute_upstream (index_t i, index_t j, index_t i_cell, index_t j_cell)
      {
        index_t up_i[]      = {j, i};
        index_t up_cell[]   = {j_cell, i_cell};

        if (is_g) UP_G      = up_i[PSI_G > 0.0];
        if (is_o) UP_O      = up_i[PSI_O > 0.0];
        if (is_w) UP_W      = up_i[PSI_W > 0.0];
        if (is_g) UP_CELL_G = up_cell[PSI_G > 0.0];
        if (is_o) UP_CELL_O = up_cell[PSI_O > 0.0];
        if (is_w) UP_CELL_W = up_cell[PSI_W > 0.0];

        if (is_w)         sw_deriv_mobility_up_w = data_[UP_CELL_W].s_deriv_mobility[d_w * n_phases + d_w];
        if (is_w && is_g) sw_deriv_mobility_up_g = data_[UP_CELL_G].s_deriv_mobility[d_w * n_phases + d_g];
        if (is_w && is_o) sw_deriv_mobility_up_o = data_[UP_CELL_O].s_deriv_mobility[d_w * n_phases + d_o];
        if (is_g && is_w) sg_deriv_mobility_up_w = data_[UP_CELL_W].s_deriv_mobility[d_g * n_phases + d_w];
        if (is_g)         sg_deriv_mobility_up_g = data_[UP_CELL_G].s_deriv_mobility[d_g * n_phases + d_g];
        if (is_g && is_o) sg_deriv_mobility_up_o = data_[UP_CELL_O].s_deriv_mobility[d_g * n_phases + d_o];
        if (is_o && is_w) so_deriv_mobility_up_w = data_[UP_CELL_W].s_deriv_mobility[d_o * n_phases + d_w];
        if (is_o && is_g) so_deriv_mobility_up_g = data_[UP_CELL_G].s_deriv_mobility[d_o * n_phases + d_g];
        if (is_o)         so_deriv_mobility_up_o = data_[UP_CELL_O].s_deriv_mobility[d_o * n_phases + d_o];
        if (is_w)         p_deriv_mobility_up_w  = data_[UP_CELL_W].p_deriv_mobility[d_w];
        if (is_g)         p_deriv_mobility_up_g  = data_[UP_CELL_G].p_deriv_mobility[d_g];
        if (is_o)         p_deriv_mobility_up_o  = data_[UP_CELL_O].p_deriv_mobility[d_o];
        if (is_w)         mobility_up_w          = data_[UP_CELL_W].mobility [d_w];
        if (is_g)         mobility_up_g          = data_[UP_CELL_G].mobility [d_g];
        if (is_o)         mobility_up_o          = data_[UP_CELL_O].mobility [d_o];
        if (is_o && is_g) gor_up_o               = gas_oil_ratio_[UP_CELL_O];
        if (is_o && is_g) p_deriv_gor_up_o       = data_[UP_CELL_O].p_deriv_gas_oil_ratio;
        if (is_g)         main_var_up_g          = main_vars_[UP_CELL_G];
        if (is_o)         main_var_up_o          = main_vars_[UP_CELL_O];
      }

      void
      compute_flow (const double &dt)
      {
        if (is_g) FLOW_G = dt * MOBILITY_UP_G * PSI_G;
        if (is_o) FLOW_O = dt * MOBILITY_UP_O * PSI_O;
        if (is_w) FLOW_W = dt * MOBILITY_UP_W * PSI_W;

        if (is_g && is_o)
          FLOW_G += FLOW_O * GOR_UP_O;
      }

      void
      fill_rhs (index_t i_cell, index_t j_cell)
      {
        if (is_g)
          {
            rhs_ [i_cell * n_phases + gas_idx] += -FLOW_G;
            rhs_ [j_cell * n_phases + gas_idx] -= -FLOW_G;
          }
        if (is_o)
          {
            rhs_ [i_cell * n_phases + oil_idx] += -FLOW_O;
            rhs_ [j_cell * n_phases + oil_idx] -= -FLOW_O;
          }
        if (is_w)
          {
            rhs_ [i_cell * n_phases + wat_idx] += -FLOW_W;
            rhs_ [j_cell * n_phases + wat_idx] -= -FLOW_W;
          }
      }

      inline item_t
      wat_sw_deriv (index_t k_cell, main_var_type main_var, item_t truns)
      {
        return MOBILITY_UP_W * truns * S_DERIV_CAP_PRESSURE_K_W;
      }
      inline item_t
      wat_sw_deriv_up ()
      {
        return SW_DERIV_MOBILITY_UP_W * PSI_W;
      }
      inline item_t
      wat_po_deriv (index_t k_cell, main_var_type main_var, item_t truns)
      {
        return MOBILITY_UP_W * truns;
      }
      inline item_t
      wat_po_deriv_up ()
      {
        return P_DERIV_MOBILITY_UP_W * PSI_W;
      }
      inline item_t
      wat_sg_deriv_up ()
      {
        return SW_DERIV_MOBILITY_UP_G * PSI_W;
      }
      inline item_t
      wat_so_deriv_up ()
      {
        return SW_DERIV_MOBILITY_UP_O * PSI_W;
      }

      inline item_t
      gas_sw_deriv_up ()
      {
        return SG_DERIV_MOBILITY_UP_W * PSI_G;
      }
      inline item_t
      gas_po_deriv (index_t k_cell, main_var_type main_var, item_t truns)
      {
        return MOBILITY_UP_G * truns;
      }
      inline item_t
      gas_po_deriv_up ()
      {
        return P_DERIV_MOBILITY_UP_G * PSI_G;
      }
      inline item_t
      gas_sg_deriv (index_t k_cell, main_var_type main_var, item_t truns)
      {
        return main_var == FI_SG_VAR ? MOBILITY_UP_G * truns * S_DERIV_CAP_PRESSURE_K_G : 0;
      }
      inline item_t
      gas_sg_deriv_up ()
      {
        return MAIN_VAR_UP_G == FI_SG_VAR ? SG_DERIV_MOBILITY_UP_G * PSI_G : 0;
      }
      inline item_t
      gas_so_deriv_up ()
      {
        return SG_DERIV_MOBILITY_UP_O * PSI_G;
      }

      inline item_t
      oil_sw_deriv_up ()
      {
        return SO_DERIV_MOBILITY_UP_W * PSI_O;
      }
      inline item_t
      oil_po_deriv (index_t k_cell, main_var_type main_var, item_t truns)
      {
        return MOBILITY_UP_O * truns;
      }
      inline item_t
      oil_po_deriv_up ()
      {
        return P_DERIV_MOBILITY_UP_O * PSI_O;
      }
      inline item_t
      oil_sg_deriv_up ()
      {
        return MAIN_VAR_UP_O == FI_SG_VAR ? (SO_DERIV_MOBILITY_UP_G * PSI_O) : (data_[UP_CELL_O].gor_deriv_invers_visc_fvf * data_[UP_CELL_O].relative_perm[d_o] * PSI_O);
      }
      inline item_t
      oil_so_deriv_up ()
      {
        return SO_DERIV_MOBILITY_UP_O * PSI_O;
      }

      /**
       * \brief compute derivatives for k cell
       * \detail
       *    ps [p3_gas]     = 0;
       *    ps [p3_oil]     = 0;
       *    ps [p3_wat]     = wat_sw_deriv (k_cell, main_var, truns);
       *
       *    pp [p3_gas_sg]  = gas_sg_deriv (k_cell, main_var, truns);
       *    pp [p3_gas_so]  = 0;
       *    pp [p3_gas_po]  = gas_po_deriv (k_cell, main_var, truns);
       *
       *    pp [p3_oil_sg]  = 0;
       *    pp [p3_oil_so]  = 0;
       *    pp [p3_oil_po]  = oil_po_deriv (k_cell, main_var, truns);
       *
       *    pp [p3_wat_sg]  = 0;
       *    pp [p3_wat_so]  = 0;
       *    pp [p3_wat_po]  = wat_po_deriv (k_cell, main_var, truns);
       *
       *    pp [p3_gas_sg]  += GOR_UP_O * pp [p3_oil_sg];
       *    pp [p3_gas_so]  += GOR_UP_O * pp [p3_oil_so];
       *    pp [p3_gas_po]  += GOR_UP_O * pp [p3_oil_po];
       *    ps [p3_gas]     += GOR_UP_O * ps [p3_oil];
       * */
      void
      compute_k_derivs (index_t k, index_t k_cell, main_var_type main_var, item_t truns, index_t i_cell, index_t j_cell)
      {
         assign (pp, 0);

        // if one phase and phase is gas
        if (is_g && !is_1p) pp [gas_sg] = gas_sg_deriv (k_cell, main_var, truns);
        if (is_g)           pp [gas_po] = gas_po_deriv (k_cell, main_var, truns);
        if (is_o)           pp [oil_po] = oil_po_deriv (k_cell, main_var, truns);
        if (is_w)           pp [wat_po] = wat_po_deriv (k_cell, main_var, truns);
        if (is_w && !is_1p) pp [wat_sw] = wat_sw_deriv (k_cell, main_var, truns);

        if (is_g && is_o)   pp [gas_po] += GOR_UP_O * pp [oil_po];

        if (!is_1p)
          {
            const rhs_item_t *sp  = &sp_diag_ [k_cell * n_phases];
            const item_t sr   = s_rhs_ [k_cell];
            rhs_item_t *ri    = &rhs_ [i_cell * n_phases];
            rhs_item_t *rj    = &rhs_ [j_cell * n_phases];

            m_minus_vv_prod <n_phases>::eliminate (&pp[b_sqr], sp, pp);
            v_minus_vs_prod <n_phases>::eliminate (&pp[b_sqr], sr, ri);
            v_minus_vs_prod <n_phases>::eliminate (&pp[b_sqr], -sr, rj);
          }
        fill_jacobian_k_derivs (k);
      }

      void
      compute_up_derivs (const double &dt, index_t i_cell, index_t j_cell)
      {
        if (is_g && !is_1p) pp [gas_sg]  = dt * gas_sg_deriv_up ();
        if (is_g && !is_1p) pp [gas_so]  = dt * gas_so_deriv_up ();
        if (is_g)           pp [gas_po]  = dt * gas_po_deriv_up ();

        if (is_o && is_g)   pp [oil_sg]  = dt * oil_sg_deriv_up ();
        if (is_o && !is_1p) pp [oil_so]  = dt * oil_so_deriv_up ();
        if (is_o)           pp [oil_po]  = dt * oil_po_deriv_up ();

        if (is_w && is_g)   pp [wat_sg]  = dt * wat_sg_deriv_up ();
        if (is_w && is_o)   pp [wat_so]  = dt * wat_so_deriv_up ();
        if (is_w)           pp [wat_po]  = dt * wat_po_deriv_up ();

        if (is_g && is_w)   pp [gas_sw]  = dt * gas_sw_deriv_up ();
        if (is_o && is_w)   pp [oil_sw]  = dt * oil_sw_deriv_up ();
        if (is_w && !is_1p) pp [wat_sw]  = dt * wat_sw_deriv_up ();

        if (!is_1p)
          eliminate (i_cell, j_cell);

        fill_jacobian ();

        if (is_g && is_o)
          {
            pp [gas_sg]  = GOR_UP_O * dt * oil_sg_deriv_up ();
            pp [gas_so]  = GOR_UP_O * dt * oil_so_deriv_up ();
            pp [gas_po]  = GOR_UP_O * dt * oil_po_deriv_up ();
            pp [gas_sw]  = GOR_UP_O * pp [oil_sw];

            if (MAIN_VAR_UP_O == FI_SG_VAR)
              {
                pp [gas_po] += P_DERIV_GOR_UP_O * FLOW_O;
              }
            else
              {
                pp [gas_sg] += FLOW_O;
              }

            BS_ASSERT (gas_sg == 0) ((int)gas_sg);
            BS_ASSERT (p3_gas == 0) ((int)p3_gas);

            if (!is_1p)
              {
                const rhs_item_t *sp_o = &sp_diag_ [UP_CELL_O * n_phases];
                pp [gas_idx * n_phases + 0] -= pp [gas_sw] * sp_o [0];
                pp [gas_idx * n_phases + 1] -= pp [gas_sw] * sp_o [1];
                pp [gas_idx * n_phases + 2] -= pp [gas_sw] * sp_o [2];

                rhs_ [i_cell * n_phases + gas_idx] -= pp [gas_sw] * s_rhs_ [UP_CELL_O];
                rhs_ [j_cell * n_phases + gas_idx] += pp [gas_sw] * s_rhs_ [UP_CELL_O];
              }

            for (index_t i = 0; i < n_phases; ++i)
              {
                reg_values_ [i + m_mem[UP_O] * b_sqr + gas_idx * n_phases] += pp[i + gas_idx * n_phases];
                reg_values_ [i + p_mem[UP_O] * b_sqr + gas_idx * n_phases] -= pp[i + gas_idx * n_phases];
              }
          }
      }

      void eliminate (index_t i_cell, index_t j_cell)
      {
        BS_ASSERT (!is_1p);
        if (is_g)
          {
            const rhs_item_t *sp_g = &sp_diag_ [UP_CELL_G * n_phases];
            if (n_phases > 0) pp [gas_idx * n_phases + 0] -= pp [gas_sw] * sp_g [0];
            if (n_phases > 1) pp [gas_idx * n_phases + 1] -= pp [gas_sw] * sp_g [1];
            if (n_phases > 2) pp [gas_idx * n_phases + 2] -= pp [gas_sw] * sp_g [2];
          }

        if (is_o)
          {
            const rhs_item_t *sp_o = &sp_diag_ [UP_CELL_O * n_phases];
            if (n_phases > 0) pp [oil_idx * n_phases + 0] -= pp [oil_sw] * sp_o [0];
            if (n_phases > 1) pp [oil_idx * n_phases + 1] -= pp [oil_sw] * sp_o [1];
            if (n_phases > 2) pp [oil_idx * n_phases + 2] -= pp [oil_sw] * sp_o [2];
          }

        if (is_w)
          {
            const rhs_item_t *sp_w = &sp_diag_ [UP_CELL_W * n_phases];
            if (n_phases > 0) pp [wat_idx * n_phases + 0] -= pp [wat_sw] * sp_w [0];
            if (n_phases > 1) pp [wat_idx * n_phases + 1] -= pp [wat_sw] * sp_w [1];
            if (n_phases > 2) pp [wat_idx * n_phases + 2] -= pp [wat_sw] * sp_w [2];
          }

        if (is_g) rhs_ [i_cell * n_phases + gas_idx]  -= pp [gas_sw] * s_rhs_ [UP_CELL_G];
        if (is_o) rhs_ [i_cell * n_phases + oil_idx]  -= pp [oil_sw] * s_rhs_ [UP_CELL_O];
        if (is_w) rhs_ [i_cell * n_phases + wat_idx]  -= pp [wat_sw] * s_rhs_ [UP_CELL_W];

        if (is_g) rhs_ [j_cell * n_phases + gas_idx]  += pp [gas_sw] * s_rhs_ [UP_CELL_G];
        if (is_o) rhs_ [j_cell * n_phases + oil_idx]  += pp [oil_sw] * s_rhs_ [UP_CELL_O];
        if (is_w) rhs_ [j_cell * n_phases + wat_idx]  += pp [wat_sw] * s_rhs_ [UP_CELL_W];
      }

      void
      fill_jacobian ()
      {
        for (index_t i = 0; i < n_phases; ++i)
          {
            if (is_g) reg_values_ [i + m_mem[UP_G] * b_sqr + gas_idx * n_phases]  += pp[i + gas_idx * n_phases];
            if (is_o) reg_values_ [i + m_mem[UP_O] * b_sqr + oil_idx * n_phases]  += pp[i + oil_idx * n_phases];
            if (is_w) reg_values_ [i + m_mem[UP_W] * b_sqr + wat_idx * n_phases]  += pp[i + wat_idx * n_phases];

            if (is_g) reg_values_ [i + p_mem[UP_G] * b_sqr + gas_idx * n_phases]  -= pp[i + gas_idx * n_phases];
            if (is_o) reg_values_ [i + p_mem[UP_O] * b_sqr + oil_idx * n_phases]  -= pp[i + oil_idx * n_phases];
            if (is_w) reg_values_ [i + p_mem[UP_W] * b_sqr + wat_idx * n_phases]  -= pp[i + wat_idx * n_phases];
          }
      }

      // Salimgareeva EM, 1.07.09


      struct cfl_info
      {
        rhs_item_array_t f11_i;
        rhs_item_array_t f12_i;
        rhs_item_array_t f21_i;
        rhs_item_array_t f22_i;

        cfl_info (index_t n_elems)
          {
            f11_i.assign (n_elems, 0.0f);
            f12_i.assign (n_elems, 0.0f);
            f21_i.assign (n_elems, 0.0f);
            f22_i.assign (n_elems, 0.0f);
          }

        void clear (index_t n_elems)
          {
            f11_i.assign (n_elems, 0.0f);
            f12_i.assign (n_elems, 0.0f);
            f21_i.assign (n_elems, 0.0f);
            f22_i.assign (n_elems, 0.0f);
          }
        rhs_item_t calc_F_i (item_t dt, item_t i)
          {
            return dt*fabs(0.5*(f11_i[i]+f22_i[i]+sqrt((f11_i[i]+f22_i[i])*(f11_i[i]+f22_i[i])-4*(f11_i[i]*f22_i[i]-f12_i[i]*f21_i[i]))));
          }
      };


      void
      fill_cfl (cfl_info &f, const rhs_item_t *truns, index_t i_cell, index_t j_cell, rhs_item_array_t &cfl, item_array_t &saturation_3p)
      {

        int k;
        index_t i, j, k_cell;
        const data_t &data_ic = data_[i_cell];
        const data_t &data_jc = data_[j_cell];

        item_t mobility_tot = 0;
        if (is_g) mobility_tot += MOBILITY_UP_G;
        if (is_o) mobility_tot += MOBILITY_UP_O;
        if (is_w) mobility_tot += MOBILITY_UP_W;

        item_t p_cwoi_deriv = data_[i_cell].s_deriv_cap_pressure[ds_w];
        item_t p_cwoj_deriv = data_[j_cell].s_deriv_cap_pressure[ds_w];
        item_t p_cgoi_deriv = 0;
        item_t p_cgoj_deriv = 0;
        item_t p_cwoi       = data_[i_cell].cap_pressure[ds_w];
        item_t p_cwoj       = data_[j_cell].cap_pressure[ds_w];
        item_t p_cgoi       = 0;
        item_t p_cgoj       = 0;


        if (is_g == 1)
          {
            p_cgoi_deriv = data_[i_cell].s_deriv_cap_pressure[ds_g];
            p_cgoj_deriv = data_[j_cell].s_deriv_cap_pressure[ds_g];
            p_cgoi = data_[i_cell].cap_pressure[ds_g];
            p_cgoj = data_[j_cell].cap_pressure[ds_g];
          }


        item_t delta_psi_w = fabs(psi_[d_w]);
        item_t delta_psi_o = fabs(psi_[d_o]);

        item_t delta_psi_g;
        if (is_g == 1) delta_psi_g = fabs(psi_[d_g]);


        if (is_g == 0)
          {
            f.f11_i[i_cell] += /*truns[i_cell]*/(MOBILITY_UP_O*SW_DERIV_MOBILITY_UP_W*delta_psi_w
            - MOBILITY_UP_W*SW_DERIV_MOBILITY_UP_O*delta_psi_o - truns[0]*MOBILITY_UP_W*MOBILITY_UP_O
            *(p_cwoi_deriv + p_cwoj_deriv))/mobility_tot;

            f.f11_i[j_cell] = f.f11_i[i_cell];
          }
        else
          {
            f.f11_i[i_cell] += /*truns[i_cell]*/((MOBILITY_UP_G + MOBILITY_UP_O)*SW_DERIV_MOBILITY_UP_W*delta_psi_w
            - MOBILITY_UP_W*SW_DERIV_MOBILITY_UP_O*delta_psi_o - truns[0]*MOBILITY_UP_W*(MOBILITY_UP_G + MOBILITY_UP_O)
            *(p_cwoi_deriv + p_cwoj_deriv))/mobility_tot;

            f.f11_i[j_cell] = f.f11_i[i_cell];
          }

        if (is_g == 0)
          {
            f.f12_i[i_cell] = 0;
            f.f12_i[j_cell] = f.f12_i[i_cell];
          }

        else
          {
            f.f12_i[i_cell] += -/*truns[i_cell]*/((MOBILITY_UP_W*SG_DERIV_MOBILITY_UP_O*delta_psi_o) + MOBILITY_UP_W*SG_DERIV_MOBILITY_UP_G
            *delta_psi_g-(MOBILITY_UP_O + MOBILITY_UP_G)*SG_DERIV_MOBILITY_UP_W*delta_psi_w + truns[0]*MOBILITY_UP_W*MOBILITY_UP_G
            *(p_cgoi_deriv + p_cgoj_deriv))/mobility_tot;

            f.f12_i[j_cell] = f.f12_i[i_cell];
          }

        if (is_g == 0)
          {
            f.f21_i[i_cell] = 0;
            f.f21_i[j_cell] = f.f21_i[i_cell];
          }

        else
          {
            f.f21_i[i_cell] += -/*truns[i_cell]*/((MOBILITY_UP_G*SW_DERIV_MOBILITY_UP_W)*delta_psi_w + MOBILITY_UP_G*SW_DERIV_MOBILITY_UP_O
            *delta_psi_o - truns[0]*MOBILITY_UP_G*MOBILITY_UP_W*(p_cgoi_deriv + p_cgoj_deriv))/mobility_tot;

            f.f21_i[j_cell] = f.f21_i[i_cell];//((MOBILITY_UP_G*SW_DERIV_MOBILITY_UP_W)*delta_psi_w + MOBILITY_UP_G*SW_DERIV_MOBILITY_UP_O
            //*delta_psi_o - truns[1]*MOBILITY_UP_G*MOBILITY_UP_W*(p_cgoi_deriv + p_cgoj_deriv))/mobility_tot;
          }

          if (is_g == 0)
          {
            f.f22_i[i_cell] = 0;
            f.f22_i[j_cell] = f.f22_i[i_cell];
          }

        else
          {
            f.f22_i[i_cell] += /*truns[i_cell]*/(-MOBILITY_UP_G*SG_DERIV_MOBILITY_UP_O*delta_psi_o + (MOBILITY_UP_W + MOBILITY_UP_O)
            *SG_DERIV_MOBILITY_UP_G*delta_psi_g - MOBILITY_UP_G*SG_DERIV_MOBILITY_UP_W*delta_psi_w
            + truns[0]*MOBILITY_UP_G*(MOBILITY_UP_O + MOBILITY_UP_W)*(p_cgoi_deriv + p_cgoj_deriv))/mobility_tot;

            f.f22_i[j_cell] = f.f22_i[i_cell];
          }

          // check addition condition
          if (fabs (saturation_3p[i_cell * n_phases + d_w] - saturation_3p[j_cell * n_phases + d_w]) > 0.0001 ||
            fabs (saturation_3p[i_cell * n_phases + d_o] - saturation_3p[j_cell * n_phases + d_o]) > 0.0001)
            {
              cfl[i_cell] = UNKNOWN;
              cfl[j_cell] = UNKNOWN;
            }
      }

    public:

      boost::array <item_t, b_sqr + n_phases> pp;

      index_t                                 d_w;
      index_t                                 d_g;
      index_t                                 d_o;
      index_t                                 ds_w;
      index_t                                 ds_g;
      const item_array_t                      &saturation_3p_;
      const item_array_t                      &pressure_;
      const item_array_t                      &gas_oil_ratio_;
      const main_var_array_t                  &main_vars_;
      const data_array_t                      &data_;
      const item_array_t                      &depths_;
      item_t                                  gravity_;
      index_t                                 n_sec_vars;

      const rhs_item_array_t                  &sp_diag_;
      const rhs_item_array_t                  &s_rhs_;
      rhs_item_array_t                        &reg_values_;
      rhs_item_array_t                        &rhs_;

      boost::array <item_t, n_phases>         rho_;
      boost::array <item_t, n_phases>         psi_;
      boost::array <index_t, n_phases>        up_;
      boost::array <index_t, n_phases>        up_cell_;
      boost::array <item_t, n_phases>         flow_;

      const index_t                           *m_mem;
      const index_t                           *p_mem;

      item_t                                  depth_top;
      item_t                                  sw_deriv_mobility_up_w;
      item_t                                  sw_deriv_mobility_up_g;
      item_t                                  sw_deriv_mobility_up_o;
      item_t                                  sg_deriv_mobility_up_w;
      item_t                                  sg_deriv_mobility_up_g;
      item_t                                  sg_deriv_mobility_up_o;
      item_t                                  so_deriv_mobility_up_w;
      item_t                                  so_deriv_mobility_up_g;
      item_t                                  so_deriv_mobility_up_o;
      item_t                                  p_deriv_mobility_up_w;
      item_t                                  p_deriv_mobility_up_g;
      item_t                                  p_deriv_mobility_up_o;
      item_t                                  mobility_up_w;
      item_t                                  mobility_up_g;
      item_t                                  mobility_up_o;
      item_t                                  gor_up_o;
      item_t                                  p_deriv_gor_up_o;
      main_var_type                           main_var_up_g;
      main_var_type                           main_var_up_o;
      const sp_mesh_iface_t                  &mesh_;
    };

    template <typename fi_operator_t, typename mpfa_t>
    void
    block_connections_mpfa (typename mpfa_t::item_t dt, fi_operator_t &fi_operator, mpfa_t mpfa_)
    {
      typedef typename mpfa_t::index_t      index_t;
      typedef typename mpfa_t::item_t       item_t;
      typedef typename mpfa_t::rhs_item_t   rhs_item_t;
      typedef typename mpfa_t::item_array_t item_array_t;

      typedef rs_mesh_iface <typename fi_operator_t::strategy_t>    mesh_iface_t;
      typedef smart_ptr <mesh_iface_t, true>                        sp_mesh_iface_t;
      typedef rs_smesh_iface <typename fi_operator_t::strategy_t>   rs_smesh_iface_t;

#ifdef BS_BOS_CORE_USE_CSR_ILU_CFL_PREC
      index_t n_active_cells = fi_operator.n_cells_;
      bool use_cfl = fi_operator.calc_model_->ts_params->get_bool (fi_params::USE_CFL);
      static typename mpfa_t::cfl_info cfl (n_active_cells);

      if (use_cfl)
        {
          cfl.clear (n_active_cells);
          assign (fi_operator.cfl_, n_active_cells, 0);
        }
#endif

      // trns: rows - number of cell = row, cols are indies neigh. cells
      // cycle by connections
      for (index_t con = 0; con < fi_operator.n_connections_; ++con)
        {
          const index_t row[] = {fi_operator.trns_rows_ptr_[con], fi_operator.trns_rows_ptr_[con + 1]};
          BS_ASSERT (row[1] - row[0] >= 2) (row[0]) (row[1]);

          const index_t *cells = &fi_operator.trns_cols_ptr_[row[0]];
          const rhs_item_t *truns = &fi_operator.trns_values_[row[0]];

          index_t i   = cells[0];
          index_t j   = cells[1];

          mpfa_.m_mem = &fi_operator.m_array_[row[0]];
          mpfa_.p_mem = &fi_operator.p_array_[row[0]];

          mpfa_.compute_rho (cells, row[1] - row[0]);
          mpfa_.compute_psi (truns, cells, row[1] - row[0]);
          mpfa_.compute_upstream (0, 1, i, j);
          mpfa_.compute_flow (dt);

          mpfa_.fill_rhs (i, j);

          for (index_t k = 0, kcnt = row[1] - row[0]; k < kcnt; ++k)
            {
              mpfa_.compute_k_derivs (k, cells[k], fi_operator.main_vars_[cells[k]], dt * truns [k], i, j);
            }

          mpfa_.compute_up_derivs (dt, i, j);


#ifdef BS_BOS_CORE_USE_CSR_ILU_CFL_PREC
          if (use_cfl)
            {
              // count cfl vector, needs for ACPR
              mpfa_.fill_cfl (cfl, truns, i, j, fi_operator.cfl_, fi_operator.saturation_3p_);
            }
#endif
        }

#ifdef BS_BOS_CORE_USE_CSR_ILU_CFL_PREC
      if (use_cfl)
        {
          for (index_t k = 0; k < n_active_cells; ++k)
            {
              if (fi_operator.cfl_[k] == UNKNOWN)
                {
                  fi_operator.cfl_[k] = cfl.calc_F_i (dt, k) / fi_operator.volume_[k];
                }
              else
                fi_operator.cfl_[k] = IS_STABLE;
            }
        }
#endif
    }
  } // namespace mpfa

  template <typename strategy_t, bool is_w, bool is_g, bool is_o>
  bool
  fi_operator_impl <strategy_t, is_w, is_g, is_o>::block_connections_mpfa (const item_t &dt)
  {
    mpfa::block_connections_mpfa (dt, *this, mpfa::mpfa_base_impl <strategy_t, is_w, is_g, is_o> (*this, flux_rhs_, reg_values_));

    //if (n_phases == 3)
    //  {
    //    mpfa::block_connections_mpfa (dt, *this, mpfa::mpfa_base_impl <strategy_t, true, true, true> (*this, flux_rhs_, reg_values_));
    //  }
    //else if (is_w && is_o)
    //  {
    //    mpfa::block_connections_mpfa (dt, *this, mpfa::mpfa_base_impl <strategy_t, true, false, true> (*this, rhs_, reg_values_));
    //  }
    //else if (is_g && is_o)
    //  {
    //    mpfa::block_connections_mpfa (dt, *this, mpfa::mpfa_base_impl <strategy_t, false, true, true> (*this, rhs_, reg_values_));
    //  }
    //else if (is_w)
    //  {
    //    mpfa::block_connections_mpfa (dt, *this, mpfa::mpfa_base_impl <strategy_t, true, false, false> (*this, rhs_, reg_values_));
    //  }
    //else if (is_g)
    //  {
    //    mpfa::block_connections_mpfa (dt, *this, mpfa::mpfa_base_impl <strategy_t, false, true, false> (*this, rhs_, reg_values_));
    //  }
    //else if (is_o)
    //  {
    //    mpfa::block_connections_mpfa (dt, *this, mpfa::mpfa_base_impl <strategy_t, false, false, true> (*this, rhs_, reg_values_));
    //  }
    //else
    //  {
    //    throw bs_exception ("block_connections_mpfa", "unknown phase value");
    //  }

    return false;
  }



} // namespace blue_sky



#endif // BS_FI_OPERATOR_BLOCK_CONNECTION_MPFA_2_H_

