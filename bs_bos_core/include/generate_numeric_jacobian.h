/**
 *       \file  generate_numeric_jacobian.h
 *      \brief  Fills jacobian with derivs that numerically calculated
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  27.11.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Obsolete
 * */
#ifndef BS_MAIN_LOOP_GENERATE_NUMERIC_JACOBIAN_H_
#define BS_MAIN_LOOP_GENERATE_NUMERIC_JACOBIAN_H_

#include "string_formater.h"

namespace blue_sky
  {

  namespace xxx
    {

    template <typename strategy_t>
    struct clmn_holder
      {
        typedef typename strategy_t::index_t  index_t;
        typedef typename strategy_t::item_t   item_t;

        typedef typename strategy_t::item_array_t item_array_t;

        clmn_holder (item_array_t &dst, index_t index, index_t shift, index_t block)
            : dst (dst)
            , index (index)
            , shift (shift)
            , block (block)
        {
        }

        template <typename op_t>
        clmn_holder <strategy_t> &
        operator-= (const op_t &op)
        {
          index_t idx = 0;
          for (index_t i = 0, cnt = (index_t)op.size (), clmn_count = (index_t)op.size () / block; i < clmn_count; ++i)
            {
              index_t i_idx = index * block * block + shift + i * clmn_count * block * block;
              for (index_t j = 0, k = 0; j < block; ++j, ++idx, k += block)
                {
                  dst [i_idx + k] -= op[idx];
                }
            }

          return *this;
        }

        item_array_t  &dst;
        index_t       index;
        index_t       shift;
        index_t       block;
      };

    template <typename array_t>
    struct clmn_sub
      {
        clmn_sub (const array_t &lhs, const array_t &rhs)
            : lhs (lhs)
            , rhs (rhs)
        {
        }

        size_t
        size () const
          {
            return lhs.size ();
          }

        typedef typename array_t::value_type value_type;
        value_type
        operator[] (size_t i) const
          {
            return lhs[i] - rhs[i];
          }

        const array_t &lhs;
        const array_t &rhs;
      };

    template <typename lhs_t>
    struct clmn_sub_div
      {
        typedef typename lhs_t::value_type value_type;

        clmn_sub_div (const lhs_t &lhs, value_type scalar)
            : lhs (lhs)
            , scalar (scalar)
        {
        }

        size_t
        size () const
          {
            return lhs.size ();
          }

        value_type
        operator[] (size_t i) const
          {
            return lhs[i] / scalar;
          }

        const lhs_t &lhs;
        value_type  scalar;
      };
    template <typename lhs_t>
    struct clmn_sub_div_mul
      {
        typedef typename lhs_t::value_type value_type;

        clmn_sub_div_mul (const lhs_t &lhs, value_type scalar)
            : lhs (lhs)
            , scalar (scalar)
        {
        }

        size_t
        size () const
          {
            return lhs.size ();
          }

        value_type
        operator[] (size_t i) const
          {
            return lhs[i] * scalar;
          }

        const lhs_t &lhs;
        value_type scalar;
      };

    template <typename array_t>
    clmn_sub <array_t>
    operator- (const array_t &lhs, const array_t &rhs)
    {
      return clmn_sub <array_t> (lhs, rhs);
    }
    template <typename lhs_t, typename item_t>
    clmn_sub_div <lhs_t>
    operator/ (const lhs_t &lhs, const item_t &scalar)
    {
      return clmn_sub_div <lhs_t> (lhs, scalar);
    }
    template <typename lhs_t, typename item_t>
    clmn_sub_div_mul <lhs_t>
    operator* (const lhs_t &lhs, const item_t &scalar)
    {
      return clmn_sub_div_mul <lhs_t> (lhs, scalar);
    }
  } // namespace xxx


  template <bool is_w, bool is_g, bool is_o>
  inline void
  main_loop_calc <is_w, is_g, is_o>::generate_numeric_jacobian (int /*init*/)
  {
    throw bs_exception ("generate_numeric_jacobian", "NOT IMPL YET");
    //using namespace xxx;

    //locked_jmatrix_t jmatrix_ (jacobian_->get_jmatrix ());

    //calc_model_->fi_operator (dt_, init, init, number_of_small_time_steps, reservoir_, mesh_, jacobian_, false, false);
    //item_array_t init_rhs;
    //init_rhs.assign (jmatrix_->get_rhs ().begin (), jmatrix_->get_rhs ().end ());

    //index_t cells_num       = mesh_->get_n_active_elements();
    //index_t n_phases        = calc_model_->n_phases;
    //index_t b_sqr           = n_phases * n_phases;

    //item_array_t numeric_jacobian;
    //numeric_jacobian.resize (b_sqr * cells_num * cells_num);

    //item_array_t &pressure  = calc_model_->pressure;
    //item_array_t &sat       = calc_model_->saturation_3p;
    //item_array_t &gor       = calc_model_->gas_oil_ratio;

    //const typename calc_model_t::main_var_array_t &main_var = calc_model_->main_variable;

    //bool is_o = calc_model_->is_oil ();
    //bool is_w = calc_model_->is_water ();
    //bool is_g = calc_model_->is_gas ();

    //index_t d_o = n_phases == 3 ? p3_oil : (is_o && is_w ? p2ow_oil : (is_o && is_g ? p2og_oil : 0));
    //index_t d_w = n_phases == 3 ? p3_wat : (is_o && is_w ? p2ow_wat : (0));
    //index_t d_g = n_phases == 3 ? p3_gas : (is_o && is_w ? 0 : (is_o && is_g ? p2og_gas : 0));

    //index_t i_o = calc_model_->phase_d[FI_PHASE_OIL];
    //index_t i_w = calc_model_->phase_d[FI_PHASE_WATER];
    //index_t i_g = calc_model_->phase_d[FI_PHASE_GAS];

#ifdef _DEBUG
    //item_t *p_sat = &sat[0];
#endif

    //item_t eps = 1.0e-8;
    //for (index_t i = 0; i < cells_num; ++i)
    //  {
    //    if (is_g)
    //      {
    //        if (main_var[i] == FI_SG_VAR)
    //          {
    //            sat [i * n_phases + i_g] += eps;
    //            sat [i * n_phases + i_w] -= eps;
    //            calc_model_->fi_operator (dt_, init, init, number_of_small_time_steps, reservoir_, mesh_, jacobian_, false, false);
    //            item_array_t &sat_g_rhs = jmatrix_->get_rhs ();
    //            clmn_holder <strategy_t> (numeric_jacobian, i, d_g, n_phases) -= (sat_g_rhs - init_rhs) / (eps);
    //            sat [i * n_phases + i_g] -= eps;
    //            sat [i * n_phases + i_w] += eps;
    //          }
    //        else if (main_var[i] == FI_RO_VAR)
    //          {
    //            gor [i] += eps;
    //            sat [i * n_phases + i_w] -= eps;
    //            calc_model_->fi_operator (dt_, init, init, number_of_small_time_steps, reservoir_, mesh_, jacobian_, false, false);
    //            item_array_t &sat_g_rhs = jmatrix_->get_rhs ();
    //            clmn_holder <strategy_t> (numeric_jacobian, i, d_g, n_phases) -= (sat_g_rhs - init_rhs) / (eps);
    //            gor [i] -= eps;
    //            sat [i * n_phases + i_w] += eps;
    //          }
    //      }

    //    if (is_o)
    //      {
    //        sat [i * n_phases + i_o] += eps;
    //        sat [i * n_phases + i_w] -= eps;
    //        calc_model_->fi_operator (dt_, init, init, number_of_small_time_steps, reservoir_, mesh_, jacobian_, false, false);
    //        item_array_t &sat_o_rhs = jmatrix_->get_rhs ();
    //        clmn_holder <strategy_t> (numeric_jacobian, i, d_o, n_phases) -= (sat_o_rhs - init_rhs) / (eps);
    //        sat [i * n_phases + i_o] -= eps;
    //        sat [i * n_phases + i_w] += eps;
    //      }

    //    if (is_w)
    //      {
    //        pressure[i] += eps;
    //        calc_model_->fi_operator (dt_, init, init, number_of_small_time_steps, reservoir_, mesh_, jacobian_, false, false);
    //        item_array_t &pressure_rhs = jmatrix_->get_rhs ();
    //        clmn_holder <strategy_t> (numeric_jacobian, i, d_w, n_phases) -= (pressure_rhs - init_rhs) / (eps);
    //        pressure[i] -= eps;
    //      }

    //    tools::save_seq_vector ("numeric_jacobian.txt").save (numeric_jacobian);
    //  }

    //static size_t iter_counter = 0;
    //++iter_counter;
    //tools::save_seq_vector (tools::string_formater ("numeric_jacobian.bs.%d.txt", iter_counter).str).save (numeric_jacobian);
  }

} // namespace blue_sky

#endif  // #ifndef BS_MAIN_LOOP_GENERATE_NUMERIC_JACOBIAN_H_

