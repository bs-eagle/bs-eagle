/**
 *       \file  fi_operator_calc_step_dt_mult.h
 *      \brief  Calculates dt at the first newton iteration on time step
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  23.01.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_FI_OPERATOR_CALC_STEP_DT_MULT_H_
#define BS_FI_OPERATOR_CALC_STEP_DT_MULT_H_

namespace blue_sky {

  /**
   * \brief  Calculates dt at the first newton iteration on time step
   * \param  prev_mult Previous dt multiplier
   * \param  max_res
   * \return New dt multiplier
   * */
  template <class strategy_t, bool is_w, bool is_g, bool is_o>
  BS_FORCE_INLINE typename strategy_t::item_t
  fi_operator_impl <strategy_t, is_w, is_g, is_o>::calc_step_dt_mult (item_t prev_mult, item_t max_res)
  {
    item_t coef = 1.0, m_r, d, d_mult;
    index_t i;
    index_t equ_w, equ_g, equ_o;
#ifdef _MPI
    index_t n_left = 0;///mpi_decomp->get_recv_l ();
    index_t n_own = 0;///mpi_decomp->get_n_local_own () + n_left;
    item_t mpi_coef;
#endif //_MPI


#ifdef _MPI
    for (i = n_left; i < n_own; ++i)
#else //_MPI
    for (i = 0; i < n_cells_; ++i)
#endif //_MPI
      {
        if (n_phases == 3)
          {
            equ_w = i * n_phases + p3_wat;
            equ_g = i * n_phases + p3_gas;
            equ_o = i * n_phases + p3_oil;
          }
        else if (n_phases == 2 && is_w)
          {
            equ_w = i * n_phases + p2ow_wat;
            equ_o = i * n_phases + p2ow_oil;
          }
        else if (n_phases == 2 && is_g)
          {
            equ_g = i * n_phases + p2og_gas;
            equ_o = i * n_phases + p2og_oil;
          }
        else
          {
            equ_g = equ_w = equ_o = i * n_phases + 0;
          }
#ifdef _MPI
        equ_w -= n_phases * n_left;
        equ_g -= n_phases * n_left;
        equ_o -= n_phases * n_left;
#endif //_MPI
        d_mult = (calc_model_->ave_volume * data_[i].porosity);
        if (is_w)
          {
            m_r = max_res * d_mult * calc_model_->invers_fvf_average[d_w];
            if (rhs_[equ_w] * flux_rhs_[equ_w] <= EPS_DIV
                && fabs (flux_rhs_[equ_w]) > fabs (rhs_[equ_w])
                && fabs (flux_rhs_[equ_w] + rhs_[equ_w]) > m_r)
              {
                if (flux_rhs_[equ_w] > 0)
                  d = (m_r - rhs_[equ_w]) / flux_rhs_[equ_w];
                else
                  d = -(m_r + rhs_[equ_w]) / flux_rhs_[equ_w];
                if (d < coef)
                  coef = d;
              }
          }
        if (is_g)
          {
            m_r = max_res * d_mult * calc_model_->invers_fvf_average[d_g];
            if (rhs_[equ_g] * flux_rhs_[equ_g] <= EPS_DIV
                && fabs (flux_rhs_[equ_g]) > fabs (rhs_[equ_g])
                && fabs (flux_rhs_[equ_g] + rhs_[equ_g]) > m_r)
              {
                if (flux_rhs_[equ_g] > 0)
                  d = (m_r - rhs_[equ_g]) / flux_rhs_[equ_g];
                else
                  d = -(m_r + rhs_[equ_g]) / flux_rhs_[equ_g];
                if (d < coef)
                  coef = d;
              }
          }
        if (is_o)
          {
            m_r = max_res * d_mult * calc_model_->invers_fvf_average[d_o];
            if (rhs_[equ_o] * flux_rhs_[equ_o] <= EPS_DIV
                && fabs (flux_rhs_[equ_o]) > fabs (rhs_[equ_o])
                && fabs (flux_rhs_[equ_o] + rhs_[equ_o]) > m_r)
              {
                if (flux_rhs_[equ_o] > 0)
                  d = (m_r - rhs_[equ_o]) / flux_rhs_[equ_o];
                else
                  d = -(m_r + rhs_[equ_o]) / flux_rhs_[equ_o];
                if (d < coef)
                  coef = d;
              }
          }

      }
    if (coef < prev_mult)
      coef = prev_mult;

#ifdef _MPI
    MPI_Allreduce (&coef, &mpi_coef, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    coef = mpi_coef;
#endif
    return coef;
  }

}


#endif // #ifndef BS_FI_OPERATOR_CALC_STEP_DT_MULT_H_

