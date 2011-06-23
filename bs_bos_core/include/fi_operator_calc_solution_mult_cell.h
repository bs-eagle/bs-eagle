/**
 *       \file  fi_operator_calc_solution_mult_cell.h
 *      \brief  Calculates multiplier for solution
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  08.02.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_FI_OPERATOR_CALC_SOLUTION_MULT_CELL_H_
#define BS_FI_OPERATOR_CALC_SOLUTION_MULT_CELL_H_

namespace blue_sky {

  /**
   * \brief  Calculates multipler for solution
   * \param  base_norm Norm storage
   * \return Calculated multiplier
   * */
  template <bool is_w, bool is_g, bool is_o>
  BS_FORCE_INLINE strategy_t::item_t
  fi_operator_impl <is_w, is_g, is_o>::calc_solution_mult_cell (const norms_storage_t &old_norm)
  {
    item_t d = 0.0;
    item_t mult = 0.0;
    const static item_t eps = item_t (1.0e-8);
    norms_storage_t norm;

    for (index_t id = norms::C_ACPV_GAS; id <= norms::C_ACPV_OIL; ++id)
      {
        norm.clear ();
        update_norm_by_cell (old_norm.idx[id], norm);
        norm.val[id] /= calc_model_->ave_volume;

        item_t abs_old = fabs (old_norm.val[id]);
        item_t abs_new = fabs (norm.val[id]);
        if ((old_norm.val[id] * norm.val[id]) < -eps || (abs_old > abs_new))
        {
          mult += abs_new * abs_old / fabs (old_norm.val[id] - norm.val[id]);
          d    += abs_new;
        }
      }

    return (d > eps) ? mult / d : 1.0;
  }
}


#endif  // #ifndef BS_FI_OPERATOR_CALC_SOLUTION_MULT_CELL_H_

