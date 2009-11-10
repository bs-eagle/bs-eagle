/**
 * \file fi_operator_calc_porosity_and_deriv.h
 * \brief calculate porosity and derivativies, also calculate trunsmissibility multiplier
 * \author Sergey Miryanov
 * \date 23.01.2009
 * */
#ifndef BS_FI_OPERATOR_CALC_POROSITY_AND_DERIV_H_
#define BS_FI_OPERATOR_CALC_POROSITY_AND_DERIV_H_

namespace blue_sky {

  /**
  * @brief calculate porosity and derivativies, also calculate trunsmissibility multiplier
  *
  * @param i -- cell index
  * @param pvt_reg -- PVT region index
  * @param poro -- <RETURN> porosity
  * @param dp_poro -- <RETURN> porosity derivative
  * @param t_mult  -- <RETURN> trunsmissibility multiplier
  * @param dp_t_mult -- <RETURN> trunsmissibility multiplier derivative
  */
  template <class strategy_t, bool is_w, bool is_g, bool is_o>
  BS_FORCE_INLINE void
  fi_operator_impl <strategy_t, is_w, is_g, is_o>::calc_porosity_and_deriv (index_t i, 
                                                          index_t pvt_reg, 
                                                          item_t *poro, 
                                                          item_t *dp_poro,
                                                          item_t *t_mult, 
                                                          item_t *dp_t_mult)
  {
    // if ROCKTAB keyword specified use it
    if (calc_model_->rocktab.size ()) // use rocktab
      {
        item_t phi_m = 0, d_phi_m = 0;

        // get region
        index_t reg = calc_model_->rock_regions[i];
        // interpolate
        calc_model_->rocktab[reg].interpolate (pressure_[i], &phi_m, &d_phi_m, t_mult, dp_t_mult);
        // calculate porosity
        *poro = poro_array_[i] * phi_m;
        if (*poro > 1)
          {
            *poro = 1;
          }
        // calculate porosity derivative
        *dp_poro = poro_array_[i] * d_phi_m;
      }
    else // use compressibility to calculate porosity
      {
        // calculate porisity and derivate
        item_t d = rock_grid_comp_const_[pvt_reg] * (pressure_[i] - rock_grid_comp_ref_pressure_[pvt_reg]);

        //*poro = rock_grid_prop.porosity_p_ref[i] * (1.0 + d + d * d * 0.5);
        //*dp_poro = rock_grid_prop.porosity_p_ref[i] * rock_grid_prop.comp_const[pvt_reg] * (1.0 + d);
        *poro = poro_array_[i] * (1.0 + d);
        *dp_poro = poro_array_[i] * rock_grid_comp_const_[pvt_reg];
        *t_mult = 1.0;
        *dp_t_mult = 0.0;
      }
  }

} // namespace blue_sky

#endif  // #ifndef BS_FI_OPERATOR_CALC_POROSITY_AND_DERIV_H_

