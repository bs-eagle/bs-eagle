/**
 *       \file  fi_operator_calc_prev_fluid_volume.h
 *      \brief  Calculates fluid volume on previous step
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  23.01.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_FI_OPERATOR_CALC_PREV_FLUID_VOLUME_H_
#define BS_FI_OPERATOR_CALC_PREV_FLUID_VOLUME_H_

namespace blue_sky {

  /**
   * \brief  Calculates fluid volume on previous step
   * */
  template <class strategy_t, bool is_w, bool is_g, bool is_o>
  BS_FORCE_INLINE void
  fi_operator_impl <strategy_t, is_w, is_g, is_o>::calc_prev_fluid_volume ()
  {
    BS_ASSERT (calc_model_->phase_d[FI_PHASE_WATER] == calc_model_->sat_d[FI_PHASE_WATER]) (calc_model_->phase_d[FI_PHASE_WATER]) (calc_model_->sat_d[FI_PHASE_WATER]);
    BS_ASSERT (calc_model_->phase_d[FI_PHASE_GAS]   == calc_model_->sat_d[FI_PHASE_GAS]) (calc_model_->phase_d[FI_PHASE_GAS]) (calc_model_->sat_d[FI_PHASE_GAS]);

    for (index_t i = 0; i < n_cells_; ++i)
      {
        index_t i_w = FI_PH_IND (i, d_w, n_phases);
        index_t i_g = FI_PH_IND (i, d_g, n_phases);
        index_t i_o = FI_PH_IND (i, d_o, n_phases);

        item_t sat_w = 0.0, sat_g = 0.0, sat_o = 0.0;
        // set up saturation
        if (is_w && n_phases > 1)
          sat_w = saturation_3p_[i_w];
        if (is_g && n_phases > 1)
          sat_g = saturation_3p_[i_g];

        if (is_o && n_phases > 1)
          sat_o = saturation_3p_[i_o];
#if 0
        if (is_o)
          {
            if (is_g && is_w)
              sat_o = 1.0 - sat_w - sat_g;
            else if (is_g)
              sat_o = 1.0 - sat_g;
            else if (is_w)
              sat_o = 1.0 - sat_w;
          }
#endif //0

        if (is_o)
          {
            if (n_phases > 1)
              data_[i].prev_fluid_volume[d_o] = data_[i].porosity * data_[i].invers_fvf[d_o] * sat_o;
            else 
              data_[i].prev_fluid_volume[d_o] = data_[i].porosity * data_[i].invers_fvf[d_o];
          }

        if (is_w)
          {
            if (n_phases > 1)
              data_[i].prev_fluid_volume[d_w] = data_[i].porosity * data_[i].invers_fvf[d_w] * sat_w;
            else 
              data_[i].prev_fluid_volume[d_w] = data_[i].porosity * data_[i].invers_fvf[d_w];
          }

        if (is_g)
          {
            if (is_o)
              {
                data_[i].prev_fluid_volume[d_g] = data_[i].porosity * gas_oil_ratio_[i] * sat_o * data_[i].invers_fvf[d_o];
                if (FI_CHK_SG (main_vars_, i))
                  data_[i].prev_fluid_volume[d_g] += data_[i].porosity * data_[i].invers_fvf[d_g] * sat_g;
              }
            else
              data_[i].prev_fluid_volume[d_g] = data_[i].porosity * data_[i].invers_fvf[d_g];
          }
      }
  }

} // namespace blue_sky

#endif  // #ifndef BS_FI_OPERATOR_CALC_PREV_FLUID_VOLUME_H_

