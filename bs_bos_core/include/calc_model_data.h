/**
 *       \file  calc_model_data.h
 *      \brief  calc_model data holder
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  17.07.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_BOS_CORE_CALC_MODEL_DATA_H_
#define BS_BOS_CORE_CALC_MODEL_DATA_H_

#include "constants.h"

namespace blue_sky {

  /**
   * \class calc_model_data
   * \brief calc_model data holder, holds calculated
   *        values for each mesh cell
   * */
  struct calc_model_data
    {
      typedef calc_model_data                   this_t;
      typedef smart_ptr<this_t, true>           sp_this_t;
      typedef t_float                           item_t;
      typedef t_long                            index_t;

      typedef boost::array <item_t, FI_PHASE_TOT>                 item_array_N_t;   //!< type for store N-Phase values
      typedef boost::array <item_t, FI_PHASE_TOT - 1>             item_array_N_1_t; //!< type for store N-1-Phase values, for example for 3phase model stores 2phase values
      typedef boost::array <item_t, FI_PHASE_TOT * FI_PHASE_TOT>  item_array_N_N_t; //!< type for store N*N-Phase values

      item_array_N_1_t      cap_pressure;               //!< capillary pressure
      item_array_N_1_t      s_deriv_cap_pressure;       //!< deriv of capillary pressure by saturation

      item_array_N_t        relative_perm;              //!< relative permability
      item_array_N_N_t      s_deriv_relative_perm;      //!< deriv of relative permability by saturation

      item_t                p_deriv_gas_oil_ratio;      //!< deriv of gas_oil_ratio (see calc_model) by pressure

      item_array_N_t        invers_fvf;                 //!< inversed value of formation volume factor (iFVF)
      item_array_N_t        p_deriv_invers_fvf;         //!< deriv of iFVF by pressure
      item_t                gor_deriv_invers_fvf;       //!< deriv of iFVF by gas_oil_ratio

      item_array_N_t        invers_viscosity;           //!< inversed value of viscosity (iVISC)
      item_array_N_t        p_deriv_invers_viscosity;   //!< deriv of iVISC by pressure
      item_t                gor_deriv_invers_viscosity; //!< deriv of iVISC by gas_oil_ratio

      item_array_N_t        invers_visc_fvf;            //!< multiplication of iFVF and iVISC (iVISC_FVF)
      item_array_N_t        p_deriv_invers_visc_fvf;    //!< deriv of iVISC_FVF by pressure
      item_t                gor_deriv_invers_visc_fvf;  //!< deriv of iVISC_FVF by gas_oil_ratio

      item_array_N_t        density;                    //!< density
      item_array_N_t        p_deriv_density;            //!< deriv of density by pressure
      item_t                gor_deriv_density;          //!< deriv of density by gas_oil_ratio

      item_t                porosity;                   //!< porosity
      item_t                p_deriv_porosity;           //!< deriv of porosity by pressure

      item_t                truns_mult;                 //!< trunsmissibility multipliers
      item_t                p_deriv_truns_mult;         //!< deriv of truns. multipliers by pressure

      item_array_N_t        mobility;                   //!< mobility
      item_array_N_t        p_deriv_mobility;           //!< deriv of mobility by pressure
      item_array_N_N_t      s_deriv_mobility;           //!< deriv of mobility by saturation

      item_array_N_t        prev_fluid_volume;          //!< fluid volume on previous step


      bool
      operator== (const this_t & /*rhs*/)
      {
        bs_throw_exception ("calc_model_data objects are not comparable");
      }

    };

} // namespace blue_sky

#endif // #ifndef BS_BOS_CORE_CALC_MODEL_DATA_H_

