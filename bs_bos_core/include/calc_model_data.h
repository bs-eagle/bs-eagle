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
      typedef t_double                          item_t;
      typedef t_long                            index_t;

//#ifdef _DEBUG
      typedef boost::array <item_t, FI_PHASE_TOT>                 item_array_N_t;   //!< type for store N-Phase values
      typedef boost::array <item_t, FI_PHASE_TOT - 1>             item_array_N_1_t; //!< type for store N-1-Phase values, for example for 3phase model stores 2phase values
      typedef boost::array <item_t, FI_PHASE_TOT * FI_PHASE_TOT>  item_array_N_N_t; //!< type for store N*N-Phase values


      bool
      operator== (const this_t &rhs)
      {
        bs_throw_exception ("calc_model_data objects are not comparable");
      }

//#else
//      typedef array_ext <item_t>      item_array_N_t;
//      typedef array_ext <item_t>      item_array_N_1_t;
//      typedef array_ext <item_t>      item_array_N_N_t;
//#endif

//#ifndef _DEBUG
//      calc_model_data ()
//      : data_ (0)
//      {
//      }
//
//      calc_model_data (size_t n_phases)
//      : cap_pressure (0, n_phases - 1) ,
//      s_deriv_cap_pressure (0, n_phases - 1),
//      relative_perm (0, n_phases),
//      s_deriv_relative_perm (0, n_phases * n_phases),
//      invers_fvf (0, n_phases),
//      p_deriv_invers_fvf (0, n_phases),
//      invers_viscosity (0, n_phases),
//      p_deriv_invers_viscosity (0, n_phases),
//      invers_visc_fvf (0, n_phases),
//      p_deriv_invers_visc_fvf (0, n_phases),
//      density (0, n_phases),
//      p_deriv_density (0, n_phases),
//      mobility (0, n_phases),
//      p_deriv_mobility (0, n_phases),
//      s_deriv_mobility (0, n_phases * n_phases),
//      prev_fluid_volume (0, n_phases),
//      data_ (0)
//      {
//        allocate_data (size ());
//        // we should deallocate memory. but now it is not a problem.
//      }
//#endif

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

//#ifndef _DEBUG
//    private:
//      void
//      allocate_data (size_t data_size)
//      {
//        delete []data_;
//        data_ = aligned_allocator <item_t>::allocate (data_size);
//        memset (data_, 0, sizeof (item_t) * data_size);
//
//        item_t *data = data_;
//        data = cap_pressure.init (data);
//        data = s_deriv_cap_pressure.init (data);
//        data = relative_perm.init (data);
//        data = s_deriv_relative_perm.init (data);
//        data = invers_fvf.init (data);
//        data = p_deriv_invers_fvf.init (data);
//        data = invers_viscosity.init (data);
//        data = p_deriv_invers_viscosity.init (data);
//        data = invers_visc_fvf.init (data);
//        data = p_deriv_invers_visc_fvf.init (data);
//        data = density.init (data);
//        data = p_deriv_density.init (data);
//        data = mobility.init (data);
//        data = p_deriv_mobility.init (data);
//        data = s_deriv_mobility.init (data);
//        data = prev_fluid_volume.init (data);
//      }
//
//      size_t
//      size () const
//      {
//        return cap_pressure.size () 
//          + s_deriv_cap_pressure.size ()
//          + relative_perm.size ()
//          + s_deriv_relative_perm.size ()
//          + invers_fvf.size ()
//          + p_deriv_invers_fvf.size ()
//          + invers_viscosity.size ()
//          + p_deriv_invers_viscosity.size ()
//          + invers_visc_fvf.size ()
//          + p_deriv_invers_visc_fvf.size ()
//          + density.size ()
//          + p_deriv_density.size ()
//          + mobility.size ()
//          + p_deriv_mobility.size ()
//          + s_deriv_mobility.size ()
//          + prev_fluid_volume.size ()
//          ;
//      }
//
//    private:
//      item_t *data_;
//#endif
    };

} // namespace blue_sky

#endif // #ifndef BS_BOS_CORE_CALC_MODEL_DATA_H_

