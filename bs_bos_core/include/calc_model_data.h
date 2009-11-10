/**
 * \file calc_model_data.h
 * \brief calc_model data holder
 * \author Sergey Miryanov
 * \date 17.07.2009
 * */
#ifndef BS_BOS_CORE_CALC_MODEL_DATA_H_
#define BS_BOS_CORE_CALC_MODEL_DATA_H_

namespace blue_sky {

  template <typename strategy_t>
  struct calc_model_data
    {
      typedef calc_model_data<strategy_t>       this_t;
      typedef smart_ptr<this_t, true>           sp_this_t;
      typedef typename strategy_t::item_t       item_t;
      typedef typename strategy_t::index_t      index_t;

//#ifdef _DEBUG
      typedef boost::array <item_t, FI_PHASE_TOT>									item_array_N_t;
      typedef boost::array <item_t, FI_PHASE_TOT - 1>							item_array_N_1_t;
      typedef boost::array <item_t, FI_PHASE_TOT * FI_PHASE_TOT>	item_array_N_N_t;
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

      item_array_N_1_t      cap_pressure;
      item_array_N_1_t      s_deriv_cap_pressure;

      item_array_N_t        relative_perm;
      item_array_N_N_t      s_deriv_relative_perm;

      item_t                p_deriv_gas_oil_ratio;

      item_array_N_t        invers_fvf;
      item_array_N_t        p_deriv_invers_fvf;
      item_t                gor_deriv_invers_fvf;

      item_array_N_t        invers_viscosity;
      item_array_N_t        p_deriv_invers_viscosity;
      item_t                gor_deriv_invers_viscosity;

      item_array_N_t        invers_visc_fvf;
      item_array_N_t        p_deriv_invers_visc_fvf;
      item_t                gor_deriv_invers_visc_fvf;

      item_array_N_t        density;
      item_array_N_t        p_deriv_density;
      item_t                gor_deriv_density;

      item_t                porosity;
      item_t                p_deriv_porosity;

      item_t                truns_mult;
      item_t                p_deriv_truns_mult;

      item_array_N_t        mobility;
      item_array_N_t        p_deriv_mobility;
      item_array_N_N_t      s_deriv_mobility;

      item_array_N_t        prev_fluid_volume;

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

