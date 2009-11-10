#include "stdafx.h"

#include "calc_model.h"
#include "py_calc_model.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "py_scal_wrapper.h"
#include "py_data_class.h"
#include BS_STOP_PLUGIN_IMPORT ()

// WTF??
#include "well_results_storage.h"
#include "fip_results_storage.h"

#include "export_python_wrapper.h"

using namespace boost::python;

namespace blue_sky
  {
  namespace python
    {

    template < class strategy_t >
    typename py_calc_model_data<strategy_t>::item_array_t
    py_calc_model_data<strategy_t>::get_cap_pressure () const
      {
        item_array_t t;
        t.insert(t.end(), this->cap_pressure.begin(), this->cap_pressure.end());
        return t;
      }

    template < class strategy_t >
    typename py_calc_model_data<strategy_t>::item_array_t
    py_calc_model_data<strategy_t>::get_s_deriv_cap_pressure () const
      {
        item_array_t t;
        t.insert(t.end(), this->s_deriv_cap_pressure.begin(), this->s_deriv_cap_pressure.end());
        return t;
      }

    template < class strategy_t >
    typename py_calc_model_data<strategy_t>::item_array_t
    py_calc_model_data<strategy_t>::get_relative_perm () const
      {
        item_array_t t;
        t.insert(t.end(), this->relative_perm.begin(), this->relative_perm.end());
        return t;
      }

    template < class strategy_t >
    typename py_calc_model_data<strategy_t>::item_array_t
    py_calc_model_data<strategy_t>::get_s_deriv_relative_perm () const
      {
        item_array_t t;
        t.insert(t.end(), this->s_deriv_relative_perm.begin(), this->s_deriv_relative_perm.end());
        return t;
      }


    template < class strategy_t >
    typename py_calc_model<strategy_t>::item_array_t
    py_calc_model<strategy_t>::get_gas_oil_ratio ()
    {
      return this->get_spx<wrapped_t> ()->gas_oil_ratio;
    }

    template < class strategy_t >
    typename py_calc_model_data<strategy_t>::item_t
    py_calc_model_data<strategy_t>::get_p_deriv_gas_oil_ratio () const
      {
        return this->p_deriv_gas_oil_ratio;
      }


    template < class strategy_t >
    typename py_calc_model_data<strategy_t>::item_array_t
    py_calc_model_data<strategy_t>::get_invers_fvf () const
      {
        item_array_t t;
        t.insert(t.end(), this->invers_fvf.begin(), this->invers_fvf.end());
        return t;
      }

    template < class strategy_t >
    typename py_calc_model_data<strategy_t>::item_array_t
    py_calc_model_data<strategy_t>::get_p_deriv_invers_fvf () const
      {
        item_array_t t;
        t.insert(t.end(), this->p_deriv_invers_fvf.begin(), this->p_deriv_invers_fvf.end());
        return t;
      }

    template < class strategy_t >
    typename py_calc_model_data<strategy_t>::item_t
    py_calc_model_data<strategy_t>::get_gor_deriv_invers_fvf () const
      {
        return this->gor_deriv_invers_fvf;
      }


    template < class strategy_t >
    typename py_calc_model_data<strategy_t>::item_array_t
    py_calc_model_data<strategy_t>::get_invers_viscosity () const
      {
        item_array_t t;
        t.insert(t.end(), this->invers_viscosity.begin(), this->invers_viscosity.end());
        return t;
      }

    template < class strategy_t >
    typename py_calc_model_data<strategy_t>::item_array_t
    py_calc_model_data<strategy_t>::get_p_deriv_invers_viscosity () const
      {
        item_array_t t;
        t.insert(t.end(), this->p_deriv_invers_viscosity.begin(), this->p_deriv_invers_viscosity.end());
        return t;
      }

    template < class strategy_t >
    typename py_calc_model_data<strategy_t>::item_t
    py_calc_model_data<strategy_t>::get_gor_deriv_invers_viscosity () const
      {
        return this->gor_deriv_invers_viscosity;
      }


    template < class strategy_t >
    typename py_calc_model_data<strategy_t>::item_array_t
    py_calc_model_data<strategy_t>::get_invers_visc_fvf () const
      {
        item_array_t t;
        t.insert(t.end(), this->invers_visc_fvf.begin(), this->invers_visc_fvf.end());
        return t;
      }

    template < class strategy_t >
    typename py_calc_model_data<strategy_t>::item_array_t
    py_calc_model_data<strategy_t>::get_p_deriv_invers_visc_fvf () const
      {
        item_array_t t;
        t.insert(t.end(), this->p_deriv_invers_visc_fvf.begin(), this->p_deriv_invers_visc_fvf.end());
        return t;
      }

    template < class strategy_t >
    typename py_calc_model_data<strategy_t>::item_t
    py_calc_model_data<strategy_t>::get_gor_deriv_invers_visc_fvf () const
      {
        return this->gor_deriv_invers_visc_fvf;
      }


    template < class strategy_t >
    typename py_calc_model_data<strategy_t>::item_array_t
    py_calc_model_data<strategy_t>::get_density () const
      {
        item_array_t t;
        t.insert(t.end(), this->density.begin(), this->density.end());
        return t;
      }

    template < class strategy_t >
    typename py_calc_model_data<strategy_t>::item_array_t
    py_calc_model_data<strategy_t>::get_p_deriv_density () const
      {
        item_array_t t;
        t.insert(t.end(), this->p_deriv_density.begin(), this->p_deriv_density.end());
        return t;
      }

    template < class strategy_t >
    typename py_calc_model_data<strategy_t>::item_t
    py_calc_model_data<strategy_t>::get_gor_deriv_density () const
      {
        return this->gor_deriv_density;
      }


    template < class strategy_t >
    typename py_calc_model_data<strategy_t>::item_t
    py_calc_model_data<strategy_t>::get_porosity () const
      {
        return this->porosity;
      }

    template < class strategy_t >
    typename py_calc_model_data<strategy_t>::item_t
    py_calc_model_data<strategy_t>::get_p_deriv_porosity () const
      {
        return this->p_deriv_porosity;
      }


    template < class strategy_t >
    typename py_calc_model_data<strategy_t>::item_t
    py_calc_model_data<strategy_t>::get_truns_mult () const
      {
        return this->truns_mult;
      }

    template < class strategy_t >
    typename py_calc_model_data<strategy_t>::item_t
    py_calc_model_data<strategy_t>::get_p_deriv_truns_mult () const
      {
        return this->p_deriv_truns_mult;
      }


    template < class strategy_t >
    typename py_calc_model_data<strategy_t>::item_array_t
    py_calc_model_data<strategy_t>::get_mobility () const
      {
        item_array_t t;
        t.insert(t.end(), this->mobility.begin(), this->mobility.end());
        return t;
      }

    template < class strategy_t >
    typename py_calc_model_data<strategy_t>::item_array_t
    py_calc_model_data<strategy_t>::get_p_deriv_mobility () const
      {
        item_array_t t;
        t.insert(t.end(), this->p_deriv_mobility.begin(), this->p_deriv_mobility.end());
        return t;
      }

    template < class strategy_t >
    typename py_calc_model_data<strategy_t>::item_array_t
    py_calc_model_data<strategy_t>::get_s_deriv_mobility () const
      {
        item_array_t t;
        t.insert(t.end(), this->s_deriv_mobility.begin(), this->s_deriv_mobility.end());
        return t;
      }


    template < class strategy_t >
    typename py_calc_model_data<strategy_t>::item_array_t
    py_calc_model_data<strategy_t>::get_prev_fluid_volume () const
      {
        item_array_t t;
        t.insert(t.end(), this->prev_fluid_volume.begin(), this->prev_fluid_volume.end());
        return t;
      }

    template < class strategy_t >
    py_calc_model<strategy_t>::py_calc_model()
        : py_objbase(wrapped_t::bs_type())
    {}

    template<class strategy_t>
    py_calc_model<strategy_t>::py_calc_model(const sp_cm_t &src)
        : py_objbase(src)
    {
    }

    template < class strategy_t >
    typename py_calc_model<strategy_t>::item_array_t
    py_calc_model<strategy_t>::get_saturation ()
    {
      return this->get_spx<wrapped_t> ()->saturation_3p;
    }

    template < class strategy_t >
    typename py_calc_model<strategy_t>::item_array_t
    py_calc_model<strategy_t>::get_pressure ()
    {
      return this->get_spx<wrapped_t> ()->pressure;
    }

    /*******************************************************/

    template < class strategy_t >
    typename py_calc_model<strategy_t>::item_array_t
    py_calc_model<strategy_t>::get_plane_flow_rate ()
    {
      return this->get_spx<wrapped_t> ()->plane_flow_rate;
    }

    template < class strategy_t >
    typename py_calc_model<strategy_t>::item_array_t
    py_calc_model<strategy_t>::get_full_step_plane_flow_rate ()
    {
      return this->get_spx<wrapped_t> ()->full_step_plane_flow_rate;
    }

    template < class strategy_t >
    typename py_calc_model<strategy_t>::index_array_t
    py_calc_model<strategy_t>::get_pvt_regions ()
    {
      return this->get_spx<wrapped_t> ()->pvt_regions;
    }

    template < class strategy_t >
    typename py_calc_model<strategy_t>::index_array_t
    py_calc_model<strategy_t>::get_sat_regions ()
    {
      return this->get_spx<wrapped_t> ()->sat_regions;
    }

    template < class strategy_t >
    typename py_calc_model<strategy_t>::index_array_t
    py_calc_model<strategy_t>::get_fip_regions ()
    {
      return this->get_spx<wrapped_t> ()->fip_regions;
    }

    template < class strategy_t >
    typename py_calc_model<strategy_t>::index_array_t
    py_calc_model<strategy_t>::get_rock_regions ()
    {
      return this->get_spx<wrapped_t> ()->rock_regions;
    }

    template < class strategy_t >
    typename py_calc_model<strategy_t>::index_array_t
    py_calc_model<strategy_t>::get_bconn_mainvar ()
    {
      return this->get_spx<wrapped_t> ()->bconn_mainvar;
    }

    template < class strategy_t >
    typename py_calc_model<strategy_t>::item_array_t
    py_calc_model<strategy_t>::get_bconn_pressure ()
    {
      return this->get_spx<wrapped_t> ()->bconn_pressure;
    }

    template < class strategy_t >
    typename py_calc_model<strategy_t>::item_array_t
    py_calc_model<strategy_t>::get_bconn_saturation ()
    {
      return this->get_spx<wrapped_t> ()->bconn_saturation;
    }

    template < class strategy_t >
    typename py_calc_model<strategy_t>::item_array_t
    py_calc_model<strategy_t>::get_bconn_gor ()
    {
      return this->get_spx<wrapped_t> ()->bconn_gor;
    }

    /*******************************************************/

    template < class strategy_t >
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::get_pvt_num ()
    {
      return this->get_spx<wrapped_t> ()->n_pvt_regions;
    }

    template < class strategy_t >
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::get_sat_num ()
    {
      return this->get_spx<wrapped_t> ()->n_sat_regions;
    }

    template < class strategy_t >
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::get_fip_num ()
    {
      return this->get_spx<wrapped_t> ()->n_fip_regions;
    }

    template <typename strategy_t>
    py_scal_3p<strategy_t>
    py_calc_model<strategy_t>::get_scal () const
      {
        return py_scal_3p<strategy_t> (this->get_spx <wrapped_t> ()->scal_prop);
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::py_calc_model_data_t
    py_calc_model<strategy_t>::get_data (index_t index) const
      {
        return py_calc_model_data_t (this->get_spx<wrapped_t> ()->get_data (index));
      }


    /*******************************************************/

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::cmd_num () const
      {
        return (index_t)this->get_spx<wrapped_t> ()->data.size ();
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::all_data_len () const
      {
        return 0;
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::cap_pressure_len () const
      {
        return this->get_spx<wrapped_t> ()->data[0].cap_pressure.size ();
      }
    template <typename strategy_t>
    typename py_calc_model<strategy_t>::item_t
    py_calc_model<strategy_t>::get_cap_pressure (index_t i, index_t j) const
      {
        return this->get_spx<wrapped_t> ()->data[j].cap_pressure[i];
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::s_deriv_cap_pressure_len () const
      {
        return this->get_spx<wrapped_t> ()->data[0].s_deriv_cap_pressure.size ();
      }
    template <typename strategy_t>
    typename py_calc_model<strategy_t>::item_t
    py_calc_model<strategy_t>::get_s_deriv_cap_pressure (index_t i, index_t j) const
      {
        return this->get_spx<wrapped_t> ()->data[j].s_deriv_cap_pressure[i];
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::relative_perm_len () const
      {
        return this->get_spx<wrapped_t> ()->data[0].relative_perm.size ();
      }
    template <typename strategy_t>
    typename py_calc_model<strategy_t>::item_t
    py_calc_model<strategy_t>::get_relative_perm (index_t i, index_t j) const
      {
        return this->get_spx<wrapped_t> ()->data[j].relative_perm[i];
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::s_deriv_relative_perm_len () const
      {
        return this->get_spx<wrapped_t> ()->data[0].s_deriv_relative_perm.size ();
      }
    template <typename strategy_t>
    typename py_calc_model<strategy_t>::item_t
    py_calc_model<strategy_t>::get_s_deriv_relative_perm (index_t i, index_t j) const
      {
        return this->get_spx<wrapped_t> ()->data[j].s_deriv_relative_perm[i];
      }

    //template <typename strategy_t>
    //typename py_calc_model<strategy_t>::index_t
    //py_calc_model<strategy_t>::gas_oil_ratio_len () const {
    //	return 1;
    //}
    //template <typename strategy_t>
    //typename py_calc_model<strategy_t>::item_t
    //py_calc_model<strategy_t>::get_gas_oil_ratio (index_t i, index_t j) const {
    //	return this->get_spx<wrapped_t> ()->gas_oil_ratio[j];
    //}

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::p_deriv_gas_oil_ratio_len () const
      {
        return 1;
      }
    template <typename strategy_t>
    typename py_calc_model<strategy_t>::item_t
    py_calc_model<strategy_t>::get_p_deriv_gas_oil_ratio (index_t /*i*/, index_t j) const
      {
        return this->get_spx<wrapped_t> ()->data[j].p_deriv_gas_oil_ratio;
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::invers_fvf_len () const
      {
        return this->get_spx<wrapped_t> ()->data[0].invers_fvf.size ();
      }
    template <typename strategy_t>
    typename py_calc_model<strategy_t>::item_t
    py_calc_model<strategy_t>::get_invers_fvf (index_t i, index_t j) const
      {
        return this->get_spx<wrapped_t> ()->data[j].invers_fvf[i];
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::p_deriv_invers_fvf_len () const
      {
        return this->get_spx<wrapped_t> ()->data[0].p_deriv_invers_fvf.size ();
      }
    template <typename strategy_t>
    typename py_calc_model<strategy_t>::item_t
    py_calc_model<strategy_t>::get_p_deriv_invers_fvf (index_t i, index_t j) const
      {
        return this->get_spx<wrapped_t> ()->data[j].p_deriv_invers_fvf[i];
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::gor_deriv_invers_fvf_len () const
      {
        return 1;
      }
    template <typename strategy_t>
    typename py_calc_model<strategy_t>::item_t
    py_calc_model<strategy_t>::get_gor_deriv_invers_fvf (index_t /*i*/, index_t j) const
      {
        return this->get_spx<wrapped_t> ()->data[j].gor_deriv_invers_fvf;
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::invers_viscosity_len () const
      {
        return this->get_spx<wrapped_t> ()->data[0].invers_viscosity.size ();
      }
    template <typename strategy_t>
    typename py_calc_model<strategy_t>::item_t
    py_calc_model<strategy_t>::get_invers_viscosity (index_t i, index_t j) const
      {
        return this->get_spx<wrapped_t> ()->data[j].invers_viscosity[i];
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::p_deriv_invers_viscosity_len () const
      {
        return this->get_spx<wrapped_t> ()->data[0].p_deriv_invers_viscosity.size ();
      }
    template <typename strategy_t>
    typename py_calc_model<strategy_t>::item_t
    py_calc_model<strategy_t>::get_p_deriv_invers_viscosity (index_t i, index_t j) const
      {
        return this->get_spx<wrapped_t> ()->data[j].p_deriv_invers_viscosity[i];
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::gor_deriv_invers_viscosity_len () const
      {
        return 1;
      }
    template <typename strategy_t>
    typename py_calc_model<strategy_t>::item_t
    py_calc_model<strategy_t>::get_gor_deriv_invers_viscosity (index_t /*i*/, index_t j) const
      {
        return this->get_spx<wrapped_t> ()->data[j].gor_deriv_invers_viscosity;
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::invers_visc_fvf_len () const
      {
        return this->get_spx<wrapped_t> ()->data[0].invers_visc_fvf.size ();
      }
    template <typename strategy_t>
    typename py_calc_model<strategy_t>::item_t
    py_calc_model<strategy_t>::get_invers_visc_fvf (index_t i, index_t j) const
      {
        return this->get_spx<wrapped_t> ()->data[j].invers_visc_fvf[i];
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::p_deriv_invers_visc_fvf_len () const
      {
        return this->get_spx<wrapped_t> ()->data[0].p_deriv_invers_visc_fvf.size ();
      }
    template <typename strategy_t>
    typename py_calc_model<strategy_t>::item_t
    py_calc_model<strategy_t>::get_p_deriv_invers_visc_fvf (index_t i, index_t j) const
      {
        return this->get_spx<wrapped_t> ()->data[j].p_deriv_invers_visc_fvf[i];
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::gor_deriv_invers_visc_fvf_len () const
      {
        return 1;
      }
    template <typename strategy_t>
    typename py_calc_model<strategy_t>::item_t
    py_calc_model<strategy_t>::get_gor_deriv_invers_visc_fvf (index_t /*i*/, index_t j) const
      {
        return this->get_spx<wrapped_t> ()->data[j].gor_deriv_invers_visc_fvf;
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::density_len () const
      {
        return this->get_spx<wrapped_t> ()->data[0].density.size ();
      }
    template <typename strategy_t>
    typename py_calc_model<strategy_t>::item_t
    py_calc_model<strategy_t>::get_density (index_t i, index_t j) const
      {
        return this->get_spx<wrapped_t> ()->data[j].density[i];
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::p_deriv_density_len () const
      {
        return this->get_spx<wrapped_t> ()->data[0].p_deriv_density.size ();
      }
    template <typename strategy_t>
    typename py_calc_model<strategy_t>::item_t
    py_calc_model<strategy_t>::get_p_deriv_density (index_t i, index_t j) const
      {
        return this->get_spx<wrapped_t> ()->data[j].p_deriv_density[i];
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::gor_deriv_density_len () const
      {
        return 1;
      }
    template <typename strategy_t>
    typename py_calc_model<strategy_t>::item_t
    py_calc_model<strategy_t>::get_gor_deriv_density (index_t /*i*/, index_t j) const
      {
        return this->get_spx<wrapped_t> ()->data[j].gor_deriv_density;
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::porosity_len () const
      {
        return 1;
      }
    template <typename strategy_t>
    typename py_calc_model<strategy_t>::item_t
    py_calc_model<strategy_t>::get_porosity (index_t /*i*/, index_t j) const
      {
        return this->get_spx<wrapped_t> ()->data[j].porosity;
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::p_deriv_porosity_len () const
      {
        return 1;
      }
    template <typename strategy_t>
    typename py_calc_model<strategy_t>::item_t
    py_calc_model<strategy_t>::get_p_deriv_porosity (index_t /*i*/, index_t j) const
      {
        return this->get_spx<wrapped_t> ()->data[j].p_deriv_porosity;
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::truns_mult_len () const
      {
        return 1;
      }
    template <typename strategy_t>
    typename py_calc_model<strategy_t>::item_t
    py_calc_model<strategy_t>::get_truns_mult (index_t /*i*/, index_t j) const
      {
        return this->get_spx<wrapped_t> ()->data[j].truns_mult;
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::p_deriv_truns_mult_len () const
      {
        return 1;
      }
    template <typename strategy_t>
    typename py_calc_model<strategy_t>::item_t
    py_calc_model<strategy_t>::get_p_deriv_truns_mult (index_t /*i*/, index_t j) const
      {
        return this->get_spx<wrapped_t> ()->data[j].p_deriv_truns_mult;
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::mobility_len () const
      {
        return this->get_spx<wrapped_t> ()->data[0].mobility.size ();
      }
    template <typename strategy_t>
    typename py_calc_model<strategy_t>::item_t
    py_calc_model<strategy_t>::get_mobility (index_t i, index_t j) const
      {
        return this->get_spx<wrapped_t> ()->data[j].mobility[i];
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::p_deriv_mobility_len () const
      {
        return this->get_spx<wrapped_t> ()->data[0].p_deriv_mobility.size ();
      }
    template <typename strategy_t>
    typename py_calc_model<strategy_t>::item_t
    py_calc_model<strategy_t>::get_p_deriv_mobility (index_t i, index_t j) const
      {
        return this->get_spx<wrapped_t> ()->data[j].p_deriv_mobility[i];
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::s_deriv_mobility_len () const
      {
        return this->get_spx<wrapped_t> ()->data[0].s_deriv_mobility.size ();
      }
    template <typename strategy_t>
    typename py_calc_model<strategy_t>::item_t
    py_calc_model<strategy_t>::get_s_deriv_mobility (index_t i, index_t j) const
      {
        return this->get_spx<wrapped_t> ()->data[j].s_deriv_mobility[i];
      }

    template <typename strategy_t>
    typename py_calc_model<strategy_t>::index_t
    py_calc_model<strategy_t>::prev_fluid_volume_len () const
      {
        return this->get_spx<wrapped_t> ()->data[0].prev_fluid_volume.size ();
      }
    template <typename strategy_t>
    typename py_calc_model<strategy_t>::item_t
    py_calc_model<strategy_t>::get_prev_fluid_volume (index_t i, index_t j) const
      {
        return this->get_spx<wrapped_t> ()->data[j].prev_fluid_volume[i];
      }

    template <typename strategy_t>
    void py_calc_model<strategy_t>::initialize_datas ()
    {
      int data_size = 3;
      this->get_spx<wrapped_t> () ->data.resize (data_size);

      for (int j = 0; j < data_size; ++j)
        {
          BOSOUT (section::init_data, level::low) << "cap_pressure len = " << cap_pressure_len () << ", initialization: ";
          for (int i = 0; i < cap_pressure_len (); ++i)
            {
              BOSOUT (section::init_data, level::low) << (item_t)(j+i) << " ";
              this->get_spx<wrapped_t> () ->data[j].cap_pressure[i] = (item_t)(j+i);
            }
          BOSOUT (section::init_data, level::low) << bs_end;

          BOSOUT (section::init_data, level::low) << "s_deriv_relative_perm len = " << s_deriv_relative_perm_len () << ", initialization: ";
          for (int i = 0; i < s_deriv_relative_perm_len (); ++i)
            {
              BOSOUT (section::init_data, level::low) << (item_t)(j+i) << " ";
              this->get_spx<wrapped_t> () ->data[j].s_deriv_relative_perm[i] = (item_t)(j+i);
            }
          BOSOUT (section::init_data, level::low) << bs_end;

          BOSOUT (section::init_data, level::low) << "p_deriv_truns_mult len = " << p_deriv_truns_mult_len () << ", initialization: ";
          for (int i = 0; i < p_deriv_truns_mult_len (); ++i)
            {
              BOSOUT (section::init_data, level::low) << (item_t)(j+i) << " ";
              this->get_spx<wrapped_t> () ->data[j].p_deriv_truns_mult = (item_t)(j+i);
            }
          BOSOUT (section::init_data, level::low) << bs_end;
        }
    }

    template <typename strategy_t>
    py_pvt_water <pvt_water <strategy_t> >
    py_calc_model <strategy_t>::get_pvt_water (index_t n_pvt_region) const
    {
      return py_pvt_water <pvt_water <strategy_t> > (get_spx (this)->pvt_water_array [n_pvt_region]);
    }
    template <typename strategy_t>
    py_pvt_gas <pvt_gas <strategy_t> >
    py_calc_model <strategy_t>::get_pvt_gas (index_t n_pvt_region) const
    {
      return py_pvt_gas <pvt_gas <strategy_t> > (get_spx (this)->pvt_gas_array [n_pvt_region]);
    }
    template <typename strategy_t>
    py_pvt_oil <pvt_oil <strategy_t> >
    py_calc_model <strategy_t>::get_pvt_oil (index_t n_pvt_region) const
    {
      return py_pvt_oil <pvt_oil <strategy_t> > (get_spx (this)->pvt_oil_array [n_pvt_region]);
    }

    template <typename strategy_t>
    int
    py_calc_model <strategy_t>::get_n_phases () const
    {
      return get_spx (this)->n_phases;
    }

    template <typename strategy_t>
    bool
    py_calc_model <strategy_t>::is_water () const
    {
      return get_spx (this)->is_water ();
    }
    template <typename strategy_t>
    bool
    py_calc_model <strategy_t>::is_gas () const
    {
      return get_spx (this)->is_gas ();
    }
    template <typename strategy_t>
    bool
    py_calc_model <strategy_t>::is_oil () const
    {
      return get_spx (this)->is_oil ();
    }

    template <typename strategy_t>
    typename py_calc_model <strategy_t>::sp_fi_params_t
    py_calc_model <strategy_t>::get_fi_params () const
    {
      return get_spx (this)->ts_params;
    }

    /*********************************************/



    template <class strategy_t>
    void py_export_calc_model_data_t(const char *name)
    {
      class_< py_calc_model_data<strategy_t> >(name)
      .def("get_cap_pressure",&py_calc_model_data<strategy_t>::get_cap_pressure)
      .def("get_s_deriv_cap_pressure",&py_calc_model_data<strategy_t>::get_s_deriv_cap_pressure)

      .def("get_relative_perm",&py_calc_model_data<strategy_t>::get_relative_perm)
      .def("get_s_deriv_relative_perm",&py_calc_model_data<strategy_t>::get_s_deriv_relative_perm)

      //.def("get_gas_oil_ratio",&py_calc_model_data<strategy_t>::get_gas_oil_ratio)
      .def("get_p_deriv_gas_oil_ratio",&py_calc_model_data<strategy_t>::get_p_deriv_gas_oil_ratio)

      .def("get_invers_fvf",&py_calc_model_data<strategy_t>::get_invers_fvf)
      .def("get_p_deriv_invers_fvf",&py_calc_model_data<strategy_t>::get_p_deriv_invers_fvf)
      .def("get_gor_deriv_invers_fvf",&py_calc_model_data<strategy_t>::get_gor_deriv_invers_fvf)

      .def("get_invers_viscosity",&py_calc_model_data<strategy_t>::get_invers_viscosity)
      .def("get_p_deriv_invers_viscosity",&py_calc_model_data<strategy_t>::get_p_deriv_invers_viscosity)
      .def("get_gor_deriv_invers_viscosity",&py_calc_model_data<strategy_t>::get_gor_deriv_invers_viscosity)

      .def("get_invers_visc_fvf",&py_calc_model_data<strategy_t>::get_invers_visc_fvf)
      .def("get_p_deriv_invers_visc_fvf",&py_calc_model_data<strategy_t>::get_p_deriv_invers_visc_fvf)
      .def("get_gor_deriv_invers_visc_fvf",&py_calc_model_data<strategy_t>::get_gor_deriv_invers_visc_fvf)

      .def("get_density",&py_calc_model_data<strategy_t>::get_density)
      .def("get_p_deriv_density",&py_calc_model_data<strategy_t>::get_p_deriv_density)
      .def("get_gor_deriv_density",&py_calc_model_data<strategy_t>::get_gor_deriv_density)

      .def("get_porosity",&py_calc_model_data<strategy_t>::get_porosity)
      .def("get_p_deriv_porosity",&py_calc_model_data<strategy_t>::get_p_deriv_porosity)

      .def("get_truns_mult",&py_calc_model_data<strategy_t>::get_truns_mult)
      .def("get_p_deriv_truns_mult",&py_calc_model_data<strategy_t>::get_p_deriv_truns_mult)

      .def("get_mobility",&py_calc_model_data<strategy_t>::get_mobility)
      .def("get_p_deriv_mobility",&py_calc_model_data<strategy_t>::get_p_deriv_mobility)
      .def("get_s_deriv_mobility",&py_calc_model_data<strategy_t>::get_s_deriv_mobility)

      .def("get_prev_fluid_volume",&py_calc_model_data<strategy_t>::get_prev_fluid_volume)
      ;
    }

    template <class strategy_t>
    void py_export_calc_model_t(const char *name)
    {
      typedef py_calc_model <strategy_t> py_calc_model_t;

      class_<py_calc_model<strategy_t>, bases<py_bs_node> >(name)
      .def("get_saturation",&py_calc_model<strategy_t>::get_saturation)
      .def("get_pressure",&py_calc_model<strategy_t>::get_pressure)
      .def("get_pvt_num",&py_calc_model<strategy_t>::get_pvt_num)
      .def("get_sat_num",&py_calc_model<strategy_t>::get_sat_num)
      .def("get_fip_num",&py_calc_model<strategy_t>::get_fip_num)
      .def("get_scal", &py_calc_model<strategy_t>::get_scal)
      .add_property ("scal", make_function (&py_calc_model <strategy_t>::get_scal))
      .def("get_plane_flow_rate",&py_calc_model<strategy_t>::get_plane_flow_rate)
      .def("get_full_step_plane_flow_rate",&py_calc_model<strategy_t>::get_full_step_plane_flow_rate)
      .def("get_pvt_regions",&py_calc_model<strategy_t>::get_pvt_regions)
      .def("get_sat_regions",&py_calc_model<strategy_t>::get_sat_regions)
      .def("get_fip_regions",&py_calc_model<strategy_t>::get_fip_regions)
      .def("get_rock_regions",&py_calc_model<strategy_t>::get_rock_regions)
      .def("get_bconn_mainvar",&py_calc_model<strategy_t>::get_bconn_mainvar)
      .def("get_bconn_pressure",&py_calc_model<strategy_t>::get_bconn_pressure)
      .def("get_bconn_saturation",&py_calc_model<strategy_t>::get_bconn_saturation)
      .def("get_bconn_gor",&py_calc_model<strategy_t>::get_bconn_gor)
      .def("get_data",&py_calc_model<strategy_t>::get_data)

      .def("cmd_num",&py_calc_model<strategy_t>::cmd_num)
      .def("all_data_len",&py_calc_model<strategy_t>::all_data_len)

      .def("cap_pressure_len",&py_calc_model<strategy_t>::cap_pressure_len)
      .def("get_cap_pressure",&py_calc_model<strategy_t>::get_cap_pressure)

      .def("s_deriv_cap_pressure_len",&py_calc_model<strategy_t>::s_deriv_cap_pressure_len)
      .def("get_s_deriv_cap_pressure",&py_calc_model<strategy_t>::get_s_deriv_cap_pressure)

      .def("relative_perm_len",&py_calc_model<strategy_t>::relative_perm_len)
      .def("get_relative_perm",&py_calc_model<strategy_t>::get_relative_perm)

      .def("s_deriv_relative_perm_len",&py_calc_model<strategy_t>::s_deriv_relative_perm_len)
      .def("get_s_deriv_relative_perm",&py_calc_model<strategy_t>::get_s_deriv_relative_perm)

      //.def("gas_oil_ratio_len",&py_calc_model<strategy_t>::gas_oil_ratio_len)
      .def("get_gas_oil_ratio",&py_calc_model<strategy_t>::get_gas_oil_ratio)

      .def("p_deriv_gas_oil_ratio_len",&py_calc_model<strategy_t>::p_deriv_gas_oil_ratio_len)
      .def("get_p_deriv_gas_oil_ratio",&py_calc_model<strategy_t>::get_p_deriv_gas_oil_ratio)

      .def("invers_fvf_len",&py_calc_model<strategy_t>::invers_fvf_len)
      .def("get_invers_fvf",&py_calc_model<strategy_t>::get_invers_fvf)

      .def("p_deriv_invers_fvf_len",&py_calc_model<strategy_t>::p_deriv_invers_fvf_len)
      .def("get_p_deriv_invers_fvf",&py_calc_model<strategy_t>::get_p_deriv_invers_fvf)

      .def("gor_deriv_invers_fvf_len",&py_calc_model<strategy_t>::gor_deriv_invers_fvf_len)
      .def("get_gor_deriv_invers_fvf",&py_calc_model<strategy_t>::get_gor_deriv_invers_fvf)

      .def("invers_viscosity_len",&py_calc_model<strategy_t>::invers_viscosity_len)
      .def("get_invers_viscosity",&py_calc_model<strategy_t>::get_invers_viscosity)

      .def("p_deriv_invers_viscosity_len",&py_calc_model<strategy_t>::p_deriv_invers_viscosity_len)
      .def("get_p_deriv_invers_viscosity",&py_calc_model<strategy_t>::get_p_deriv_invers_viscosity)

      .def("gor_deriv_invers_viscosity_len",&py_calc_model<strategy_t>::gor_deriv_invers_viscosity_len)
      .def("get_gor_deriv_invers_viscosity",&py_calc_model<strategy_t>::get_gor_deriv_invers_viscosity)

      .def("invers_visc_fvf_len",&py_calc_model<strategy_t>::invers_visc_fvf_len)
      .def("get_invers_visc_fvf",&py_calc_model<strategy_t>::get_invers_visc_fvf)

      .def("p_deriv_invers_visc_fvf_len",&py_calc_model<strategy_t>::p_deriv_invers_visc_fvf_len)
      .def("get_p_deriv_invers_visc_fvf",&py_calc_model<strategy_t>::get_p_deriv_invers_visc_fvf)

      .def("gor_deriv_invers_visc_fvf_len",&py_calc_model<strategy_t>::gor_deriv_invers_visc_fvf_len)
      .def("get_gor_deriv_invers_visc_fvf",&py_calc_model<strategy_t>::get_gor_deriv_invers_visc_fvf)

      .def("density_len",&py_calc_model<strategy_t>::density_len)
      .def("get_density",&py_calc_model<strategy_t>::get_density)

      .def("p_deriv_density_len",&py_calc_model<strategy_t>::p_deriv_density_len)
      .def("get_p_deriv_density",&py_calc_model<strategy_t>::get_p_deriv_density)

      .def("gor_deriv_density_len",&py_calc_model<strategy_t>::gor_deriv_density_len)
      .def("get_gor_deriv_density",&py_calc_model<strategy_t>::get_gor_deriv_density)

      .def("porosity_len",&py_calc_model<strategy_t>::porosity_len)
      .def("get_porosity",&py_calc_model<strategy_t>::get_porosity)

      .def("p_deriv_porosity_len",&py_calc_model<strategy_t>::p_deriv_porosity_len)
      .def("get_p_deriv_porosity",&py_calc_model<strategy_t>::get_p_deriv_porosity)

      .def("truns_mult_len",&py_calc_model<strategy_t>::truns_mult_len)
      .def("get_truns_mult",&py_calc_model<strategy_t>::get_truns_mult)

      .def("p_deriv_truns_mult_len",&py_calc_model<strategy_t>::p_deriv_truns_mult_len)
      .def("get_p_deriv_truns_mult",&py_calc_model<strategy_t>::get_p_deriv_truns_mult)

      .def("mobility_len",&py_calc_model<strategy_t>::mobility_len)
      .def("get_mobility",&py_calc_model<strategy_t>::get_mobility)

      .def("p_deriv_mobility_len",&py_calc_model<strategy_t>::p_deriv_mobility_len)
      .def("get_p_deriv_mobility",&py_calc_model<strategy_t>::get_p_deriv_mobility)

      .def("s_deriv_mobility_len",&py_calc_model<strategy_t>::s_deriv_mobility_len)
      .def("get_s_deriv_mobility",&py_calc_model<strategy_t>::get_s_deriv_mobility)

      .def("prev_fluid_volume_len",&py_calc_model<strategy_t>::prev_fluid_volume_len)
      .def("get_prev_fluid_volume",&py_calc_model<strategy_t>::get_prev_fluid_volume)
      .def("initialize_datas",&py_calc_model<strategy_t>::initialize_datas)

      .def ("get_pvt_water", &py_calc_model_t::get_pvt_water)
      .def ("get_pvt_gas",   &py_calc_model_t::get_pvt_gas)
      .def ("get_pvt_oil",   &py_calc_model_t::get_pvt_oil)

      .add_property ("n_phases", &py_calc_model_t::get_n_phases)
      .add_property ("is_water", &py_calc_model_t::is_water)
      .add_property ("is_gas", &py_calc_model_t::is_gas)
      .add_property ("is_oil", &py_calc_model_t::is_oil)
      .add_property ("params", &py_calc_model_t::get_fi_params)
      ;
    }

    template <typename T> int get_n_phases  (T *t) { return t->n_phases; }
    template <typename T> bool get_is_water (T *t) { return t->is_water (); }
    template <typename T> bool get_is_gas   (T *t) { return t->is_gas (); }
    template <typename T> bool get_is_oil   (T *t) { return t->is_oil (); }

    template <typename T>
    smart_ptr <named_pbase, true>
    get_params (T *t)
    {
      return t->ts_params;
    }

    PY_EXPORTER (calc_model_exporter, default_exporter)
      .add_property ("n_phases",  get_n_phases <T>)
      .add_property ("is_water",  get_is_water <T>)
      .add_property ("is_gas",    get_is_gas <T>)
      .add_property ("is_oil",    get_is_oil <T>)
      .add_property ("params",    get_params <T>)
    PY_EXPORTER_END;

    void py_export_calc_model()
    {
      py_export_calc_model_data_t<base_strategy_di>("calc_model_data_d");
      py_export_calc_model_data_t<base_strategy_fi>("calc_model_data_f");
      py_export_calc_model_data_t<base_strategy_mixi>("calc_model_data_m");

      strategy_exporter::export_base <calc_model, calc_model_exporter> ("calc_model");
      //py_export_calc_model_t<base_strategy_di>("calc_model_d");
      //py_export_calc_model_t<base_strategy_fi>("calc_model_f");
      //py_export_calc_model_t<base_strategy_mixi>("calc_model_m");
    }

    template class py_calc_model<base_strategy_fi>;
    template class py_calc_model<base_strategy_di>;
    template class py_calc_model<base_strategy_mixi>;
  } //ns python
} //ns bs

