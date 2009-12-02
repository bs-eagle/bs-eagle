/**
 *       \file  py_calc_model.cpp
 *      \brief  Python wrappers for calc_model, calc_model_data
 *     \author  Sergey Miryanov
 *       \date  02.12.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "boost_array_adapter.h"
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
      .add_property ("n_phases",    get_n_phases <T>)
      .add_property ("is_water",    get_is_water <T>)
      .add_property ("is_gas",      get_is_gas <T>)
      .add_property ("is_oil",      get_is_oil <T>)
      .add_property ("params",      get_params <T>)
      .add_property ("data",        &T::data)
      .add_property ("saturation",  &T::saturation_3p)
      .add_property ("pressure",    &T::pressure)
      .add_property ("pvt_num",     &T::n_pvt_regions)
      .add_property ("sat_num",     &T::n_sat_regions)
      .add_property ("fip_num",     &T::n_fip_regions)
      .add_property ("scal",        &T::scal_prop)
      .add_property ("pvt_regions", &T::pvt_regions)
      .add_property ("sat_regions", &T::sat_regions)
      .add_property ("fip_regions", &T::fip_regions)
      .add_property ("rock_regions",&T::rock_regions)
      .add_property ("pvt_water",   &T::pvt_water_array)
      .add_property ("pvt_gas",     &T::pvt_gas_array)
      .add_property ("pvt_oil",     &T::pvt_oil_array)
    PY_EXPORTER_END;

    PY_EXPORTER (calc_model_data_exporter, empty_exporter)
      .add_property ("cap_pressure",              detail::boost_array_adapter (&T::cap_pressure))
      .add_property ("s_deriv_cap_pressure",      detail::boost_array_adapter (&T::s_deriv_cap_pressure))
      .add_property ("relative_perm",             detail::boost_array_adapter (&T::relative_perm))
      .add_property ("s_deriv_relative_perm",     detail::boost_array_adapter (&T::s_deriv_relative_perm))
      .add_property ("p_deriv_gas_oil_ratio",     &T::p_deriv_gas_oil_ratio)
      .add_property ("invers_fvf",                detail::boost_array_adapter (&T::invers_fvf))
      .add_property ("p_deriv_invers_fvf",        detail::boost_array_adapter (&T::p_deriv_invers_fvf))
      .add_property ("gor_deriv_invers_fvf",      &T::gor_deriv_invers_fvf)
      .add_property ("invers_visc",               detail::boost_array_adapter (&T::invers_viscosity))
      .add_property ("p_deriv_invers_visc",       detail::boost_array_adapter (&T::p_deriv_invers_viscosity))
      .add_property ("gor_deriv_invers_visc",     &T::gor_deriv_invers_viscosity)
      .add_property ("invers_visc_fvf",           detail::boost_array_adapter (&T::invers_visc_fvf))
      .add_property ("p_deriv_invers_visc_fvf",   detail::boost_array_adapter (&T::p_deriv_invers_visc_fvf))
      .add_property ("gor_deriv_invers_visc_fvf", &T::gor_deriv_invers_visc_fvf)
      .add_property ("density",                   detail::boost_array_adapter (&T::density))
      .add_property ("p_deriv_density",           detail::boost_array_adapter (&T::p_deriv_density))
      .add_property ("gor_deriv_density",         &T::gor_deriv_density)
      .add_property ("porosity",                  &T::porosity)
      .add_property ("p_deriv_porosity",          &T::p_deriv_porosity)
      .add_property ("truns_mult",                &T::truns_mult)
      .add_property ("p_deriv_truns_mult",        &T::p_deriv_truns_mult)
      .add_property ("mobility",                  detail::boost_array_adapter (&T::mobility))
      .add_property ("p_deriv_mobility",          detail::boost_array_adapter (&T::p_deriv_mobility))
      .add_property ("s_deriv_mobility",          detail::boost_array_adapter (&T::s_deriv_mobility))
      .add_property ("prev_fluid_volume",         detail::boost_array_adapter (&T::prev_fluid_volume))
    PY_EXPORTER_END;

    template <typename T>
    void
    export_calc_model_data_vector (const char *name)
    {
      class_ <T> (name)
        .def (vector_indexing_suite <T> ())
        ;
    }

    void py_export_calc_model()
    {
      strategy_exporter::export_base_ext <calc_model_data, calc_model_data_exporter, class_type::concrete_class> ("calc_model_data");

      export_calc_model_data_vector <shared_vector <calc_model_data <base_strategy_di> > >   ("calc_model_data_vector_di");
      export_calc_model_data_vector <shared_vector <calc_model_data <base_strategy_fi> > >   ("calc_model_data_vector_fi");
      export_calc_model_data_vector <shared_vector <calc_model_data <base_strategy_mixi> > > ("calc_model_data_vector_mixi");

      strategy_exporter::export_base <calc_model, calc_model_exporter> ("calc_model");
    }

  } //ns python
} //ns bs

