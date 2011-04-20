/**
 *       \file  py_calc_model.cpp
 *      \brief  Python wrappers for calc_model, calc_model_data
 *     \author  Sergey Miryanov
 *       \date  02.12.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"

#ifdef BSPY_EXPORTING_PLUGIN
#include "boost_array_adapter.h"
#include "calc_model.h"
#include "py_calc_model.h"
#include BS_FORCE_PLUGIN_IMPORT ()
#include "py_scal_wrapper.h"
#include "py_data_class.h"
#include BS_STOP_PLUGIN_IMPORT ()

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

    template <typename T>
    smart_ptr <typename T::scal_3p_t, true>
    get_scal_3p (T *t)
    {
      return t->scal_prop;
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
      .add_property ("scal",        get_scal_3p <T>)
      .add_property ("pvt_regions", &T::pvt_regions)
      .add_property ("sat_regions", &T::sat_regions)
      .add_property ("fip_regions", &T::fip_regions)
      .add_property ("rock_regions",&T::rock_regions)
      .add_property ("pvt_water",   &T::pvt_water_array)
      .add_property ("pvt_gas",     &T::pvt_gas_array)
      .add_property ("pvt_oil",     &T::pvt_oil_array)
    PY_EXPORTER_END;

    PY_EXPORTER (calc_model_data_exporter, empty_exporter)
#if 0
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
#endif //0
    PY_EXPORTER_END;

    template <typename T>
    void
    export_calc_model_data_vector (const char *name)
    {
      class_ <T> (name)
        .def (vector_indexing_suite <T> ())
        ;
    }

    template <typename T>
    typename T::value_type
    get_pvt_item (T *t, size_t index)
    {
      return t->operator[] (index);
    }

    template <typename T>
    struct pvt_array_iterator
    {
      pvt_array_iterator (T *t)
      : array_ (t)
      , iterator_ (t->begin ())
      , iterator_end_ (t->end ())
      {
      }

      typename T::value_type
      next ()
      {
#ifdef _DEBUG
        if (iterator_end_ != array_->end ())
          {
            bs_throw_exception ("PVT array iterator not more valid");
          }
#endif

        if (iterator_ == iterator_end_)
          {
            PyErr_SetString (PyExc_StopIteration, "No more data");
            boost::python::throw_error_already_set ();
          }

        return *(iterator_++);
      }

      T                     *array_;
      typename T::iterator  iterator_;
      typename T::iterator  iterator_end_;
    };

    template <typename T>
    pvt_array_iterator <T>
    get_pvt_iterator (T *t)
    {
      return pvt_array_iterator <T> (t);
    }

    template <typename T>
    void
    export_pvt_iterator (const char *name)
    {
      using namespace boost::python;

      class_ <pvt_array_iterator <T> > (name, init <T *> ())
        .def ("next",     &pvt_array_iterator <T>::next)
        .def ("__iter__", pass_through)
        ;
    }

    template <typename T>
    void
    export_pvt_array (const char *name)
    {
      typedef std::vector <smart_ptr <T> > T_v;

      class_ <T_v> (name)
        .def ("__getitem__",  get_pvt_item <T_v>)
        .def ("__len__",      &T_v::size)
        .def ("__iter__",     get_pvt_iterator <T_v>)
        ;

      export_pvt_iterator <T_v> ("pvt_array_iter");
    }

    void py_export_calc_model()
    {
      base_exporter_nobs <calc_model_data, calc_model_data_exporter, class_type::concrete_class>::export_class ("calc_model_data");
      export_calc_model_data_vector <shared_vector <calc_model_data> >   ("calc_model_data_vector");

      export_pvt_array <pvt_dead_oil> ("pvt_dead_oil_array");
      export_pvt_array <pvt_water> ("pvt_water_array");
      export_pvt_array <pvt_gas> ("pvt_gas");

      base_exporter <calc_model, calc_model_exporter>::export_class ("calc_model");
    }

  } //ns python
} //ns bs
#endif
