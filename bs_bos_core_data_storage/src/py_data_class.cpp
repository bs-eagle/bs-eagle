#include "bs_bos_core_data_storage_stdafx.h"

#include "py_data_class.h"
#include "numpy_tools.h"
#include "py_object_handler.h"

#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>

#include "export_python_wrapper.h"

#ifdef BSPY_EXPORTING_PLUGIN
using namespace boost::python;
namespace bp = boost::python;
using namespace blue_sky::tools;

namespace blue_sky {
namespace python {

  //template <typename strategy_t, typename strategy_t::i_type_t i>
  //struct idata_int_array_getter
  //{
  //  static 
  //  array_uint8_t
  //  get_array (py_idata  *self)
  //  {
  //    return self->get_int_array (i);
  //  }
  //  static 
  //  void 
  //  set_array (py_idata  *self, const boost::python::object &obj)
  //  {
  //    self->set_int_array (i, obj);
  //  }
  //};
  //template <typename strategy_t, typename strategy_t::i_type_t i>
  //struct idata_float_array_getter
  //{
  //  static 
  //  array_float16_t
  //  get_array (py_idata  *self)
  //  {
  //    return self->get_float_array (i);
  //  }
  //  static 
  //  void 
  //  set_array (py_idata  *self, const boost::python::object &obj)
  //  {
  //    self->set_float_array (i, obj);
  //  }
  //};

#define REGISTER_PROPERTIES_I(z, n_, t_)                                              \
  BOOST_PP_TUPLE_ELEM (4, 0, t_).add_property (BOOST_PP_TUPLE_ELEM (4, 1, t_)[n_],    \
    &BOOST_PP_TUPLE_ELEM (4, 2, t_) <BOOST_PP_TUPLE_ELEM (4, 3, t_), n_>::get_array,  \
    &BOOST_PP_TUPLE_ELEM (4, 2, t_) <BOOST_PP_TUPLE_ELEM (4, 3, t_), n_>::set_array);

#define REGISTER_PROPERTIES(from_, to_, t_) \
  BOOST_PP_REPEAT_FROM_TO (from_, to_, REGISTER_PROPERTIES_I, t_)

  template <typename T>
  typename T::vec_i &
  get_equil_regions (T *t)
  {
    return t->equil_regions;
  }

  template <typename T>
  int
  get_units_in (T *t)
  {
    return t->input_units_converter.get_input_units ();
  }
  template <typename T>
  int
  get_units_out (T *t)
  {
    return t->output_units_converter.get_output_units ();
  }
  template <typename T>
  void
  set_units_in (T *t, int v)
  {
    t->input_units_converter.set_input_units (v);
  }
  template <typename T>
  void
  set_units_out (T *t, int v)
  {
    t->output_units_converter.set_output_units (v);
  }
/*
  template <typename T>
  void 
  set_int_array (T *t, int cur_index, typename T::sp_arr_i i_array)
  {
    (*t->i_map)[cur_index].array = i_array;
  }
  template <typename T>
  void 
  set_fp_array (T *t, int cur_index, typename T::sp_arr_fp fp_array)
  {
    (*t->fp_map)[cur_index].array = fp_array;
  }

  template <typename T>
  typename T::sp_arr_i
  get_int_array (T *t, int cur_index)
  {
    return t->get_int_non_empty_array (cur_index);
  }
  template <typename T>
  typename T::sp_arr_fp
  get_fp_array (T *t, int cur_index)
  {
    return t->get_fp_non_empty_array (cur_index);
  }

  template <typename T>
  long long
  get_nx (T *t)
  {
    return t->dimens.nx;
  }
  template <typename T>
  long long
  get_ny (T *t)
  {
    return t->dimens.ny;
  }
  template <typename T>
  long long
  get_nz (T *t)
  {
    return t->dimens.nz;
  }
*/

  PY_EXPORTER (idata_exporter, default_exporter)
/*
    .def("init",                          &T::init)
    .add_property ("rpo_model",           &T::rpo_model)
    .add_property ("minimal_pore_volume", &T::minimal_pore_volume)
    .add_property ("nx",                  get_nx <T>)
    .add_property ("ny",                  get_ny <T>)
    .add_property ("nz",                  get_nz <T>)
    .add_property ("dimens",              &T::dimens)
    .add_property ("pvt_region",          &T::pvt_region)
    .add_property ("sat_region",          &T::sat_region)
    .add_property ("eql_region",          &T::eql_region)
    .add_property ("fip_region",          &T::fip_region)
    .add_property ("fi_n_phase",          &T::fi_n_phase)
    .add_property ("fi_phases",           &T::fi_phases)
    .add_property ("rock_region",         &T::rock_region)
    //.add_property ("equil_regions",       make_function (get_equil_regions <T>, return_value_policy<reference_existing_object> ()))
    //.add_property ("units_in",            get_units_in <T>,  set_units_in <T>)
    //.add_property ("units_out",           get_units_out <T>, set_units_out <T>)
    .add_property ("title",               &T::title)
    .add_property ("int_array",           get_int_array   <T>, set_int_array <T>)
    .add_property ("fp_array",            get_fp_array <T>, get_fp_array <T>)
    .def ("get_int_array",                get_int_array   <T>)
    .def ("get_fp_array",                 get_fp_array <T>)
    .def ("set_int_array",                set_int_array   <T>)
    .def ("set_fp_array",                 set_fp_array <T>)
    //.def ("get_float_array",              &T::get_float_array)
    //.def ("get_int_array",                &T::get_int_array)
*/    
  PY_EXPORTER_END;

    //
    //void py_export_idata_t(const char *name)
    //{
    //  typedef py_idata  T;

    //  class_<T, bases<py_bs_node> > idata_export(name);

    //  idata_export
    //    .def("init",&T::init)
    //    .add_property ("rpo_model",           &T::get_rpo_model, &T::set_rpo_model)
    //    .add_property ("minimal_pore_volume", &T::get_minimal_pore_volume, &T::set_minimal_pore_volume)
    //    .add_property ("restart",             &T::get_restart, &T::set_restart)
    //    .add_property ("nx",                  &T::get_nx, &T::set_nx)
    //    .add_property ("ny",                  &T::get_ny, &T::set_ny)
    //    .add_property ("nz",                  &T::get_nz, &T::set_nz)
    //    .add_property ("pvt_region",          &T::get_pvt_region, &T::set_pvt_region)
    //    .add_property ("sat_region",          &T::get_sat_region, &T::set_sat_region)
    //    .add_property ("eql_region",          &T::get_eql_region, &T::set_eql_region)
    //    .add_property ("fip_region",          &T::get_fip_region, &T::set_fip_region)
    //    .add_property ("fi_n_phase",          &T::get_fi_n_phase, &T::set_fi_n_phase)
    //    .add_property ("fi_phases",           &T::get_fi_phases, &T::set_fi_phases)
    //    .add_property ("rock_region",         &T::get_rock_region, &T::set_rock_region)
    //    .add_property ("equil_regions",       make_function (&T::get_equil_regions, return_value_policy<copy_const_reference> ()))
    //    .add_property ("units_in",            &T::get_units_in, &T::set_units_in)
    //    .add_property ("units_out",           &T::get_units_out, &T::set_units_out)
    //    .add_property ("title",               &T::get_title, &T::set_title)
    //    .add_property ("int_array",           &T::get_int_array, &T::set_int_array)
    //    .add_property ("float_array",         &T::get_float_array, &T::set_float_array)
    //    .def ("get_int_array",                &T::get_int_array)
    //    .def ("get_float_array",              &T::get_float_array)
    //    .def ("set_int_array",                &T::set_int_array)
    //    .def ("set_float_array",              &T::set_float_array)
    //    .def ("__getitem__",                  &T::get_float_array_by_name)
    //    ;

    //  REGISTER_PROPERTIES (0, 9,  (idata_export, int_names_table,    idata_int_array_getter,   strategy_t));
    //  REGISTER_PROPERTIES (0, 16, (idata_export, double_names_table, idata_float_array_getter, strategy_t));
    //}


  void py_export_idata()
  {
    //strategy_exporter::export_base<idata, idata_exporter> ("idata");
  }

} // namespace python
} // namespace blue_sky
#endif