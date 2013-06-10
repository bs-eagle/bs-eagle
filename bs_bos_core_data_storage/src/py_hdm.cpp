#include "bs_bos_core_data_storage_stdafx.h"

#include "py_hdm.h"
#include "hdm.h"
#include "py_data_class.h"

#include "export_python_wrapper.h"

#ifdef BSPY_EXPORTING_PLUGIN
using namespace boost::python;

namespace blue_sky {
namespace python {

  PY_EXPORTER (hdm_exporter, default_exporter)
    //.add_property ("data", get_data <T, strategy_t>)
    .def("init_fluids", &T::init_fluids)
    .def("init_equil", &T::init_equil)
    .def("init", &T::init)
    .def("read_keyword_file", &T::read_keyword_file)
    .def("get_pool", &T::get_pool)
    .def("get_reader", &T::get_reader)
    .def("get_mesh", &T::get_mesh)
    .def("get_scal", &T::get_scal)
    .def("get_pvt", &T::get_pvt)
    .def("get_prop", &T::get_prop)
    .def("get_proc_params", &T::get_proc_params)
    .def("get_well_pool", &T::get_well_pool)
    .def("get_keyword_manager", &T::get_keyword_manager)
    .def("get_equil_model", &T::get_equil_model)
    .def("set_comp_ref_pressure", &T::set_comp_ref_pressure)
    .def("get_comp_ref_pressure", &T::get_comp_ref_pressure)
    .def("set_comp_const", &T::set_comp_const)
    .def("get_comp_const", &T::get_comp_const)
    .def_readwrite("scal", &T::scal_3p_)
    .def_readwrite("pvt", &T::pvt_3p_)
    .def_readwrite("init_model", &T::init_model_)
    .def_readwrite("event_manager", &T::event_manager_)
    .def_readwrite("well_pool", &T::well_pool_)
    .def_readwrite("equil_model", &T::equil_model_)
    .def_readwrite("mesh", &T::mesh)
    .def_readwrite("keyword_manager", &T::km)
    .def_readwrite("reader", &T::reader)
    .def_readwrite("data", &T::data)
  PY_EXPORTER_END;

  void py_export_hdm ()
  {
    base_exporter<hdm_iface, empty_exporter>::export_class ("hdm_iface");
    class_exporter <hdm, hdm_iface, hdm_exporter>::export_class ("hdm");
  }

} // namespace python
} // namespace blue_sky
#endif
