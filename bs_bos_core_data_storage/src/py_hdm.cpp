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
    .def("init", &T::init)
    .def("read_keyword_file", &T::read_keyword_file)
    .def("get_pool", &T::get_pool)
    .def("get_reader", &T::get_reader)
    .def("get_mesh", &T::get_mesh)
    .def("get_scal", &T::get_scal)
    .def("get_pvt", &T::get_pvt)
    .def("get_prop", &T::get_prop)
    .def("get_keyword_manager", &T::get_keyword_manager)
  PY_EXPORTER_END;

  void py_export_hdm ()
  {
    base_exporter<hdm_iface, empty_exporter>::export_class ("hdm_iface");
    class_exporter <hdm, hdm_iface, hdm_exporter>::export_class ("hdm");
  }

} // namespace python
} // namespace blue_sky
#endif