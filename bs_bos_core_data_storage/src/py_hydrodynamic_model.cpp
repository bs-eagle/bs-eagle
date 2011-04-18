#include "bs_bos_core_data_storage_stdafx.h"

#include "py_hydrodynamic_model.h"
#include "hydrodynamic_model.h"
#include "py_data_class.h"

#include "export_python_wrapper.h"

#ifdef BSPY_EXPORTING_PLUGIN
using namespace boost::python;

namespace blue_sky {
namespace python {

  PY_EXPORTER (hydrodynamic_model_exporter, default_exporter)
    //.add_property ("data", get_data <T, strategy_t>)
    .def("init", &T::init)
    .def("read_keyword_file", &T::read_keyword_file)
    .def("get_data", &T::get_data)
    .def("get_reader", &T::get_reader)
    .def("get_mesh", &T::get_mesh)
    .def("get_keyword_manager", &T::get_keyword_manager)
  PY_EXPORTER_END;

  void py_export_hydrodynamic_model ()
  {
    base_exporter<hydrodynamic_model_iface, empty_exporter>::export_class ("hydrodynamic_model_iface");
    class_exporter <hydrodynamic_model, hydrodynamic_model_iface, hydrodynamic_model_exporter>::export_class ("hydrodynamic_model");
  }

} // namespace python
} // namespace blue_sky
#endif