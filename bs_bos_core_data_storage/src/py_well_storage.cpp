#include "bs_bos_core_data_storage_stdafx.h"

#include "py_well_storage.h"
#include "well_storage.h"

#include "export_python_wrapper.h"

#ifdef BSPY_EXPORTING_PLUGIN
using namespace boost::python;

namespace blue_sky {
namespace python {

  PY_EXPORTER (params_exporter, default_exporter)
  PY_EXPORTER_END;
  
  PY_EXPORTER (ws_exporter, default_exporter)
    .def("add_group", &T::add_group)
    .def("add_well", &T::add_well)
    .def("add_conng", &T::add_conng)
    .def("get_well_dates", &T::get_well_dates)
    .def("get_well_fvalues", &T::get_well_fvalues)
    .def("get_well_ivalues", &T::get_well_ivalues)
    .def("add_well_to_group", &T::add_well_to_group)
    .def("set_well_fparam", &T::set_well_fparam)
    .def("set_well_iparam", &T::set_well_iparam)
    .def("set_well_fparams", &T::set_well_fparams)
    .def("set_well_iparams", &T::set_well_iparams)
    .def("get_group_dates", &T::get_group_dates)
    .def("get_group_fvalues", &T::get_group_fvalues)
    .def("get_group_ivalues", &T::get_group_ivalues)
  PY_EXPORTER_END;

  void py_export_well_storage ()
  {
    base_exporter<well_storage_iface, empty_exporter>::export_class ("well_storage_iface");
    class_exporter <well_storage, well_storage_iface, ws_exporter>::export_class ("well_storage");
    //base_exporter<well_params, params_exporter>::export_class ("well_params");
  }

} // namespace python
} // namespace blue_sky
#endif