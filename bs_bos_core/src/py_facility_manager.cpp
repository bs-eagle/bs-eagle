#include "stdafx.h"

#include "py_facility_manager.h"
#include "reservoir.h"
#include "facility_manager.h"

#include "export_python_wrapper.h"

using namespace boost::python;

namespace blue_sky {
namespace python {

  //template <class strategy_t>
  //void py_export_facility_manager_t(const char *name)
  //{
  //  class_<py_facility_manager<strategy_t>, bases<py_objbase> >(name)
  //  //.def("__iter__", range(&py_facility_manager<strategy_t>::wells_begin,
  //  //									 	 	 &py_facility_manager<strategy_t>::wells_end))
  //  .def("add_well", &py_facility_manager<strategy_t>::add_well)
  //  .def("wells", &py_facility_manager<strategy_t>::wells_begin)
  //  .def("wells_end", &py_facility_manager<strategy_t>::wells_end)
  //  //.def("__iter__", pass_throught) // boost::python::iterator< py_facility_manager<strategy_t> >())
  //  //.def("next", &py_facility_manager<strategy_t>::next)
  //  ;
  //}

  PY_EXPORTER (facility_manager_exporter, default_exporter)
    .def ("add_well",   &T::add_well)
    .def ("wells",      &T::wells_begin)
    .def ("wells_end",  &T::wells_end)
  PY_EXPORTER_END;

  void py_export_facility_manager()
  {
    strategy_exporter::export_base <facility_manager, facility_manager_exporter> ("facility_manager");
  }

} // namespace python
} // namespace blue_sky
