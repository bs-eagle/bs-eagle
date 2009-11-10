#include "stdafx.h"

#include "py_reservoir.h"
#include "reservoir.h"
#include "well_connection.h"
#include "py_facility_manager.h"

using namespace boost::python;

namespace blue_sky {
namespace python {

  PY_EXPORTER (reservoir_exporter, default_exporter)
    .def ("add_filter_well",                        &T::add_filter_well)
    .add_property ("facility_list",                 &T::get_facility_list)
    .add_property ("well_factory",                  make_function (&T::get_well_factory),                 make_function (&T::set_well_factory))
    .add_property ("well_controller_factory",       make_function (&T::get_well_controller_factory),      make_function (&T::set_well_controller_factory))
    .add_property ("well_limit_operation_factory",  make_function (&T::get_well_limit_operation_factory), make_function (&T::set_well_limit_operation_factory))
  PY_EXPORTER_END;

  void py_export_reservoir()
  {
    strategy_exporter::export_base <reservoir, reservoir_exporter> ("reservoir");
  }

} // namespace python
} // namespace blue_sky
