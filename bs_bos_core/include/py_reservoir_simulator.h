#ifndef PY_RESERVOIR_SIMULATOR_H
#define PY_RESERVOIR_SIMULATOR_H

#include "py_event_manager.h"
#include "py_calc_model.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "py_data_manager.h"
#include "py_jacobian_matrix.h"
#include "py_linear_solvers.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky {
namespace python {

  void 
  py_export_reservoir_simulator ();

} // namespace python
} // namespace blue_sky

#endif // PY_RESERVOIR_SIMULATOR_H
