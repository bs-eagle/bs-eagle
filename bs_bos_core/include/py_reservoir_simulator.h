/**
 *       \file  py_reservoir_simulator.h
 *      \brief  Exports python wrappers for reservoir_simulator,
 *              see reservoir_simulator.h
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  18.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef PY_RESERVOIR_SIMULATOR_H
#define PY_RESERVOIR_SIMULATOR_H

#ifdef BSPY_EXPORTING_PLUGIN
#include "py_event_manager.h"
#include "py_calc_model.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "py_data_manager.h"
#include "py_jacobian_matrix.h"
#include "py_linear_solvers.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky {
namespace python {

  /**
   * \brief  Exports wrappers to python
   * */
  void 
  py_export_reservoir_simulator ();

} // namespace python
} // namespace blue_sky

#endif
#endif // PY_RESERVOIR_SIMULATOR_H
