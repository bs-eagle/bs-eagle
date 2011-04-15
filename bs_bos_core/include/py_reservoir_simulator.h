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
