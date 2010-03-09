/**
 *       \file  py_reservoir.h
 *      \brief  Export python wrappers for reservoir,
 *              see reservoir.h
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  17.10.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef PY_RESERVOIR_H
#define PY_RESERVOIR_H

#ifdef BSPY_EXPORTING_PLUGIN
#include "reservoir.h"
#include "facility_manager.h"
#include "well_connection.h"

namespace blue_sky {
namespace python {

  /**
   * \brief  Export wrappers to python
   * */
  void 
  py_export_reservoir();

} // namespace python
} // namespace blue_sky

#endif
#endif // PY_RESERVOIR_H

