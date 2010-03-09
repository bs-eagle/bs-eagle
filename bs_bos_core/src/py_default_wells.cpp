/**
 *       \file  py_default_wells.cpp
 *      \brief  Exports python wrappers for default_well
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  22.05.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "stdafx.h"
#include "py_default_wells.h"
#include "export_python_wrapper.h"
#include "default_well.h"
#include "default_connection.h"
#include "py_calc_well.h"

#ifdef BSPY_EXPORTING_PLUGIN
using namespace blue_sky::wells;

namespace blue_sky {
namespace python {

  void 
  py_export_default_wells ()
  {
    strategy_exporter::export_class <default_well, well, well_exporter> ("default_well");
    strategy_exporter::export_class <default_connection, blue_sky::wells::connection, connection_exporter> ("default_connection");
  }

} // namespace python
} // namespace blue_sky
#endif
