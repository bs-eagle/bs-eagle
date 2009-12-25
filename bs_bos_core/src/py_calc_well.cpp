/**
 *       \file  py_calc_well.cpp
 *      \brief  Python wrapper for calc_well, see calc_well.h
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  23.06.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  A bit outdate
 * */
#include "stdafx.h"

#ifdef BSPY_EXPORTING_PLUGIN

#include "py_calc_well.h"
#include "reservoir.h"
#include "facility_manager.h"

namespace blue_sky {
namespace python {

  void
  py_export_calc_well ()
  {
    strategy_exporter::export_base <well, well_exporter> ("well_seq");
    strategy_exporter::export_base <wells::connection, connection_exporter> ("connection_seq");

    strategy_exporter::export_class <py_well, well, well_exporter> ("py_well_seq");
    strategy_exporter::export_class <py_connection, wells::connection, connection_exporter> ("py_connection_seq");
  }

#endif  // #ifdef BSPY_EXPORTING_PLUGIN
} // namespace python
} // namespace blue_sky

