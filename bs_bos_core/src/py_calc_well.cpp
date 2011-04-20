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
#include "py_calc_well.h"
#include "reservoir.h"
#include "facility_manager.h"

#ifdef BSPY_EXPORTING_PLUGIN
namespace blue_sky {
namespace python {

  void
  py_export_calc_well ()
  {
    base_exporter <well, well_exporter>::export_class ("well_seq");
    base_exporter <wells::connection, connection_exporter>::export_class ("connection_seq");

    class_exporter <py_well, well, well_exporter>::export_class ("py_well_seq");
    class_exporter <py_connection, wells::connection, connection_exporter>::export_class ("py_connection_seq");
  }
} // namespace python
} // namespace blue_sky
#endif  // #ifdef BSPY_EXPORTING_PLUGIN
