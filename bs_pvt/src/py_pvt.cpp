/**
 *       \file  py_pvt.cpp
 *      \brief  python wrappers for pvt_
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  08.05.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "bs_pvt_stdafx.h"

#include "py_pvt.h"

#ifdef BSPY_EXPORTING_PLUGIN

#include "pvt_base.h"
#include "pvt_dead_oil.h"
#include "pvt_gas.h"
#include "pvt_oil.h"
#include "pvt_water.h"

#include "export_python_wrapper.h"

namespace blue_sky  {
namespace python    {

  PY_EXPORTER (pvt_exporter, default_exporter)
    .add_property ("p_step",          &T::get_p_step)
    .add_property ("surface_density", &T::get_surface_density)
    .def ("get_gor_for_pressure",     &T::get_gor_for_pressure)
    .def ("print",                    &T::print)
  PY_EXPORTER_END;

  void
  py_export_pvt ()
  {
    using namespace boost::python;

    base_exporter<pvt_base, pvt_exporter>::export_class ("pvt_base");
    class_exporter<pvt_oil, pvt_base, pvt_exporter>::export_class  ("pvt_oil");
    class_exporter<pvt_gas, pvt_base, pvt_exporter>::export_class ("pvt_gas");
    class_exporter<pvt_water, pvt_base, pvt_exporter>::export_class ("pvt_water");
    class_exporter<pvt_dead_oil, pvt_base, pvt_exporter>::export_class ("pvt_dead_oil");

  }

} // namespace python
} // namespace blue_sky

#endif  // #ifdef BSPY_EXPORTING_PLUGIN

