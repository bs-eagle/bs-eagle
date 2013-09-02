/**
 * @file py_dt_tools.h
 * @brief python wraper for #dt_tools
 * @author Oleg Borschuk
 * @version
 * @date 2012-03-01
 */
#ifndef PY_DT_TOOLS_DEU83ZYI

#define PY_DT_TOOLS_DEU83ZYI

#include <string>
#include "dt_tools_iface.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "dummy_base.h"
#include "construct_python_object.h"
#include "make_me_happy.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include "export_python_wrapper.h"

#ifdef BSPY_EXPORTING_PLUGIN
namespace blue_sky
  {
  namespace python
    {

  PY_EXPORTER (py_dt_tools_exporter, default_exporter)
    .def ("d2date",                             &T::d2date,
        args ("d"), "convert d to list (year, month, day, hour, minute, second)")
    .def ("date2d",                             &T::date2d,
        args ("year", "month", "day", "hour", "minute", "second"), "date and time to double")
    .def ("d2str",                              &T::d2str,
        args ("d"), "convert date from double to str")
    .def ("t2str",                              &T::t2str,
        args ("d"), "convert time from double to str")
    .def ("__str__",                            &T::py_str)
  PY_EXPORTER_END;

  //! export matrices to python
  void py_export_dt_tools ();


  } // namespace python
} // namespace blue_sky

#endif // #ifdef BSPY_EXPORTING_PLUGIN

#endif /* end of include guard: PY_DT_TOOLS_DEU83ZYI */
