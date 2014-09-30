/**
 *       \file  py_event_base.cpp
 *      \brief  Python wrappers for event_base, python events iterator,
 *              for event_base see event_base.h
 *     \author  Nikonov Max
 *       \date  17.10.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  A bit outdate
 * */


#ifdef BSPY_EXPORTING_PLUGIN
#include "calc_model.h"
#include "reservoir.h"

#include "py_event_base.h"

using namespace boost::python;

namespace blue_sky  {
namespace python    {

  PY_EXPORTER (event_base_exporter, default_exporter)
    .def ("apply", &T::apply)
  PY_EXPORTER_END;

  void py_export_events ()
  {
    base_exporter <event_base, event_base_exporter>::export_class ("event_base");
    class_exporter <py_event_base, event_base, event_base_exporter>::export_class ("py_event_base");
  }

} // namespace python
} // namespace blue_sky
#endif
