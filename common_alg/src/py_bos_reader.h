/**
 * @file py_bos_reader.h
 * @brief python wraper for #bos_reader
 * @author Oleg Borschuk
 * @version
 * @date 2012-03-01
 */
#ifndef PY_BOS_READER_O4VRW03K

#define PY_BOS_READER_O4VRW03K


#include <string>
#include "bos_reader_iface.h"

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

  PY_EXPORTER (py_bos_reader_exporter, default_exporter)
    .def ("read_line",                          &T::read_line_str)
    .def ("get_next_keyword",                   &T::get_next_keyword)
    .def ("init",                               &T::init)
    .def ("__str__",                            &T::py_str)
  PY_EXPORTER_END;

  //! export matrices to python
  void py_export_bos_reader ();


  } // namespace python
} // namespace blue_sky

#endif // #ifdef BSPY_EXPORTING_PLUGIN
#endif /* end of include guard: PY_BOS_READER_O4VRW03K */

