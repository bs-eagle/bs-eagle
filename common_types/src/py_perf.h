/** 
 * @file py_perf.h
 * @brief python wraper for well perf storage
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-29
 */
#ifndef PY_PERF_O4VRW03K

#define PY_PERF_O4VRW03K


#include <string>
#include "perf_iface.h"

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

  PY_EXPORTER (py_perf_exporter, default_exporter)
    .def ("get_prop",                         &T::get_prop, 
        args (""), "Return perf properties")
    .def ("__str__",                            &T::py_str)
  PY_EXPORTER_END;                               

  //! export matrices to python                  
  void py_export_perf ();                    

    
  } // namespace python
} // namespace blue_sky

#endif // #ifdef BSPY_EXPORTING_PLUGIN
#endif /* end of include guard: PY_FRAC_O4VRW03K */

