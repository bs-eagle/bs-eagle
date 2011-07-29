/** 
 * @file py_frac.h
 * @brief python wraper for well fracture storage
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-29
 */
#ifndef PY_FRAC_O4VRW03K

#define PY_FRAC_O4VRW03K


#include <string>
#include "frac_iface.h"

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

  PY_EXPORTER (py_frac_exporter, default_exporter)
    .def ("get_prop",                         &T::get_prop, 
        args (""), "Return frac properties")
    .def ("__str__",                            &T::py_str)
  PY_EXPORTER_END;                               

  //! export matrices to python                  
  void py_export_frac ();                    

    
  } // namespace python
} // namespace blue_sky

#endif // #ifdef BSPY_EXPORTING_PLUGIN
#endif /* end of include guard: PY_FRAC_O4VRW03K */

