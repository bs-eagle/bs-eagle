/** 
 * @file py_well.h
 * @brief python wraper for well wellture storage
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-29
 */
#ifndef PY_WELL_O4VRW03K
#define PY_WELL_O4VRW03K

#include "well_iface.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "dummy_base.h"
#include "construct_python_object.h"
#include "make_me_happy.h"
#include "export_python_wrapper.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include <string>

#ifdef BSPY_EXPORTING_PLUGIN
namespace blue_sky
  {
  namespace python
    {

  PY_EXPORTER (py_well_exporter, default_exporter)
    .def ("get_prop",                         &T::get_prop, 
        args (""), "Return well properties")
    .def ("get_branch_names",                 &T::get_branch_names, 
        args (""), "Return list of branches")
    .def ("get_branch",                       &T::get_branch, 
        args ("Name"), "Return branch by name")
    .def ("__str__",                          &T::py_str)
  PY_EXPORTER_END;                               

  //! export matrices to python                  
  void py_export_well ();                    

    
  } // namespace python
} // namespace blue_sky

#endif // #ifdef BSPY_EXPORTING_PLUGIN
#endif /* end of include guard: PY_well_O4VRW03K */

