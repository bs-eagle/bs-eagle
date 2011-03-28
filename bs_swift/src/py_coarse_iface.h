#ifndef PY_COARSE_IFACE_H_
#define PY_COARSE_IFACE_H_
#include "amg_coarse_iface.h"

//#include "bdiag_matrix.h"
//#include "jac_matrix.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "dummy_base.h"
#include "construct_python_object.h"
#include "make_me_happy.h"
#include "python/py_bs_object_base.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include "export_python_wrapper.h"

#ifdef BSPY_EXPORTING_PLUGIN
namespace blue_sky
  {
  namespace python
    {

  PY_EXPORTER (py_coarse_exporter, default_exporter)
    .def ("build",              
          &T::build, 
          args ("matrix", "wksp", "cf_markers", "s_markers"), "Create coarse grid")
    .def ("__str__",                            
          &T::py_str)
                                                 
  PY_EXPORTER_END;                               
 

  //! export matrices to python                  
  void py_export_coarse ();                    

    
  } // namespace python
} // namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
#endif // #ifndef PY_COARSE_BASE_H_
