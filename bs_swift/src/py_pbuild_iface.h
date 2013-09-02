#ifndef PY_PBUILD_IFACE_H_
#define PY_PBUILD_IFACE_H_
#include "amg_pbuild_iface.h"

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

  PY_EXPORTER (py_pbuild_exporter, default_exporter)
    .def ("build",              
          &T::build, 
          args ("A_matrix", "n_coarse_size", "max_connections", "cf_markers", "s_markers", "P_matrix"), "Build prolangation matrix")
    .def ("__str__",                            
          &T::py_str)
                                                 
  PY_EXPORTER_END;                               
 

  //! export matrices to python                  
  void py_export_pbuild ();                    

    
  } // namespace python
} // namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
#endif // #ifndef PY_PBUILD_BASE_H_
