/** 
 * @file py_traj.h
 * @brief python wraper for #traj
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-28
 */
#ifndef PY_TRAJ_KJZ43FXM

#define PY_TRAJ_KJZ43FXM


#include <string>
#include "traj_iface.h"

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

  PY_EXPORTER (py_traj_exporter, default_exporter)
    .def ("get_table",                       &T::get_table, 
        args (""), "Return table with depth and DEV curves")
    .def ("read_from_dev_file",               &T::read_from_dev_file, 
        args ("fname"), "Read wellbore trajectory from DEV file format")
    .def ("__str__",                            &T::py_str)
  PY_EXPORTER_END;                               

  //PY_EXPORTER (py_dummy_exporter, default_exporter)
                                                 
  //PY_EXPORTER_END;                               
                                                 
  //! export matrices to python                  
  void py_export_traj ();                    

    
  } // namespace python
} // namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN

#endif /* end of include guard: PY_TRAJ_KJZ43FXM */
