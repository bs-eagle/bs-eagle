#ifndef PY_AMG_SOLVER_IFACE_H_
#define PY_AMG_SOLVER_IFACE_H_
#include "amg_solver_iface.h"

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

  PY_EXPORTER (py_amg_solver_iface_exporter, default_exporter)
    .def ("setup",                &T::setup, args ("matrix"), "Setup AMG")
    .def ("solve",                &T::solve, args ("matrix", "rhs", "solution"), "Solve linear system")
    .def ("set_prop",             &T::set_prop, args ("sp_prop"), "Set up properties")
    .def ("get_prop",             &T::get_prop, args (""), "Return smart pointer to the properties")
    .def ("get_final_residual",   &T::get_final_residual, args (""), "Return final residual (for exact methods 0)")
    .def ("get_niters",           &T::get_niters, args (""), "Return number of iterations (for exact methods 1)")
    .def ("__str__", &T::py_str)
  PY_EXPORTER_END;

  //! export to python
  void py_export_amg_solver ();


  } // namespace python
} // namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
#endif // #ifndef PY_AMG_SOLVER_BASE_H_

