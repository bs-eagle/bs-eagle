/**
* \file   py_amg_solver_iface.cpp
* \brief  Python wrapper for amg solver
* \author Sayfullin Ilshat
* \date 2011-04-01
*/

#include "py_amg_solver_iface.h"
#include "amg_solver.h"

using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //! export to python
  void py_export_amg_solver ()
  {
    using namespace boost::python;

    base_exporter<amg_solver_iface, py_amg_solver_iface_exporter>::export_class ("amg_solver_iface");
    class_exporter<amg_solver, amg_solver_iface, py_amg_solver_iface_exporter>::export_class ("amg_solver");

  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN

