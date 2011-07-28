/** 
 * @file py_traj.cpp
 * @brief python wraper for #traj
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-28
 */

#include "py_traj.h"
#include "traj.h"

using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //! export matrices to python
  void py_export_traj ()
  {
    using namespace boost::python;

    base_exporter <traj_iface, py_traj_exporter>::export_class ("traj_iface");

    class_exporter <traj, traj_iface, py_traj_exporter>::export_class ("traj");

  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
