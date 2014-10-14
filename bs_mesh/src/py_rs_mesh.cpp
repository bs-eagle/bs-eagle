#include "rs_mesh_iface.h"
#include "flux_connections_iface.h"
#include "py_rs_mesh.h"

#ifdef BSPY_EXPORTING_PLUGIN
using namespace boost::python;

namespace blue_sky {
namespace python {

  void py_export_mesh ()
  {
    base_exporter< rs_mesh_iface, rs_mesh_iface_exporter >::export_class("rs_mesh");
  }

} // namespace python
} // namespace blue_sky
#endif
