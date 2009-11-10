#include "bs_mesh_stdafx.h"

#include "rs_mesh_iface.h"
#include "flux_connections_iface.h"
#include "py_rs_mesh.h"
#include "py_flux_connections.h"

using namespace boost::python;

namespace blue_sky {
namespace python {

  void py_export_mesh ()
  {
    strategy_exporter::export_base <rs_mesh_iface, rs_mesh_iface_exporter> ("rs_mesh");
  }

} // namespace python
} // namespace blue_sky
