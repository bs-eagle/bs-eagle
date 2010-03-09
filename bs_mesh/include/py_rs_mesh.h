#ifndef PY_RS_MESH_H
#define PY_RS_MESH_H

#ifdef BSPY_EXPORTING_PLUGIN
#include "bcsr_matrix.h"
#include "py_flux_connections.h"
#include "rs_mesh_iface.h"

#include "export_python_wrapper.h"

namespace blue_sky {
namespace python {
   
  PY_EXPORTER (rs_mesh_iface_exporter, default_exporter)
    .def ("init_props", &T::init_props)
    .def ("init_ext_to_int", &T::init_ext_to_int)
    .def ("build_jacobian_and_flux_connections", &T::build_jacobian_and_flux_connections)
  PY_EXPORTER_END;

  void 
  py_export_mesh ();

} // namespace python
} // namespace blue_sky

#endif
#endif // PY_RS_MESH_H
