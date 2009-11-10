#ifndef PY_MESH_GRDECL_H
#define PY_MESH_GRDECL_H

#include "py_rs_mesh.h"
#include "bs_mesh_grdecl.h"

namespace blue_sky {
namespace python {

  PY_EXPORTER (mesh_grdecl_exporter, rs_mesh_iface_exporter)
    //.def ("get_ext_to_int", &T::get_ext_to_int)
    //.def ("get_int_to_ext", &T::get_int_to_ext)
  PY_EXPORTER_END;

  void 
  py_export_mesh_grdecl ();

} // namespace python
} // namespace blue_sky

#endif // PY_MESH_GRDECL_H
