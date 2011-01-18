#ifndef PY_MESH_GRDECL_H
#define PY_MESH_GRDECL_H

#ifdef BSPY_EXPORTING_PLUGIN
#include "py_rs_mesh.h"
#include "bs_mesh_grdecl.h"

namespace blue_sky {
namespace python {

  PY_EXPORTER (mesh_grdecl_exporter, rs_mesh_iface_exporter)
    .def ("get_ext_to_int", &T::get_ext_to_int, args(""), "Return reference to external-to-internal mesh index")
    .def ("get_int_to_ext", &T::get_int_to_ext, args(""), "Return reference to internal-to-external mesh index")
    .def ("get_volumes", &T::get_volumes, args(""), "Return reference to volumes vector")
    .def ("get_dimensions_range", &T::get_dimensions_range, args("dim1_max, dim1_min, dim2_max, dim2_min, dim3_max, dim3_min"), "Get dimensions ranges")
    .def ("get_element_size", &T::get_element_size, args ("n_elem, dx, dy, dz"), "get elements sizes")
    .def ("get_element_ijk_to_int", &T::get_element_ijk_to_int, args ("i, j, k"), "get elements sizes")
    .def ("get_n_active_elements", &T::get_n_active_elements, args (""), "Get elements sizes")
    .def ("calc_element_tops", &T::calc_element_tops, args (""), "Calc element tops");
  PY_EXPORTER_END;

  void 
  py_export_mesh_grdecl ();

} // namespace python
} // namespace blue_sky

#endif
#endif // PY_MESH_GRDECL_H
