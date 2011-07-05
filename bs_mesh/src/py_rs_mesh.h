#ifndef PY_RS_MESH_H
#define PY_RS_MESH_H

#ifdef BSPY_EXPORTING_PLUGIN
//#include "bcsr_matrix.h"
#include "rs_mesh_iface.h"

#include "export_python_wrapper.h"

namespace blue_sky {
namespace python {
   
  template< class T >
  struct rs_mesh_iface_exporter {
    template <typename class_t>
    static class_t &
    export_class (class_t &class__) {
      using namespace boost::python;

      // resolve ambiguity in resolving init_props
      void (T::*init_props_hdm)(const typename T::sp_hdm_t) = &T::init_props;
      default_exporter< T >::export_class(class__)
        .def ("init_props", init_props_hdm)
        .def ("init_ext_to_int", &T::init_ext_to_int)
        .def ("build_jacobian_and_flux_connections", &T::build_jacobian_and_flux_connections,
            args("jacobian", "flux_connections", "boundary_array"), "Build jacobian structure, calculate transmissibilities")
        ;

      return class__;
    }
  };

  //PY_EXPORTER (rs_mesh_iface_exporter, default_exporter)
  //  .def ("init_props", &T::init_props, args("hdm"), "Initialize mesh properties")
  //  .def ("init_ext_to_int", &T::init_ext_to_int, args(""), "Initialize internal indexation")
  //  //.def ("build_jacobian_and_flux_connections", &T::build_jacobian_and_flux_connections, 
  //  //      args("jacobian", "flux_connections", "boundary_array"), "Build jacobian structure, calculate transmissibilities")
  //PY_EXPORTER_END;

  void 
  py_export_mesh ();

} // namespace python
} // namespace blue_sky

#endif
#endif // PY_RS_MESH_H
