#include "bs_mesh_stdafx.h"
#include "py_flux_connections.h"

#ifdef BSPY_EXPORTING_PLUGIN
using namespace boost::python;

namespace blue_sky
  {
  namespace python
    {
    py_flux_connections::py_flux_connections ()
    : py_objbase (BS_KERNEL.create_object (wrapped_t::bs_type ()))
    {}

    py_flux_connections::sp_bcsr_matrix_t
    py_flux_connections::get_conn_trans () const
    {
      return this-> get_spx (this)->get_conn_trans();
    }
    void py_export_flux_connections_t (const char *class_name)
    {
      class_< py_flux_connections, bases<py_objbase> >(class_name)
        .def ("get_conn_trans",&py_flux_connections::get_conn_trans)
        ;
    }

    void py_export_flux_connections ()
    {
      py_export_flux_connections_t ("flux_connections");
    }
  }
}
#endif
