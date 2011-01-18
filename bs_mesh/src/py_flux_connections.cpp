#include "bs_mesh_stdafx.h"

#include "py_flux_connections.h"
#include "bs_flux_connections.h"
#include "export_python_wrapper.h"

#ifdef BSPY_EXPORTING_PLUGIN
using namespace boost::python;

namespace blue_sky
  {
  namespace python
    {
    template <class strategy_t>
    py_flux_connections<strategy_t>::py_flux_connections ()
    : py_objbase (BS_KERNEL.create_object (wrapped_t::bs_type ()))
    {}

    template<class strategy_t>
    typename py_flux_connections<strategy_t>::sp_bcsr_matrix_t
    py_flux_connections<strategy_t>::get_conn_trans () const
    {
      return this-> get_spx (this)->get_conn_trans();
    }
    template <class strategy_t>
    void py_export_flux_connections_t (const char *class_name)
    {
      class_< py_flux_connections<strategy_t>, bases<py_objbase> >(class_name)
        .def ("get_conn_trans",&py_flux_connections<strategy_t>::get_conn_trans)
        ;
    }

    void py_export_flux_connections ()
    {
      py_export_flux_connections_t <base_strategy_fif> ("flux_connections_fif");
      py_export_flux_connections_t <base_strategy_did> ("flux_connections_did");
      py_export_flux_connections_t <base_strategy_dif> ("flux_connections_dif");
      
      strategy_exporter::export_base <flux_connections_iface, flux_conn_exporter> ("flux_conn_iface");
      strategy_exporter::export_class <bs_flux_connections, flux_connections_iface, flux_conn_exporter> ("flux_conn");
    }

    template class py_flux_connections< base_strategy_fif >;
    template class py_flux_connections< base_strategy_did >;
  }
}
#endif