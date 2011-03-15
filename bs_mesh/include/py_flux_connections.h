#ifndef PY_FLUX_CONNECTIONS_H
#define PY_FLUX_CONNECTIONS_H

#ifdef BSPY_EXPORTING_PLUGIN
#include "flux_connections_iface.h"


namespace blue_sky {
namespace python {

    
    class py_flux_connections : public py_objbase
    {
    public:
      typedef flux_connections_iface                  wrapped_t;              //!< wrapped type flux_connections

      typedef smart_ptr<wrapped_t, true>                          sp_flux_connections_t;  //!< smart_ptr to flux_connections
      typedef smart_ptr <bcsr_matrix_iface, true>     sp_bcsr_matrix_t;
      
      py_flux_connections ();
      virtual ~py_flux_connections () {} 
      
      sp_bcsr_matrix_t
      get_conn_trans () const;

    };
    
    PY_EXPORTER (flux_conn_iface_exporter, default_exporter)
    PY_EXPORTER_END;
    
    PY_EXPORTER (flux_conn_exporter, flux_conn_iface_exporter)
    PY_EXPORTER_END;

    void 
    py_export_flux_connections ();

} // namespace python
} // namespace boost

#endif
#endif // PY_FLUX_CONNECTIONS_H
