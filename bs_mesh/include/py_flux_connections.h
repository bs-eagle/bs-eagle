#ifndef PY_FLUX_CONNECTIONS_H
#define PY_FLUX_CONNECTIONS_H

#include "flux_connections_iface.h"

namespace blue_sky {
namespace python {

    template <typename strategy_t>
    class py_flux_connections : public py_objbase
    {
    public:
      typedef typename strategy_t::item_array_t                   item_array_t;
      typedef typename strategy_t::index_array_t                  index_array_t;
      typedef typename strategy_t::csr_matrix_t                   bcsr_matrix_t;
      typedef flux_connections_iface<strategy_t>                  wrapped_t;              //!< wrapped type flux_connections

      typedef smart_ptr<wrapped_t, true>                          sp_flux_connections_t;  //!< smart_ptr to flux_connections
      typedef smart_ptr <bcsr_matrix_t, true>                     sp_bcsr_matrix_t;
      
      py_flux_connections ();
      virtual ~py_flux_connections () {} 
      
      sp_bcsr_matrix_t
      get_conn_trans () const;

    };

    void 
    py_export_flux_connections ();

} // namespace python
} // namespace boost

#endif // PY_FLUX_CONNECTIONS_H
