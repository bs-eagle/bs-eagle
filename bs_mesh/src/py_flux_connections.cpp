#include "bs_kernel.h"
#include "flux_connections_iface.h"
#include "bs_flux_connections.h"
#include "export_python_wrapper.h"
#include "py_bs_object_base.h"

#ifdef BSPY_EXPORTING_PLUGIN
using namespace boost::python;

namespace blue_sky { namespace python {

class py_flux_connections : public py_objbase {
public:
	typedef flux_connections_iface                  wrapped_t;              //!< wrapped type flux_connections

	typedef smart_ptr<wrapped_t, true>              sp_flux_connections_t;  //!< smart_ptr to flux_connections
	typedef smart_ptr <bcsr_matrix_iface, true>     sp_bcsr_matrix_t;

	py_flux_connections ()
		: py_objbase (BS_KERNEL.create_object (wrapped_t::bs_type ()))
	{}

	virtual ~py_flux_connections () {}

	sp_bcsr_matrix_t
		get_conn_trans () const;

};

PY_EXPORTER (flux_conn_iface_exporter, default_exporter)
PY_EXPORTER_END;

PY_EXPORTER (flux_conn_exporter, flux_conn_iface_exporter)
PY_EXPORTER_END;

py_flux_connections::sp_bcsr_matrix_t
py_flux_connections::get_conn_trans () const {
	return this-> get_spx (this)->get_conn_trans();
}

void py_export_flux_connections_t (const char *class_name) {
	class_< py_flux_connections, bases<py_objbase> >(class_name)
		.def ("get_conn_trans",&py_flux_connections::get_conn_trans)
		;
}

void py_export_flux_connections () {
	py_export_flux_connections_t ("flux_connections");

	base_exporter< flux_connections_iface, flux_conn_exporter >::export_class("flux_conn_iface");
	class_exporter< bs_flux_connections, flux_connections_iface, flux_conn_exporter >::export_class("flux_conn");
}

}}  // eof blue_sky::python
#endif
