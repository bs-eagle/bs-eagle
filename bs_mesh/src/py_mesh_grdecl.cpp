#include "bs_mesh_stdafx.h"

#include "py_mesh_grdecl.h"
#include "export_python_wrapper.h"
#include "py_pair_converter.h"

#ifdef BSPY_EXPORTING_PLUGIN
using namespace boost::python;

namespace {
using namespace blue_sky;
using namespace blue_sky::python;

// same as grdecl exporter, but also export gen_coord_zcorn
template <typename T>
struct mesh_grdecl_exporter_plus {
	template <typename class_t>
	static class_t &
	export_class (class_t &class__) {
		using namespace boost::python;
		typedef typename T::fp_type_t fp_t;
		typedef bs_array< fp_t > fp_array_t;
		typedef smart_ptr< fp_array_t > spfp_array_t;
		//typedef typename T::fp_array_t fp_array_t;
		//typedef typename T::spfp_array_t spfp_array_t;

		mesh_grdecl_exporter<T>::export_class (class__)
			.def("gen_coord_zcorn", &T::gen_coord_zcorn, args("nx, ny, nz, dx, dy, dz"), "Generate COORD & ZCORN from given dimensions")
			.staticmethod("gen_coord_zcorn")
			//.def ("gen_coord_zcorn", &T::template gen_coord_zcorn< fp_t, fp_t, fp_t>, args("nx, ny, nz, dx, dy, dz"), "Generate COORD & ZCORN from given dimensions")
			;
		return class__;
	}
};

template< class fp_type >
void reg_sparray_pair() {
	// register converter from pair of returned arrayes to python list
	typedef smart_ptr< bs_array< fp_type > > spfp_array;
	typedef std::pair< spfp_array, spfp_array > array_pair;
	typedef bspy_converter< pair_traits< array_pair > > array_pair_converter;
	array_pair_converter::register_to_py();
	array_pair_converter::register_from_py();
}

} // eof hidden namespace

namespace blue_sky {
namespace python {

void py_export_mesh_grdecl () {
	strategy_exporter::export_class <bs_mesh_grdecl, rs_mesh_iface, mesh_grdecl_exporter_plus> ("mesh_grdecl");
	reg_sparray_pair< float >();
	reg_sparray_pair< double >();
}

} // namespace python
} // namespace blue_sky
#endif
