#include "bs_mesh_stdafx.h"

#include "py_mesh_grdecl.h"
#include "export_python_wrapper.h"
#include "py_pair_converter.h"
#include <boost/python/tuple.hpp>

#ifdef BSPY_EXPORTING_PLUGIN
using namespace boost::python;

namespace {
using namespace blue_sky;
using namespace blue_sky::python;

// same as grdecl exporter, but also export gen_coord_zcorn
template <typename T>
struct mesh_grdecl_exporter_plus {
	typedef t_long int_t;
	typedef t_double fp_t;
	typedef spv_float spfp_storarr_t;
	typedef spv_long spi_arr_t;
	typedef typename spi_arr_t::pure_pointed_t int_arr_t;

	static tuple refine_mesh(int_t nx, int_t ny, spfp_storarr_t coord, spfp_storarr_t zcorn, spfp_storarr_t points,
			fp_t m_thresh = DEF_CELL_MERGE_THRESHOLD, fp_t b_thresh = DEF_BAND_THRESHOLD)
	{
		spi_arr_t hit_idx = BS_KERNEL.create_object(int_arr_t::bs_type());
		std::pair< spfp_storarr_t, spfp_storarr_t > r = T::refine_mesh(nx, ny, coord, zcorn, points, hit_idx, m_thresh, b_thresh);
		return make_tuple(r.first, r.second, nx, ny, hit_idx);
	}

	// overloads
	static tuple refine_mesh1(int_t nx, int_t ny, spfp_storarr_t coord, spfp_storarr_t zcorn, spfp_storarr_t points) {
		return refine_mesh(nx, ny, coord, zcorn, points);
	}
	static tuple refine_mesh2(int_t nx, int_t ny, spfp_storarr_t coord, spfp_storarr_t zcorn, spfp_storarr_t points, fp_t m_thresh) {
		return refine_mesh(nx, ny, coord, zcorn, points, m_thresh);
	}
	

	template <typename class_t>
	static class_t &
	export_class (class_t &class__) {
		using namespace boost::python;

		mesh_grdecl_exporter<T>::export_class (class__)
			.def("gen_coord_zcorn", &T::gen_coord_zcorn, args("nx, ny, nz, dx, dy, dz"), "Generate COORD & ZCORN from given dimensions")
			.staticmethod("gen_coord_zcorn")
			.def("refine_mesh", &refine_mesh, "Refine existing mesh in given points")
			.def("refine_mesh", &refine_mesh1, "Refine existing mesh in given points")
			.def("refine_mesh", &refine_mesh2, "Refine existing mesh in given points")
			//.def("refine_mesh", &T::refine_mesh, args("nx, ny, coord, zcorn, points"), "Refine existing mesh in given points")
			.staticmethod("refine_mesh")
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
	class_exporter< bs_mesh_grdecl, rs_mesh_iface, mesh_grdecl_exporter_plus >::export_class ("mesh_grdecl");
	//reg_sparray_pair< float >();
	reg_sparray_pair< double >();
}

} // namespace python
} // namespace blue_sky
#endif
