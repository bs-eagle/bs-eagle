#include "py_rs_mesh.h"
#include "bs_mesh_grdecl.h"

#include "export_python_wrapper.h"
#include "py_pair_converter.h"
#include <boost/python/tuple.hpp>

#ifdef BSPY_EXPORTING_PLUGIN
using namespace boost::python;

namespace blue_sky { namespace python {

void init_props1(bs_mesh_grdecl& m, t_long nx, t_long ny, t_long nz,
	spv_float dx, spv_float dy, spv_float dz)
{
	m.init_props(nx, ny, nz, dx, dy, dz);
}

void init_props2(bs_mesh_grdecl& m, t_long nx, t_long ny, spv_float coord, spv_float zcorn) {
	m.init_props(nx, ny, coord, zcorn);
}

PY_EXPORTER (mesh_grdecl_exporter, rs_mesh_iface_exporter)
	.def ("get_ext_to_int", &T::get_ext_to_int, args(""), "Return reference to external-to-internal mesh index")
	.def ("get_int_to_ext", &T::get_int_to_ext, args(""), "Return reference to internal-to-external mesh index")
	.def ("get_volumes", &T::get_volumes, args(""), "Return reference to volumes vector")
	.def ("get_cell_volumes", &T::get_cell_volumes, args(""), "Return reference to volumes vector")
	.def ("get_dimensions_range", &T::get_dimensions_range, args("dim1_max, dim1_min, dim2_max, dim2_min, dim3_max, dim3_min"), "Get dimensions ranges")
	.def ("get_element_size", &T::get_element_size, args ("n_elem, dx, dy, dz"), "get elements sizes")
	.def ("get_element_dz", &T::get_element_dim3_size, args ("n_elem"), "get elements dz by internal index")
	.def ("get_element_dz_ext", &T::get_element_dim3_size_ext, args ("i, j, k"), "get elements dz by external index")
	.def ("get_element_ijk_to_int", &T::get_element_ijk_to_int, args ("i, j, k"), "get elements sizes")
	.def ("get_n_active_elements", &T::get_n_active_elements, args (""), "Get elements sizes")
	.def ("calc_element_tops", &T::calc_element_tops, args (""), "Calc element tops")
	.def ("calc_element_center", &T::calc_element_center, args (""), "Calc element center")
	.def ("calc_cells_vertices", &T::calc_cells_vertices)
	.def ("calc_cells_vertices_xyz", &T::calc_cells_vertices_xyz)
	.def ("init_props", &init_props1)
	.def ("init_props", &init_props2)
  .def ("check_data", &T::check_data)
	.def ("get_element_sizes", &T::get_element_sizes)
  .def ("intersect_trajectories", &T::intersect_trajectories)
	;
PY_EXPORTER_END;

}}  // eof blue_sky::python

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
	typedef std::pair< spv_float, spv_float > coord_zcorn_pair;

	// gen_coord_zcorn overloads
	static coord_zcorn_pair gen_coord_zcorn1(int_t nx, int_t ny, int_t nz, spfp_storarr_t dx, spfp_storarr_t dy, spfp_storarr_t dz,
			fp_t x0, fp_t y0)
	{
		return T::gen_coord_zcorn(nx, ny, nz, dx, dy, dz, x0, y0);
	}
	static coord_zcorn_pair gen_coord_zcorn2(int_t nx, int_t ny, int_t nz, spfp_storarr_t dx, spfp_storarr_t dy, spfp_storarr_t dz,
			fp_t x0)
	{
		return T::gen_coord_zcorn(nx, ny, nz, dx, dy, dz, x0);
	}
	static coord_zcorn_pair gen_coord_zcorn3(int_t nx, int_t ny, int_t nz, spfp_storarr_t dx, spfp_storarr_t dy, spfp_storarr_t dz)
	{
		return T::gen_coord_zcorn(nx, ny, nz, dx, dy, dz);
	}

	// refine_mesh_deltas with overloads
	static tuple refine_mesh_deltas(int_t nx, int_t ny, spfp_storarr_t coord, spfp_storarr_t points,
			fp_t m_thresh = DEF_CELL_MERGE_THRESHOLD, fp_t b_thresh = DEF_BAND_THRESHOLD)
	{
		spi_arr_t hit_idx = BS_KERNEL.create_object(int_arr_t::bs_type());
		std::pair< spfp_storarr_t, spfp_storarr_t > r = T::refine_mesh_deltas(nx, ny, coord, points, hit_idx, m_thresh, b_thresh);
		return make_tuple(r.first, r.second, nx, ny, hit_idx);
	}

	static tuple refine_mesh_deltas1(int_t nx, int_t ny, spfp_storarr_t coord, spfp_storarr_t points) {
		return refine_mesh_deltas(nx, ny, coord, points);
	}
	static tuple refine_mesh_deltas2(int_t nx, int_t ny, spfp_storarr_t coord, spfp_storarr_t points, fp_t m_thresh) {
		return refine_mesh_deltas(nx, ny, coord, points, m_thresh);
	}

	// refine_mesh with overloads
	static tuple refine_mesh(int_t nx, int_t ny, spfp_storarr_t coord, spfp_storarr_t zcorn, spfp_storarr_t points,
			fp_t m_thresh = DEF_CELL_MERGE_THRESHOLD, fp_t b_thresh = DEF_BAND_THRESHOLD)
	{
		spi_arr_t hit_idx = BS_KERNEL.create_object(int_arr_t::bs_type());
		std::pair< spfp_storarr_t, spfp_storarr_t > r = T::refine_mesh(nx, ny, coord, zcorn, points, hit_idx, m_thresh, b_thresh);
		return make_tuple(r.first, r.second, nx, ny, hit_idx);
	}

	static tuple refine_mesh1(int_t nx, int_t ny, spfp_storarr_t coord, spfp_storarr_t zcorn, spfp_storarr_t points) {
		return refine_mesh(nx, ny, coord, zcorn, points);
	}
	static tuple refine_mesh2(int_t nx, int_t ny, spfp_storarr_t coord, spfp_storarr_t zcorn, spfp_storarr_t points, fp_t m_thresh) {
		return refine_mesh(nx, ny, coord, zcorn, points, m_thresh);
	}

	// refine_mesh_deltas with overloads
	// (i,j) points format
	static tuple refine_mesh_deltas_ij(int_t nx, int_t ny, spfp_storarr_t coord,
			spi_arr_t points_pos, spfp_storarr_t points_param,
			fp_t m_thresh = DEF_CELL_MERGE_THRESHOLD, fp_t b_thresh = DEF_BAND_THRESHOLD)
	{
		spi_arr_t hit_idx = BS_KERNEL.create_object(int_arr_t::bs_type());
		std::pair< spfp_storarr_t, spfp_storarr_t > r = T::refine_mesh_deltas(nx, ny, coord, points_pos, points_param,
				hit_idx, m_thresh, b_thresh);
		return make_tuple(r.first, r.second, nx, ny, hit_idx);
	}

	static tuple refine_mesh_deltas_ij1(int_t nx, int_t ny, spfp_storarr_t coord,
			spi_arr_t points_pos, spfp_storarr_t points_param)
	{
		return refine_mesh_deltas_ij(nx, ny, coord, points_pos, points_param);
	}

	static tuple refine_mesh_deltas_ij2(int_t nx, int_t ny, spfp_storarr_t coord,
			spi_arr_t points_pos, spfp_storarr_t points_param, fp_t m_thresh)
	{
		return refine_mesh_deltas_ij(nx, ny, coord, points_pos, points_param, m_thresh);
	}

	// refine_mesh with overloads
	// (i, j) points format
	static tuple refine_mesh_ij(int_t nx, int_t ny, spfp_storarr_t coord, spfp_storarr_t zcorn,
			spi_arr_t points_pos, spfp_storarr_t points_param,
			fp_t m_thresh = DEF_CELL_MERGE_THRESHOLD, fp_t b_thresh = DEF_BAND_THRESHOLD)
	{
		spi_arr_t hit_idx = BS_KERNEL.create_object(int_arr_t::bs_type());
		std::pair< spfp_storarr_t, spfp_storarr_t > r = T::refine_mesh(nx, ny, coord, zcorn, points_pos, points_param,
				hit_idx, m_thresh, b_thresh);
		return make_tuple(r.first, r.second, nx, ny, hit_idx);
	}

	static tuple refine_mesh_ij1(int_t nx, int_t ny, spfp_storarr_t coord, spfp_storarr_t zcorn,
			spi_arr_t points_pos, spfp_storarr_t points_param)
	{
		return refine_mesh_ij(nx, ny, coord, zcorn, points_pos, points_param);
	}

	static tuple refine_mesh_ij2(int_t nx, int_t ny, spfp_storarr_t coord, spfp_storarr_t zcorn,
			spi_arr_t points_pos, spfp_storarr_t points_param, fp_t m_thresh)
	{
		return refine_mesh_ij(nx, ny, coord, zcorn, points_pos, points_param, m_thresh);
	}

	template <typename class_t>
	static class_t &
	export_class (class_t &class__) {
		using namespace boost::python;

		mesh_grdecl_exporter<T>::export_class (class__)
			.def("gen_coord_zcorn", &T::gen_coord_zcorn, "Generate COORD & ZCORN from given dimensions")
			.def("gen_coord_zcorn", &gen_coord_zcorn1, "Generate COORD & ZCORN from given dimensions")
			.def("gen_coord_zcorn", &gen_coord_zcorn2, "Generate COORD & ZCORN from given dimensions")
			.def("gen_coord_zcorn", &gen_coord_zcorn3, "Generate COORD & ZCORN from given dimensions")
      .def("gen_coord2", &T::gen_coord2, "Generate COORD")
      .staticmethod("gen_coord2")
			.staticmethod("gen_coord_zcorn")
			.def("refine_mesh_deltas", &refine_mesh_deltas, "Calc dx and dy arrays for refined mesh in given points")
			.def("refine_mesh_deltas", &refine_mesh_deltas1, "Calc dx and dy arrays for refined mesh in given points")
			.def("refine_mesh_deltas", &refine_mesh_deltas2, "Calc dx and dy arrays for refined mesh in given points")
			.def("refine_mesh_deltas", &refine_mesh_deltas_ij, "Calc dx and dy arrays for refined mesh in given points")
			.def("refine_mesh_deltas", &refine_mesh_deltas_ij1, "Calc dx and dy arrays for refined mesh in given points")
			.def("refine_mesh_deltas", &refine_mesh_deltas_ij2, "Calc dx and dy arrays for refined mesh in given points")
			.staticmethod("refine_mesh_deltas")
			.def("refine_mesh", &refine_mesh, "Refine existing mesh in given points")
			.def("refine_mesh", &refine_mesh1, "Refine existing mesh in given points")
			.def("refine_mesh", &refine_mesh2, "Refine existing mesh in given points")
			.def("refine_mesh", &refine_mesh_ij, "Refine existing mesh in given points")
			.def("refine_mesh", &refine_mesh_ij1, "Refine existing mesh in given points")
			.def("refine_mesh", &refine_mesh_ij2, "Refine existing mesh in given points")
			.staticmethod("refine_mesh")
			;
		return class__;
	}
};

template< class fp_type >
void reg_sparray_pair() {
	// register converter from pair of returned arrayes to python list
	typedef std::pair< spv_float, spv_float > array_pair;
	typedef bspy_converter< pair_traits< array_pair > > array_pair_converter;
	array_pair_converter::register_to_py();
	array_pair_converter::register_from_py();
}

} // eof hidden namespace

namespace blue_sky { namespace python {

void py_export_mesh_grdecl () {
	class_exporter< bs_mesh_grdecl, rs_mesh_iface, mesh_grdecl_exporter_plus >::export_class ("mesh_grdecl");
	//reg_sparray_pair< float >();
	reg_sparray_pair< double >();
}

}} // namespace blue_sky::python
#endif
