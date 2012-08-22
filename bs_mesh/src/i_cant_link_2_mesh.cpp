/// @file i_cant_link_2_mesh.cpp
/// @brief Implementation of handy interface to mesh
/// @author uentity
/// @version 
/// @date 14.10.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "bs_mesh_stdafx.h"
#include "i_cant_link_2_mesh.h"
#include "bs_mesh_grdecl.h"
#include "well_path_ident.h"

namespace blue_sky {

namespace {

// separate function for Python export
spv_float calc_cells_vertices_xyz_impl(t_long nx, t_long ny, spv_float coord, spv_float zcorn) {
	typedef smart_ptr< bs_mesh_grdecl, true > sp_grd_mesh;
	// build mesh_grdecl around given mesh
	sp_grd_mesh grd_src = BS_KERNEL.create_object(bs_mesh_grdecl::bs_type());
	if(!grd_src) return NULL;
	grd_src->init_props(nx, ny, coord, zcorn);
	// obtain coordinates for all vertices of all cells
	return grd_src->calc_cells_vertices_xyz();
}

class BS_API_PLUGIN handy_object : public handy_mesh_iface {
public:
	spv_float calc_cells_vertices_xyz(t_long nx, t_long ny, spv_float coord, spv_float zcorn) {
		return calc_cells_vertices_xyz_impl(nx, ny, coord, zcorn);
	}

	// 3D
	spv_uint where_is_points(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_float points) {
		return blue_sky::where_is_points(nx, ny, coord, zcorn, points);
	}
	// 2D
	spv_uint where_is_points_2d(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_float points) {
		return blue_sky::where_is_points_2d(nx, ny, coord, zcorn, points);
	}

	// 3D
	t_uint where_is_point(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_float point) {
		return blue_sky::where_is_point(nx, ny, coord, zcorn, point);
	}
	// 2D
	t_uint where_is_point_2d(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_float point) {
		return blue_sky::where_is_point_2d(nx, ny, coord, zcorn, point);
	}

	BLUE_SKY_TYPE_DECL(handy_object)
};

// std and copy ctors
handy_object::handy_object(bs_type_ctor_param)
{}

handy_object::handy_object(const handy_object& rhs)
	: bs_refcounter()
{}

}  // eof hidden namespace

BLUE_SKY_TYPE_IMPL(handy_object, objbase, "handy_mesh_iface", "Interface to everything needed from mesh", "")
BLUE_SKY_TYPE_STD_CREATE(handy_object)
BLUE_SKY_TYPE_STD_COPY(handy_object)

bool register_handy_mesh_iface(const plugin_descriptor& pd) {
	return BS_KERNEL.register_type(pd, handy_object::bs_type());
}

#ifdef BSPY_EXPORTING_PLUGIN
namespace python {
	using namespace boost::python;

	void py_export_handymesh() {
		def("calc_cells_vertices_xyz", &calc_cells_vertices_xyz_impl);
	}
} /* python */
#endif

} /* blue_sky */
