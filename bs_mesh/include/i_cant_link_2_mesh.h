/// @file i_cant_link_2_mesh.h
/// @brief Stupid interface just to prevent linking bs_mesh
/// @author uentity
/// @version 
/// @date 14.10.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef I_CANT_LINK_2_MESH_KXZEOZXG
#define I_CANT_LINK_2_MESH_KXZEOZXG

#include "rs_smesh_iface.h"

namespace blue_sky {

class handy_mesh_iface : public objbase {
public:
	typedef smart_ptr< rs_smesh_iface, true > sp_smesh;

	virtual spv_float calc_cells_vertices_xyz(t_long nx, t_long ny, spv_float coord, spv_float zcorn) = 0;
	//virtual spv_ulong where_is_points(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_float points) = 0;
	//virtual spv_ulong where_is_points_2d(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_float points) = 0;
	//virtual t_ulong where_is_point(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_float point) = 0;
	//virtual t_ulong where_is_point_2d(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_float point) = 0;
	virtual sp_smesh make_mesh_grdecl(t_long nx, t_long ny, spv_float coord, spv_float zcorn) = 0;
};
typedef smart_ptr< handy_mesh_iface > sp_himesh;

} /* blue_sky */

#endif /* end of include guard: I_CANT_LINK_2_MESH_KXZEOZXG */

