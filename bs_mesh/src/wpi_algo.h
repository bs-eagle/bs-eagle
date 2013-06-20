/// @file wpi_algo.h
/// @brief Algorithms implementing well path identification using given strategy
/// @author uentity
/// @version 
/// @date 15.09.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_ALGO_I21Y0RBS
#define WPI_ALGO_I21Y0RBS

#include <iterator>
#include <cmath>
#include <algorithm>

#include "wpi_common.h"
#include "wpi_algo_pod.h"
#include "wpi_trimesh_impl.h"
#include "wpi_algo_meshp.h"
#include "wpi_algo_xaction.h"
#include "wpi_algo_xaction_build.h"
#include "wpi_algo_xaction_build2.h"
#include "wpi_algo_xaction_build3.h"
//#include "wpi_algo_vtk.h"

#include "conf.h"
#include "i_cant_link_2_mesh.h"

// DEBUG
//#include <iostream>

namespace blue_sky { namespace wpi {

/*-----------------------------------------------------------------
 * implement well path identification algos depending on strategy
 *----------------------------------------------------------------*/
template< class strat_t >
struct algo : public helpers< strat_t > {
	// import strategy typedefs
	typedef typename strat_t::Point    Point;
	typedef typename strat_t::Segment  Segment;
	typedef typename strat_t::Bbox     Bbox;
	typedef typename strat_t::Iso_bbox Iso_bbox;

	typedef typename strat_t::vertex_pos   vertex_pos;
	typedef typename strat_t::vertex_pos_i vertex_pos_i;

	typedef typename strat_t::cell_vertex_iterator cell_vertex_iterator;
	typedef typename strat_t::well_traj_iterator   well_traj_iterator;

	// import global consts
	enum { D = strat_t::D };

	// import pods
	typedef pods< strat_t > pods_t;
	typedef typename pods_t::cell_data cell_data;
	typedef typename pods_t::sp_cell_data sp_cell_data;
	typedef typename pods_t::trimesh trimesh;
	//typedef typename pods_t::trim_iterator trim_iterator;
	//typedef typename pods_t::ctrim_iterator ctrim_iterator;

	typedef typename pods_t::well_data well_data;
	typedef typename pods_t::well_path well_path;
	typedef typename pods_t::wp_iterator wp_iterator;
	typedef typename pods_t::cwp_iterator cwp_iterator;

	typedef typename pods_t::well_hit_cell well_hit_cell;
	typedef typename pods_t::intersect_path intersect_path;

	// import mesh_part
	typedef mesh_tools< strat_t > mesh_tools_t;
	typedef typename mesh_tools_t::mesh_part mesh_part;

	// import intersect_action
	typedef intersect_base< strat_t > xbase;
	typedef typename xbase::hit_idx_t hit_idx_t;
	typedef intersect_builder2< strat_t > xbuilder;

	static ulong fill_well_path(spv_float well_info, well_path& W) {
		ulong well_node_num = well_info->size() >> 2;
		if(well_node_num < 2) return 0;

		// storage
		W.resize(well_node_num - 1);

		// walk along well
		v_float::iterator pw = well_info->begin();
		W[0] = well_data(pw);
		//well_data wd;
		for(ulong i = 1; i < well_node_num - 1; ++i) {
			pw += 4;
			//if(i)
			//	wd = well_data(pw, &W[i - 1]);
			//else
			//	wd = well_data(pw);
	
			// insert well segment
			W[i] = well_data(pw, &W[i - 1]);
			//pw += 4;
		}

		return well_node_num;
	}

	/*-----------------------------------------------------------------
	* implementation of main routine
	*----------------------------------------------------------------*/
	template< bool pythonish, class = void >
	struct wpi_return {
		typedef spv_float type;

		static type make(xbase& A) {
			return A.export_1d();
		}
	};

	template< class unused >
	struct wpi_return< false, unused > {
		typedef std::vector< well_hit_cell > type;

		static type make(xbase& A) {
			type res(A.path().size());
			ulong i = 0;
			for(typename intersect_path::const_iterator px = A.path().begin(), end = A.path().end(); px != end; ++px)
				res[i++] = *px;

			return res;
		}
	};

	template< bool pythonish >
	static typename wpi_return< pythonish >::type well_path_ident_d(
		t_long nx, t_long ny, sp_obj trim_backend,
		spv_float well_info, bool include_well_nodes)
	{
		typedef typename wpi_return< pythonish >::type ret_t;

		// 1) calculate mesh nodes coordinates and build initial trimesh
		trimesh M(nx, ny, trim_backend);
		//vertex_pos_i mesh_size;
		//spv_float tops = coord_zcorn2trimesh(nx, ny, coord, zcorn, M, mesh_size);
		// DEBUG
		//std::cout << "trimesh built" << std::endl;

		// 2) create well path description
		well_path W;
		if(!fill_well_path(well_info, W)) return ret_t();
		// DEBUG
		//std::cout << "well_path created" << std::endl;

		// 3) construct main object
		xbuilder A(M, W);
		// DEBUG
		//std::cout << "hit_idx found" << std::endl;

		// 4) narrow search space via branch & bound algo
		hit_idx_t& hit_idx = A.build();
		// DEBUG
		//std::cout << "build() done" << std::endl;

		// 5) remove duplicates in X,Y,Z directions
		//A.remove_dups2();
		// DEBUG
		//std::cout << "remove_dups2 done" << std::endl;

		// 6) finalize intersection
		if(include_well_nodes)
			A.append_wp_nodes(hit_idx);
		// DEBUG
		//std::cout << "well nodes inserted" << std::endl;

		return wpi_return< pythonish >::make(A);
	}

	template< bool pythonish >
	static typename wpi_return< pythonish >::type well_path_ident_d(
		t_long nx, t_long ny, spv_float coord, spv_float zcorn,
		spv_float well_info, bool include_well_nodes)
	{
		// calculate mesh nodes coordinates and build initial trimesh
		return well_path_ident_d< pythonish >(
			nx, ny, trimesh::create_backend(nx, ny, coord, zcorn),
			well_info, include_well_nodes
		);
	}

}; // algo

}} /* blue_sky::wpi */

#endif /* end of include guard: WPI_ALGO_I21Y0RBS */

