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
#include "wpi_algo_xaction_build_mp.h"
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
	typedef typename pods_t::well_paths well_paths;
	typedef typename pods_t::wp_iterator wp_iterator;
	typedef typename pods_t::cwp_iterator cwp_iterator;

	typedef typename pods_t::well_hit_cell well_hit_cell;
	typedef typename pods_t::intersect_path intersect_path;
	typedef typename pods_t::intersect_paths intersect_paths;

	// import mesh_part
	typedef mesh_tools< strat_t > mesh_tools_t;
	typedef typename mesh_tools_t::mesh_part mesh_part;

	// import intersect_action
	typedef intersect_base< strat_t > xbase;
	typedef typename xbase::hit_idx_t hit_idx_t;
	typedef intersect_builder2< strat_t > xbuilder;
	// multiple wells
	typedef intersect_builder_mp< strat_t > xbuilder_mp;
	typedef typename xbuilder_mp::hit_idxs_t hit_idxs_t;

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
			// insert well segment
			W[i] = well_data(pw, &W[i - 1]);
		}

		return well_node_num;
	}

	/*-----------------------------------------------------------------
	 * path identification result constructor
	 *----------------------------------------------------------------*/
	template< bool pythonish, class = void >
	struct wpi_return {
		typedef spv_float type;

		template< class intersector_t >
		static type make(const intersector_t& A) {
			return A.export_1d();
		}
	};

	template< class unused >
	struct wpi_return< false, unused > {
		typedef std::vector< well_hit_cell > type;

		template< class intersector_t >
		static type make(const intersector_t& A) {
			type res(A.path().size());
			std::copy(A.path().begin(), A.path().end(), res.begin());

			return res;
		}
	};

	// the same for multiple wells
	template< bool pythonish, class = void >
	struct wpi_return_mp {
		typedef std::vector< spv_float > type;

		template< class intersector_t >
		static type make(const intersector_t& A) {
			return A.export_1d();
		}
	};

	template< class unused >
	struct wpi_return_mp< false, unused > {
		typedef wpi_return< false > wpi_ret;
		typedef std::vector< std::vector< well_hit_cell > > type;

		template< class intersector_t >
		static type make(const intersector_t& A) {
			type res(A.xbricks().size());
			for(ulong i = 0; i < A.xbricks().size(); ++i)
				res[i] = wpi_ret::make(A.xbricks()[i]);

			return res;
		}
	};

	/*-----------------------------------------------------------------
	* implementation of main routine
	*----------------------------------------------------------------*/
	template< bool pythonish >
	static typename wpi_return< pythonish >::type well_path_ident_d(
		ulong nx, ulong ny, sp_obj trim_backend,
		spv_float well_info, bool include_well_nodes, spv_ulong H = NULL)
	{
		typedef typename wpi_return< pythonish >::type ret_t;

		// 1) calculate mesh nodes coordinates and build initial trimesh
		trimesh M(nx, ny, trim_backend);

		// 2) create well path description
		well_path W;
		if(!fill_well_path(well_info, W)) return ret_t();

		// 3) construct main object
		xbuilder A(M, W);

		// 4) narrow search space via branch & bound algo
		const hit_idx_t& hit_idx = A.build();

		// 5) remove duplicates in X,Y,Z directions
		//A.remove_dups2();

		// 6) finalize intersection
		if(include_well_nodes)
			A.append_wp_nodes(hit_idx);

		// return well nodes hit indexes
		if(H) {
			H->resize(hit_idx.size());
			if(H->size() == hit_idx.size())
				std::copy(hit_idx.begin(), hit_idx.end(), H->begin());
		}
		// return intersections
		return wpi_return< pythonish >::make(A);
	}

	template< bool pythonish >
	static typename wpi_return< pythonish >::type well_path_ident_d(
		ulong nx, ulong ny, spv_float coord, spv_float zcorn,
		spv_float well_info, bool include_well_nodes, spv_ulong H = NULL)
	{
		// calculate mesh nodes coordinates and build initial trimesh
		return well_path_ident_d< pythonish >(
			nx, ny, trimesh::create_backend(nx, ny, coord, zcorn),
			well_info, include_well_nodes, H
		);
	}

	/*-----------------------------------------------------------------
	 * well path identification for multiple wells
	 *----------------------------------------------------------------*/
	template< bool pythonish >
	static typename wpi_return_mp< pythonish >::type well_paths_ident_d(
		ulong nx, ulong ny, sp_obj trim_backend,
		std::vector< spv_float > well_info, bool include_well_nodes,
		std::vector< spv_ulong > H = std::vector< spv_ulong >()
	) {
		typedef typename wpi_return< pythonish >::type ret_t;

		// 1) calculate mesh nodes coordinates and build initial trimesh
		trimesh M(nx, ny, trim_backend);

		// 2) create well path description
		well_paths W(well_info.size());
		for(ulong i = 0; i < well_info.size(); ++i) {
			if(!fill_well_path(well_info[i], W[i]))
				return ret_t();
		}

		// 3) construct main object
		xbuilder_mp A(M, W);

		// 4) narrow search space via branch & bound algo
		const hit_idxs_t& hit_idxs = A.build(include_well_nodes);

		// return well nodes hit indexes
		for(ulong i = 0; i < std::min(H.size(), hit_idxs.size()); ++i) {
			H[i]->resize(hit_idxs[i].size());
			if(H[i]->size() == hit_idxs[i].size())
				std::copy(hit_idxs[i].begin(), hit_idxs[i].end(), H[i]->begin());
		}
		// return intersections
		return wpi_return< pythonish >::make(A);
	}

	template< bool pythonish >
	static typename wpi_return_mp< pythonish >::type well_paths_ident_d(
		ulong nx, ulong ny, spv_float coord, spv_float zcorn,
		spv_float well_info, bool include_well_nodes, spv_ulong H = NULL)
	{
		// calculate mesh nodes coordinates and build initial trimesh
		return well_paths_ident_d< pythonish >(
			nx, ny, trimesh::create_backend(nx, ny, coord, zcorn),
			well_info, include_well_nodes, H
		);
	}

}; // algo

}} /* blue_sky::wpi */

#endif /* end of include guard: WPI_ALGO_I21Y0RBS */

