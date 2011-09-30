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

#include "wpi_common.h"
#include "wpi_algo_pod.h"
#include "wpi_algo_meshp.h"
#include "wpi_algo_xaction.h"

#include "conf.h"
#include "bs_mesh_grdecl.h"
// DEBUG
//#include <iostream>

namespace blue_sky { namespace wpi {

/*-----------------------------------------------------------------
 * implement well path identification algos depending on strategy
 *----------------------------------------------------------------*/
template< class strat_t >
struct wpi_algo : public wpi_algo_helpers< strat_t > {
	// import strategy typedefs
	typedef typename strat_t::Point    Point;
	typedef typename strat_t::Segment  Segment;
	typedef typename strat_t::Bbox     Bbox;
	typedef typename strat_t::Iso_bbox Iso_bbox;

	typedef typename strat_t::vertex_pos   vertex_pos;
	typedef typename strat_t::vertex_pos_i vertex_pos_i;

	// import global consts
	enum { D = strat_t::D };

	// import pods
	typedef wpi_algo_pod< strat_t > wpi_pod;
	typedef typename wpi_pod::cell_data cell_data;
	typedef typename wpi_pod::sp_cell_data sp_cell_data;
	typedef typename wpi_pod::trimesh trimesh;
	typedef typename wpi_pod::trim_iterator trim_iterator;
	typedef typename wpi_pod::ctrim_iterator ctrim_iterator;

	typedef typename wpi_pod::well_data well_data;
	typedef typename wpi_pod::well_path well_path;
	typedef typename wpi_pod::wp_iterator wp_iterator;
	typedef typename wpi_pod::cwp_iterator cwp_iterator;

	typedef typename wpi_pod::well_hit_cell well_hit_cell;
	typedef typename wpi_pod::intersect_path intersect_path;

	// import mesh_part
	typedef wpi_algo_meshp< strat_t > wpi_meshp;
	typedef typename wpi_meshp::mesh_part mesh_part;

	// import intersect_action
	typedef wpi_algo_xaction< strat_t > wpi_xaction;
	typedef typename wpi_xaction::intersect_action intersect_action;

  typedef smart_ptr< bs_mesh_grdecl, true > sp_grd_mesh;

	// helper to create initial cell_data for each cell
	static spv_float coord_zcorn2trimesh(t_long nx, t_long ny, spv_float coord, spv_float zcorn,
			trimesh& res, vertex_pos_i& mesh_size)
	{
		// build mesh_grdecl around given mesh
		sp_grd_mesh grd_src = BS_KERNEL.create_object(bs_mesh_grdecl::bs_type());
		grd_src->init_props(nx, ny, coord, zcorn);
		// init mesh size
		const ulong full_sz[] = {nx, ny, (zcorn->size() / (nx * ny)) >> 3};
		const ulong n_cells = ulong(full_sz[0] * full_sz[1] * full_sz[2]);
		std::copy(full_sz, full_sz + D, mesh_size);

		// obtain coordinates for all vertices of all cells
		spv_float tops = grd_src->calc_cells_vertices_xyz();
		// clear COORD & ZCORN arrays
		grd_src->clear();

		// fill trimesh with triangles corresponding to each cell
		res.resize(n_cells);
		v_float::iterator pv = tops->begin();
		for(ulong i = 0; i < n_cells; ++i) {
			// DEBUG
			//if(i < 100) {
			//	std::cout << std::fixed << std::setprecision(2);
			//	for(uint j = 0; j < 24; ++j)
			//		std::cout << *(pv + j) << ' ';
			//	std::cout << std::endl;
			//}
			res[i] = cell_data(&*pv);
			pv += 3*8;
		}

		return tops;
	}

	/*-----------------------------------------------------------------
	* implementation of main routine
	*----------------------------------------------------------------*/
	template< bool pythonish, class = void >
	struct wpi_return {
		typedef spv_float type;

		static type make(intersect_action& A) {
			return A.export_1d();
		}
	};

	template< class unused >
	struct wpi_return< false, unused > {
		typedef std::vector< well_hit_cell > type;

		static type make(intersect_action& A) {
			type res(A.path().size());
			ulong i = 0;
			for(typename intersect_path::const_iterator px = A.path().begin(), end = A.path().end(); px != end; ++px)
				res[i++] = *px;

			return res;
		}
	};

	template< bool pythonish >
	static typename wpi_return< pythonish >::type well_path_ident_d(
		t_long nx, t_long ny, spv_float coord, spv_float zcorn,
		spv_float well_info, bool include_well_nodes)
	{
		typedef typename wpi_return< pythonish >::type ret_t;

		// 1) calculate mesh nodes coordinates and build initial trimesh
		trimesh M;
		vertex_pos_i mesh_size;
		spv_float tops = coord_zcorn2trimesh(nx, ny, coord, zcorn, M, mesh_size);
		// free memory, don't need mesh any more
		coord.release(); zcorn.release();
		// DEBUG
		//std::cout << "trimesh built" << std::endl;

		// 2) create well path description and
		// bounding boxes for line segments representing well trajectory
		ulong well_node_num = well_info->size() >> 2;
		if(well_node_num < 2) return ret_t();

		// storage
		well_path W(well_node_num - 1);
		//std::vector< Box > well_boxes(well_node_num - 1);
		// build array of well nodes as Point_2
		std::vector< Point > wnodes(well_node_num);

		// walk along well
		v_float::iterator pw = well_info->begin();
		//double md = 0;
		well_data wd;
		for(ulong i = 0; i < well_node_num - 1; ++i) {
			if(i)
				wd = well_data(pw, &W[i - 1]);
			else
				wd = well_data(pw);
	
			// insert well segment
			W[i] = wd;
			wnodes[i] = wd.start();

			pw += 4;
		}
		// put last node to array
		wnodes[well_node_num - 1] = W[well_node_num - 2].finish();
		// DEBUG
		//std::cout << "well_path created" << std::endl;

		// 3) find where each node of well is located
		// to restrict search area

		// find where well path nodes are located
		intersect_action A(M, W, mesh_size);
		const std::vector< ulong >& hit_idx = wpi_meshp::where_is_point(M, mesh_size, wnodes);
		// DEBUG
		//std::cout << "hit_idx found" << std::endl;
		// dump hit_idx
		//for(ulong i = 0; i < hit_idx.size(); ++i)
		//	std::cout << hit_idx[i] << ' ';
		//std::cout << std::endl;

		// narrow search space via branch & bound algo
		//const std::vector< ulong > hit_idx = A.build2();
		A.build(hit_idx);
		// DEBUG
		//std::cout << "build() done" << std::endl;

		// remove duplicates in X,Y,Z directions
		A.remove_dups2();
		// DEBUG
		//std::cout << "remove_dups2 done" << std::endl;

		// finalize intersection
		if(include_well_nodes)
			A.append_wp_nodes(hit_idx);
		// DEBUG
		//std::cout << "well nodes inserted" << std::endl;

		return wpi_return< pythonish >::make(A);
	}
};

}} /* blue_sky::wpi */

#endif /* end of include guard: WPI_ALGO_I21Y0RBS */
