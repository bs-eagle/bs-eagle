/// @file wpi_algo.h
/// @brief Algorithms implementing well path identification using given strategy
/// @author uentity
/// @version 
/// @date 15.09.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_ALGO_I21Y0RBS
#define WPI_ALGO_I21Y0RBS

#include "conf.h"
#include "bs_mesh_grdecl.h"

#include "wpi_algo_pod.h"
#include "wpi_algo_meshp.h"
#include "wpi_algo_xaction.h"

#include <iterator>
#include <cmath>

namespace blue_sky { namespace wpi {

/*-----------------------------------------------------------------
 * implement well path identification algos depending on strategy
 *----------------------------------------------------------------*/
template< class strat_t >
struct wpi_algo : public wpi_algo_helpers< strat_t > {
	// common typedefs
	typedef t_ulong ulong;
	typedef t_uint uint;
	typedef v_float::iterator vf_iterator;
	typedef smart_ptr< bs_mesh_grdecl > sp_grd_mesh;

	// import strategy typedefs
	typedef typename strat_t::Object   Object;
	typedef typename strat_t::Kernel   Kernel;
	typedef typename strat_t::Point    Point;
	typedef typename strat_t::Segment  Segment;
	typedef typename strat_t::Bbox     Bbox;
	typedef typename strat_t::Iso_bbox Iso_bbox;

	typedef typename strat_t::vertex_pos   vertex_pos;
	typedef typename strat_t::vertex_pos_i vertex_pos_i;
	typedef typename strat_t::cell_pos     cell_pos;

	// import global consts
	enum { D = strat_t::D, CVN = strat_t::CVN, inner_point_id = strat_t::inner_point_id };
	//using strat_t::D;
	//using strat_t::CVN;
	//using strat_t::inner_point_id;

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
	typedef typename wpi_xaction::box_handle box_handle;
	typedef typename wpi_xaction::cell_box_handle cell_box_handle;
	typedef typename wpi_xaction::well_box_handle well_box_handle;
	typedef typename wpi_xaction::Box Box;
	typedef typename wpi_xaction::intersect_action intersect_action;

	/*-----------------------------------------------------------------
	 * branch & bound algorithm for finding intersections
	 *----------------------------------------------------------------*/
	struct branch_bound_intersect {
		
	};

	// helper to create initial cell_data for each cell
	static spv_float coord_zcorn2trimesh(t_long nx, t_long ny, spv_float coord, spv_float zcorn, trimesh& res) {
		// build mesh_grdecl around given mesh
		sp_grd_mesh grd_src = BS_KERNEL.create_object(bs_mesh_grdecl::bs_type());
		grd_src->init_props(nx, ny, coord, zcorn);
		t_long nz = (zcorn->size() / nx / ny) >> 3;
		ulong n_cells = ulong(nx * ny * nz);

		// obtain coordinates for all vertices of all cells
		spv_float tops = grd_src->calc_cells_vertices_xyz();
		v_float::iterator pv = tops->begin();

		// fill trimesh with triangles corresponding to each cell
		for(ulong i = 0; i < n_cells; ++i) {
			res[i] = cell_data(&*pv);
			pv += 3*8;
		}

		return tops;
	}

	/*-----------------------------------------------------------------
	* implementation of main routine
	*----------------------------------------------------------------*/
	static spv_float well_path_ident_d(t_long nx, t_long ny, spv_float coord, spv_float zcorn,
		spv_float well_info, bool include_well_nodes)
	{
		// 1) calculate mesh nodes coordinates and build initial trimesh
		trimesh M;
		spv_float tops;
		tops = coord_zcorn2trimesh(nx, ny, coord, zcorn, M);
		t_long nz = (zcorn->size() / nx / ny) >> 3;

		// 2) create well path description and
		// bounding boxes for line segments representing well trajectory
		ulong well_node_num = well_info->size() >> 2;
		if(well_node_num < 2) return spv_float();

		// storage
		well_path W;
		std::vector< Box > well_boxes(well_node_num - 1);
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
			wnodes[i] = wd.start();

			// make bbox
			well_boxes[i] = Box(
				wd.bbox(),
				new well_box_handle(W.insert(std::make_pair(i, wd)).first)
			);

			pw += 4;
		}
		// put last node to array
		wnodes[well_node_num - 1] = W[well_node_num - 2].finish();

		// 3) find where each node of well is located
		// to restrict search area
		// intersections storage
		intersect_path X;
		// TODO: init mesh_size in better way
		ulong full_mesh_size[] = {nx, ny, nz};
		vertex_pos_i mesh_size;
		std::copy(full_mesh_size, full_mesh_size + D, mesh_size);

		// find where well path nodes are located
		intersect_action A(M, W, X, mesh_size);
		const std::vector< ulong >& hit_idx = wpi_meshp::where_is_point(M, mesh_size, wnodes);

		// create list of mesh parts for each well segment
		//std::list< mesh_part > search_space;
		//for(ulong i = 0; i < hit_idx.size() - 1; ++i) {
		//	mesh_part seg_m(M, mesh_size);
		//	seg_m.init(hit_idx[i], hit_idx[i + 1]);
		//	search_space.push_back(seg_m);
		//}

		// split current mesh parts and find intersections on each split


		// create part of mesh to process based on these cells
		mesh_part hot_mesh(M, mesh_size);
		hot_mesh.init(hit_idx);

		// create bounding box for each cell in given mesh
		std::vector< Box > mesh_boxes(hot_mesh.size());
		ulong cnt = 0;
		trim_iterator pm;
		for(ulong i = 0; i < hot_mesh.size(); ++i) {
			pm = hot_mesh.ss_iter(i);
			const cell_data& d = pm->second;
			mesh_boxes[cnt++] = Box(d.bbox(), new cell_box_handle(pm));
		}


		// Run the intersection algorithm with all defaults on the
		// indirect pointers to cell bounding boxes. Avoids copying the boxes
		CGAL::box_intersection_d(
			mesh_boxes.begin(), mesh_boxes.end(),
			well_boxes.begin(), well_boxes.end(),
			A
		);

		// remove duplicates in X,Y,Z directions
		A.remove_dups2();

		// finalize intersection
		if(include_well_nodes)
			A.append_wp_nodes(hit_idx);

		return A.export_1d();
	}
};

}} /* blue_sky::wpi */

#endif /* end of include guard: WPI_ALGO_I21Y0RBS */

