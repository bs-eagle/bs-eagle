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

#include "wpi_common.h"
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
	// import strategy typedefs
	typedef typename strat_t::Point    Point;
	typedef typename strat_t::Segment  Segment;
	typedef typename strat_t::Bbox     Bbox;
	typedef typename strat_t::Iso_bbox Iso_bbox;

	typedef typename strat_t::vertex_pos   vertex_pos;
	typedef typename strat_t::vertex_pos_i vertex_pos_i;
	typedef typename strat_t::cell_pos     cell_pos;

	// import global consts
	enum { D = strat_t::D, CVN = strat_t::CVN, inner_point_id = strat_t::inner_point_id };

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

	// helper to resolve issue with CGAL that can intersect only Bbox_3 in 3D
	// and Iso_rectangle_2 in 2D! holy shit
	template< int dims, class = void >
	struct meshp2xbbox {
		static Bbox get(const mesh_part& mp) {
			return mp.bbox();
		}
	};
	template< class unused >
	struct meshp2xbbox< 2, unused > {
		static Iso_bbox get(const mesh_part& mp) {
			return mp.iso_bbox();
		}
	};

	/*-----------------------------------------------------------------
	 * branch & bound algorithm for finding cells that really intersect with well
	 *----------------------------------------------------------------*/
	struct branch_bound {
		typedef typename wpi_xaction::mesh_box_handle mesh_box_handle;
		typedef typename strat_t::xpoints_list xpoints_list;

		// well segment -> list of corresponding mesh parts
		typedef std::multimap< ulong, mesh_part > search_space;
		typedef typename search_space::iterator ss_iterator;
		typedef typename search_space::const_iterator css_iterator;

		typedef std::list< trim_iterator > result_t;

		branch_bound(const well_path& W, trimesh& M,
			const vertex_pos_i& mesh_size, const std::vector< ulong > hit_idx)
			: W_(W)
		{
			init(M, mesh_size, hit_idx);
		}

		void init(trimesh& M, const vertex_pos_i& mesh_size, const std::vector< ulong > hit_idx) {
			// create list of mesh parts for each well segment
			for(ulong i = 0; i < hit_idx.size() - 1; ++i) {
				mesh_part seg_m(M, mesh_size);
				seg_m.init(hit_idx[i], hit_idx[i + 1]);
				space_.insert(std::make_pair(i, seg_m));
			}
		}

		result_t& go() {
			typedef typename mesh_part::container_t meshp_container;

			res_.clear();
			while(space_.size()) {
				// split each mesh part and intersect splitting with well path
				search_space div_space;
				for(css_iterator pp = space_.begin(), end = space_.end(); pp != end; ++pp) {
					meshp_container kids = pp->second.divide();

					// test for intersections with corresponding well segment
					const ulong wseg_id = pp->first;
					const Segment& seg = W_.find(wseg_id)->second.segment();
					for(typename meshp_container::iterator pk = kids.begin(), kend = kids.end(); pk != kend; ++pk) {
						if(CGAL::do_intersect(seg, meshp2xbbox< D >::get(*pk))) {
							// mesh parts of only 1 cell goes to result
							if(pk->size() == 1)
								res_.push_back(const_cast< mesh_part& >(*pk).ss_iter(0));
							else
								div_space.insert(std::make_pair(wseg_id, *pk));
						}
					}
				}

				// update search space
				space_.clear();
				space_.insert(div_space.begin(), div_space.end());
				//space_ = div_space;
			}

			return res_;
		}

		// access result
		result_t& res() {
			return res_;
		}
		const result_t& res() const {
			return res_;
		}

		// bounding boxes ariund well segments
		const well_path& W_;
		// live mesh parts
		search_space space_;
		// resulting cells contained here
		result_t res_;
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

		// narrow search space via branch & bound algo
		branch_bound bb(W, M, mesh_size, hit_idx);
		typedef typename branch_bound::result_t search_space;
		search_space& s_space = bb.go();

		// create part of mesh to process based on these cells
		//mesh_part hot_mesh(M, mesh_size);
		//hot_mesh.init(hit_idx);

		// create bounding box for each cell in given mesh
		std::vector< Box > mesh_boxes(s_space.size());
		ulong cnt = 0;
		//trim_iterator pm;
		for(typename search_space::iterator ps = s_space.begin(), end = s_space.end(); ps != end; ++ps) {
			const cell_data& d = (*ps)->second;
			//pm = hot_mesh.ss_iter(i);
			//const cell_data& d = pm->second;
			mesh_boxes[cnt++] = Box(d.bbox(), new cell_box_handle(*ps));
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

