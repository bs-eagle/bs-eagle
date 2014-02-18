/// @file wpi_algo_xaction_build_mp.h
/// @brief Build intersection points for multiple well_paths in one turn
/// @author uentity
/// @version 1.0
/// @date 18.02.2014
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_ALGO_XACTION_BUILD_MP_C9I39ZHJ
#define WPI_ALGO_XACTION_BUILD_MP_C9I39ZHJ

#include "wpi_algo_xaction_build2.h"

namespace blue_sky { namespace wpi {

template< class strat_t >
class intersect_builder_mp : public pods< strat_t > {
public:
	typedef pods< strat_t > base_t;
	typedef intersect_builder2< strat_t > xbuild1;

	// inherit types
	typedef typename base_t::cell_data    cell_data;
	typedef typename base_t::well_data    well_data;
	typedef typename base_t::trimesh      trimesh;
	typedef typename base_t::well_path    well_path;
	typedef typename base_t::well_paths   well_paths;
	typedef typename xbuild1::mesh_part   mesh_part;

	typedef typename strat_t::vertex_pos_i vertex_pos_i;
	typedef typename strat_t::Point        Point;
	typedef typename strat_t::Segment      Segment;

	typedef typename xbuild1::mesh_tools_t mesh_tools_t;
	typedef typename mesh_part::container_t parts_container;
	typedef typename parts_container::iterator part_iterator;
	typedef typename parts_container::const_iterator part_citerator;

	enum { D = strat_t::D };

	// intersection builder for single well_path
	// import some usefull stuff from it
	typedef typename xbuild1::template xbbox< D > xbbox_t;
	typedef typename xbuild1::cell_box_handle cell_box_handle;
	typedef typename xbuild1::well_box_handle well_box_handle;
	typedef typename xbuild1::mesh_box_handle mesh_box_handle;
	typedef typename xbuild1::Box Box;
	typedef typename xbuild1::leafs_builder leafs_builder;

	typedef typename xbuild1::intersect_path intersect_path;
	typedef std::vector< intersect_path >    intersect_paths;

	typedef typename xbuild1::hit_idx_t       hit_idx_t;
	typedef typename std::vector< hit_idx_t > hit_idxs_t;

	typedef intersect_base< strat_t > xbuild_base;

	// comparator to sort Boxes by handle values
	struct Box_compare {
		bool operator()(const Box& lhs, const Box& rhs) {
			return *lhs.handle() < *rhs.handle();
		}
	};
	typedef std::set< Box, Box_compare > Boxes;
	typedef typename Boxes::iterator Boxes_iterator;
	typedef typename Boxes::const_iterator Boxes_citerator;

	// ctor
	intersect_builder_mp(trimesh& mesh, well_paths& wp)
		: m_(mesh), wps_(wp)
	{}

	// branch & bound algorithm for finding cells that really intersect with well
	// works using boxes intersection
	hit_idxs_t& build(bool append_wp_nodes = true) {
		typedef std::vector< Segment > Segments;
		typedef std::vector< Segments > Segments_mp;

		Segments_mp wseg(wps_.size());
		std::vector< std::vector< Box > > well_boxes(wps_.size());
		std::vector< std::vector< Box* > > well_boxes_p(wps_.size());

		// make separate base intersector for each well_path
		xbricks_.clear();
		xbricks_.reserve(wps_.size());

		// list of mesh parts
		// new parts is temp storage for intersected parts with given leafs_builder
		parts_container parts, new_parts;
		parts.insert(mesh_part(m_));

		// pools for allocating well & mesh box handles
		// they all will be detroyed automatically on exit
		boost::object_pool< well_box_handle > wbh_pool;
		boost::object_pool< mesh_box_handle > mbh_pool;

		// cache well segments
		// and make boxes for them
		for(ulong i = 0; i < wps_.size(); ++i) {
			for(ulong j = 0; j < wps_[i].size(); ++i) {
				wseg[i][j] = wps_[i][j].segment();
				// NOTE: use value of j as local well_box identificator
				well_boxes[i][j] = Box(wps_[i][j].bbox(), new(wbh_pool.malloc()) well_box_handle(j));
				//well_boxes[i] = Box(wp_[i].bbox(), new well_box_handle(i));
				well_boxes_p[i][j] = &well_boxes[i][j];
			}

			xbricks_[i] = xbuild_base(m_, wps_[i]);
			// init hit index with invalid value
			std::fill(xbricks_[i].hit_idx().begin(), xbricks_[i].hit_idx().end(), m_.size_flat());
		}

		// result - hit index
		//hit_idx_.resize(wp_.size() + 1);
		//std::fill(hit_idx_.begin(), hit_idx_.end(), m_.size_flat());

		// let's go
		while(parts.size()) {
			// we need container to hold all mesh parts
			// list can be used because all mesh_parts in parts
			// are different at this point, so results of division will also
			// be unqie
			std::list< mesh_part > leafs;

			// make boxes around parts
			// try to use single boxes container over all iterations - no profit
			// store mesh boxes instances in std::set (to remove marked for split boxes)
			//const ulong max_boxes = parts.size() * (1 << D);
			Boxes mp_boxes;
			//mp_boxes.reserve(max_boxes);

			// split every part
			//ulong i = 0;
			for(part_iterator p = parts.begin(), end = parts.end(); p != end; ++p) {
				const parts_container kids = p->divide();
				//ulong i = 0;
				for(part_citerator pk = kids.begin(), kend = kids.end(); pk != kend; ++pk) {
					mp_boxes.insert(
						Box(pk->bbox(),
							new(mbh_pool.malloc()) mesh_box_handle(
								&*leafs.insert(leafs.end(), *pk)
							)
						)
					);
					//mp_boxes_p.push_back(&mp_boxes[mp_boxes.size() - 1]);
				}
			}

			// IMPORTANT: clear parts container
			// boxes to split on next step are stored here by leafs_builder
			parts.clear();

			// do intersect boxes for each well_path
			for(ulong i = 0; i < wps_.size(); ++i) {
				// make mesh box pointers
				std::vector< Box* > mp_boxes_p(mp_boxes.size());
				ulong j = 0;
				for(
					Boxes_citerator p_box = mp_boxes.begin(), e = mp_boxes.end();
					p_box != e; ++p_box, ++j
				) {
					mp_boxes_p[j] = &*p_box;
				}

				// IMPORTANT: clear new_parts before processing next leafs_builder
				new_parts.clear();
				// intersect boxes
				CGAL::box_intersection_d(
					mp_boxes_p.begin(),   mp_boxes_p.end(),
					well_boxes_p[i].begin(), well_boxes_p[i].end(),
					leafs_builder(xbricks_[i], wseg[i], xbricks_[i].hit_idx(), new_parts)
				);

				// parts marked for further split can be excluded from mp_boxes
				for(part_citerator pk = new_parts.begin(), e = new_parts.end(); pk != e; ++pk) {
					// erase kid from boxes storage
					// use specifix Box ctor, provide only handle for searching
					const Boxes_citerator pb = mp_boxes.find(Box(false, &*pk));
					if(pb != mp_boxes.end())
						mp_boxes.erase(pb);

					// and add part for further splitting
					parts.insert(*pk);
				}
			}
		}

		// collect all hit_idx from xbricks
		hit_idxs_t res(xbricks_.size());
		for(ulong i = 0; i < xbricks_.size(); ++i) {
			res[i] = xbricks_[i].hit_idx();
			if(append_wp_nodes)
				xbricks_[i].append_wp_nodes(res[i]);
		}
		return res;
	}

	// export resulting intersections as vector of 1D arrays
	std::vector< spv_float > export_1d() const {
		std::vector< spv_float > res;
		res.reserve(xbricks_.size());
		for(ulong i = 0; i < xbricks_.size(); ++i)
			res.push_back(xbricks_[i].export_1d());
		return res;
	}

	const std::vector< xbuild_base >& path() const {
		return xbricks_;
	}

	// mesh
	trimesh& m_;
	// well paths
	well_paths& wps_;
	// base intersection builders that store resulting intersections
	// for each well
	std::vector< xbuild_base > xbricks_;
	// cell IDs of where each well node is located
	//hit_idxs_t hit_idx_;
};

}} /* blue_sky::wpi */

#endif /* end of include guard: WPI_ALGO_XACTION_BUILD_MP_C9I39ZHJ */

