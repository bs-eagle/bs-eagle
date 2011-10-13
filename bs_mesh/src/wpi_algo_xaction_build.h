/// @file wpi_algo_xaction_build.h
/// @brief Fastest, but of limited usage, algorithm for finding intersection points
/// @author uentity
/// @version 
/// @date 06.10.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_ALGO_XACTION_BUILD_8LRGTK5D
#define WPI_ALGO_XACTION_BUILD_8LRGTK5D

#include "wpi_algo_xaction.h"

namespace blue_sky { namespace wpi {

template< class strat_t >
class intersect_builder : public intersect_base< strat_t > {
public:
	typedef intersect_base< strat_t > base_t;

	// inherit types
	typedef typename base_t::trimesh      trimesh;
	typedef typename base_t::well_path    well_path;
	typedef typename base_t::mesh_part    mesh_part;
	typedef typename base_t::hit_idx_t    hit_idx_t;

	typedef typename strat_t::vertex_pos_i vertex_pos_i;
	typedef typename strat_t::Segment      Segment;

	typedef std::multimap< ulong, mesh_part >     search_space;
	typedef typename search_space::iterator       ss_iterator;
	typedef typename search_space::const_iterator css_iterator;

	typedef typename base_t::template xbbox< base_t::D > xbbox_t;

	// base functions
	using base_t::calc_hit_idx;
	using base_t::check_intersection;

	// base members
	using base_t::hit_idx_;
	using base_t::m_;
	using base_t::m_size_;
	using base_t::wp_;

	// propagate ctor
	intersect_builder(trimesh& mesh, well_path& wp, const vertex_pos_i& mesh_size)
		: base_t(mesh, wp, mesh_size)
	{}

	hit_idx_t& build() {
		typedef typename mesh_part::container_t meshp_container;
		typedef typename meshp_container::iterator meshp_iterator;

		// prepare hit index
		calc_hit_idx();

		// storage of mesh parts that really intersects with well
		search_space space;
		// create list of mesh parts for each well segment
		for(ulong i = 0; i < hit_idx_.size() - 1; ++i) {
			mesh_part seg_m(m_, m_size_);
			seg_m.init(hit_idx_[i], hit_idx_[i + 1]);
			space.insert(std::make_pair(i, seg_m));
		}

		// let's go
		while(space.size()) {
			// split each mesh part and intersect splitting with well path
			search_space div_space;
			//bool do_intersect;
			for(css_iterator pp = space.begin(), end = space.end(); pp != end; ++pp) {
				meshp_container kids = pp->second.divide();

				// test for intersections with corresponding well segment
				const ulong wseg_id = pp->first;
				//wp_iterator pw = wp_.begin() + wseg_id;
				//if(pw == wp_.end()) continue;
				const Segment& seg = wp_[wseg_id].segment();

				for(meshp_iterator pk = kids.begin(), kend = kids.end(); pk != kend; ++pk) {
					// is segment is null-length (vertical well in 2D)
					// then check if point lie on mesh rect boundary
					//if(seg.is_degenerate())
					//	do_intersect = pk->iso_bbox().has_on_boundary(seg.source());
					//// otherwise check that segment intersect with this mesh rect
					//else
					//	do_intersect = CGAL::do_intersect(seg, xbbox_t::get(*pk));

					//if(!seg.is_degenerate() && strat_t::bbox_segment_x(pk->bbox(), seg)) {
					if(!seg.is_degenerate() && CGAL::do_intersect(seg, xbbox_t::get(*pk))) {
						// mesh parts of only 1 cell goes to result
						if(pk->size() == 1)
							// find intersection points if any
							check_intersection(pk->ss_id(0), wseg_id, seg);
						else
							div_space.insert(std::make_pair(wseg_id, *pk));
					}
				}
			}

			// update search space
			space.clear();
			space.insert(div_space.begin(), div_space.end());
			//space = div_space;
		}

		return hit_idx_;
	}
};

}} /* blue_sky::wpi */

#endif /* end of include guard: WPI_ALGO_XACTION_BUILD_8LRGTK5D */

