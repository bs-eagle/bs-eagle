/// @file wpi_algo_xaction_build3.h
/// @brief Third approach to finding intersection points
/// @author uentity
/// @version 
/// @date 06.10.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_ALGO_XACTION_BUILD3_3QC754QF
#define WPI_ALGO_XACTION_BUILD3_3QC754QF

#include "wpi_algo_xaction.h"

namespace blue_sky { namespace wpi {

template< class strat_t >
class intersect_builder3 : public intersect_base< strat_t > {
public:
	typedef intersect_base< strat_t > base_t;

	// inherit types
	typedef typename base_t::cell_data    cell_data;
	typedef typename base_t::well_data    well_data;
	typedef typename base_t::trimesh      trimesh;
	typedef typename base_t::well_path    well_path;
	typedef typename base_t::mesh_part    mesh_part;
	typedef typename base_t::hit_idx_t    hit_idx_t;

	typedef typename strat_t::vertex_pos_i vertex_pos_i;
	typedef typename strat_t::Point        Point;
	typedef typename strat_t::Segment      Segment;
	typedef typename strat_t::Bbox         Bbox;

	typedef typename base_t::mesh_tools_t mesh_tools_t;

	enum { D = strat_t::D };
	typedef typename base_t::template xbbox< D > xbbox_t;

	// base functions
	using base_t::check_intersection;
	using base_t::ss_wp;

	// base members
	using base_t::hit_idx_;
	using base_t::m_;
	//using base_t::m_size_;
	using base_t::wp_;

	// propagate ctor
	intersect_builder3(trimesh& mesh, well_path& wp)
		: base_t(mesh, wp)
	{}

	// branch & bound algorithm for finding cells that really intersect with well
	// works using boxes intersection
	hit_idx_t& build() {
		typedef typename xbbox_t::type xrect_t;
		// mesh partition stored here
		typedef typename mesh_part::container_t parts_container;
		typedef typename mesh_part::container_t::iterator part_iterator;

		// list of mesh parts
		parts_container parts;
		parts.insert(mesh_part(m_));

		// cache list of normal well segments
		std::vector< Segment > wseg;
		wseg.reserve(wp_.size());
		for(ulong i = 0; i < wp_.size(); ++i) {
			//wseg[i] = wp_[i].segment();
			const Segment& s = wp_[i].segment();
			if(!s.is_degenerate())
				wseg.push_back(s);
		}

		// result - hit index
		hit_idx_.resize(wp_.size() + 1);
		std::fill(hit_idx_.begin(), hit_idx_.end(), m_.size());

		// let's got
		while(parts.size()) {
			// split every part
			parts_container leafs;
			for(part_iterator p = parts.begin(), end = parts.end(); p != end; ++p) {
				parts_container kids = p->divide();
				leafs.insert(kids.begin(), kids.end());
			}

			// collection of points inside current partition
			//std::list< ulong > catched_seg;
			//std::list< ulong > catched_points;

			// process each leaf and find points inside it
			for(part_iterator l = leafs.begin(), end = leafs.end(); l != end; ) {
				const xrect_t& cur_rect = xbbox_t::get(*l);
				//const Iso_bbox& cur_ibb = l->iso_bbox();
				const Bbox& cur_bbox = l->bbox();

				std::list< ulong > catched_points;
				for(ulong i = 0; i <= wp_.size(); ++i) {
					// skip already found points
					// check that point lies inside this part
					if(hit_idx_[i] >= m_.size() && mesh_tools_t::point_inside_bbox(cur_bbox, ss_wp(i))) {
						catched_points.push_back(i);
						// parts with > 1 cell anyway undergo splitting
						if(l->size() > 1)
							break;
					}
				}

				// check for segment intersection only if point inside
				// isn't found
				//catched_seg.clear();
				std::list< ulong > catched_seg;
				if(!catched_points.size()) {
					for(ulong i = 0; i < wp_.size(); ++i) {
						// first check if point is inside 

						// is segment is null-length (vertical well in 2D)
						// then check if point lie on mesh rect boundary
						// *update - degenerated segments checked on find points stage
						//if(wseg[i].is_degenerate()) {
						//	if(cur_ibb.has_on_boundary(wseg[i].source()))
						//		catched_seg.push_back(i);
						//}

						//const Segment& s = wp_[i].segment();
						// check that segment intersect with this mesh rect
						if(!wseg[i].is_degenerate() && CGAL::do_intersect(wseg[i], cur_rect)) {
							catched_seg.push_back(i);
							// parts with > 1 cell anyway undergo splitting
							if(l->size() > 1)
								break;
						}
					}
				}

				// if this part don't intersect any segments - remove it
				// if box contains only 1 cell - find intersection points
				if(!catched_points.size() && !catched_seg.size())
					leafs.erase(l++);
				else if(l->size() == 1) {
					// check for contained points if any
					ulong cell_id = encode_cell_id(l->lo, m_.size_());
					cell_data& cell = m_[cell_id];
					for(std::list< ulong >::iterator pp = catched_points.begin(),
						cp_end = catched_points.end();
						pp != cp_end; ++pp
						)
					{
						if(cell.contains(ss_wp(*pp)))
							hit_idx_[*pp] = cell_id;
					}

					// check for intersection points
					for(std::list< ulong >::iterator ps = catched_seg.begin(),
						cs_end = catched_seg.end(); ps != cs_end; ++ps
						)
						check_intersection(l->ss_id(0), *ps, wseg[*ps]);

					leafs.erase(l++);
				}
				else
					++l;
			}

			// leafs become the new start point for further division
			parts = leafs;
		}

		return hit_idx_;
	}
};

}} // eof blue_sky::wpi

#endif /* end of include guard: WPI_ALGO_XACTION_BUILD3_3QC754QF */

