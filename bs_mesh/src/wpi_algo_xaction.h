/// @file wpi_algo_xaction.h
/// @brief Intersect action class actually track intersections of well & cell boxes
/// @author uentity
/// @version 
/// @date 19.09.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_ALGO_XACTION_PLJRVQ8B
#define WPI_ALGO_XACTION_PLJRVQ8B

//#include <CGAL/box_intersection_d.h>
#include <CGAL/intersections.h>

#include "wpi_algo_pod.h"
#include "wpi_algo_meshp.h"

// DEBUG
//#include <iostream>

#define MD_TOL 0.000001

namespace blue_sky { namespace wpi {

template< class strat_t >
struct wpi_algo_xaction : public wpi_algo_helpers< strat_t > {
	// basic types
	typedef typename strat_t::vertex_pos   vertex_pos;
	typedef typename strat_t::vertex_pos_i vertex_pos_i;

	typedef typename strat_t::Point    Point;
	typedef typename strat_t::Segment  Segment;
	typedef typename strat_t::Bbox     Bbox;
	typedef typename strat_t::Iso_bbox Iso_bbox;

	// import global consts
	enum { D = strat_t::D, inner_point_id = strat_t::inner_point_id };

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

	// mesh_part
	typedef wpi_algo_meshp< strat_t > wpi_meshp;
	typedef typename wpi_meshp::mesh_part mesh_part;


	/*-----------------------------------------------------------------
	* holds all data and search actual intersection points
	*----------------------------------------------------------------*/
	class intersect_action {
	public:
		typedef typename intersect_path::iterator x_iterator;
		// well segment -> list of corresponding mesh parts
		typedef std::multimap< ulong, mesh_part > search_space;
		typedef typename search_space::iterator ss_iterator;
		typedef typename search_space::const_iterator css_iterator;

		//template< int N >
		struct spatial_sort {
			typedef int dirvec_t[D];
			typedef typename intersect_path::iterator x_iterator;

			spatial_sort(const dirvec_t& dir, const vertex_pos_i& mesh_size)
				: dir_(dir), m_size_(mesh_size)
			{}

			spatial_sort(const spatial_sort& rhs)
				: dir_(rhs.dir_), m_size_(rhs.m_size_)
			{}

			bool operator()(const x_iterator& r1, const x_iterator& r2) const {
				// check if r1 or r2 are node points
				if(r1->is_node) {
					if(!r2->is_node)
						return true;
					else
						// merge node points in same position
						return false;
				}
				else if(r2->is_node)
					return false;

				// cell ids
				vertex_pos_i c1, c2;
				decode_cell_id(r1->cell, c1, m_size_);
				decode_cell_id(r2->cell, c2, m_size_);

				//bool res = false;
				for(uint i = 0; i < D; ++i) {
					int f = greater(i, c1[i], c2[i]);
					if(f > 0)
						return true;
					else if(f == 0)
						return false;
				}
				return false;
			}

			int greater(uint ndim, ulong v1, ulong v2) const {
				if(v1 == v2) return -1;
				return dir_[ndim] == 0 ? v1 > v2 : v2 > v1;
			}

			const dirvec_t& dir_;
			// TODO: remove mesh ref from here if well_hit_cell store directly cell id
			//const trimesh& m_;
			const vertex_pos_i& m_size_;
		};


		// helper to resolve issue with CGAL that can intersect only Bbox_3 in 3D
		// and Iso_rectangle_2 in 2D! holy shit
		template< int dims, class = void >
		struct meshp2xbbox {
			typedef Bbox type;
			static Bbox get(const mesh_part& mp) {
				return mp.bbox();
			}
		};
		template< class unused >
		struct meshp2xbbox< 2, unused > {
			typedef Iso_bbox type;
			static Iso_bbox get(const mesh_part& mp) {
				return mp.iso_bbox();
			}
		};

		// ctor
		intersect_action(trimesh& mesh, well_path& wp, intersect_path& X, const vertex_pos_i& mesh_size)
			: m_(mesh), wp_(wp), x_(X)
		{
			ca_assign(m_size_, mesh_size);
		}

		static double distance(const Point& p1, const Point& p2) {
			return std::sqrt(Segment(p1, p2).squared_length());
		}

		static double calc_md(const wp_iterator& fish, const Point& target) {
			// walk all segments before the last one;
			double md = fish->md();
			// append tail
			md += distance(fish->start(), target);
			return md;
		}

		// branch & bound algorithm for finding cells that really intersect with well
		void build(const std::vector< ulong >& hit_idx) {
			typedef typename mesh_part::container_t meshp_container;
			typedef typename meshp_container::iterator meshp_iterator;

			// storage of mesh parts that really intersects with well
			search_space space;
			// create list of mesh parts for each well segment
			for(ulong i = 0; i < hit_idx.size() - 1; ++i) {
				mesh_part seg_m(m_, m_size_);
				seg_m.init(hit_idx[i], hit_idx[i + 1]);
				space.insert(std::make_pair(i, seg_m));
			}

			// let's go
			while(space.size()) {
				// split each mesh part and intersect splitting with well path
				search_space div_space;
				bool do_intersect;
				for(css_iterator pp = space.begin(), end = space.end(); pp != end; ++pp) {
					meshp_container kids = pp->second.divide();

					// test for intersections with corresponding well segment
					const ulong wseg_id = pp->first;
					wp_iterator pw = wp_.begin() + wseg_id;
					if(pw == wp_.end()) continue;
					const Segment& seg = pw->second.segment();

					for(meshp_iterator pk = kids.begin(), kend = kids.end(); pk != kend; ++pk) {
						// is segment is null-length (vertical well in 2D)
						// then check if point lie on mesh rect boundary
						if(seg.is_degenerate())
							do_intersect = pk->iso_bbox().has_on_boundary(seg.source());
						// otherwise check that segment intersect with this mesh rect
						else
							do_intersect = CGAL::do_intersect(seg, meshp2xbbox< D >::get(*pk));

						if(do_intersect) {
							// mesh parts of only 1 cell goes to result
							if(pk->size() == 1)
								// find intersection points if any
								check_intersection(pk->ss_id(0), pw, seg);
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
		}

		std::vector< ulong > build2() {
			typedef typename meshp2xbbox< D >::type xrect_t;
			// mesh partition stored here
			typedef typename mesh_part::container_t parts_container;
			typedef typename mesh_part::container_t::iterator part_iterator;

			// list of mesh parts
			parts_container parts;
			parts.insert(mesh_part(m_, m_size_));

			// cache list of well nodes
			std::vector< Segment > wseg(wp_.size());
			//wseg.reserve(wp_.size());
			for(ulong i = 0; i < wp_.size(); ++i) {
				wseg[i] = wp_[i].segment();
				//const Segment& s = wp_[i].segment();
				//if(!s.is_degenerate())
				//	wseg.push_back(s);
			}

			// result - hit index
			std::vector< ulong > res(wp_.size() + 1, m_.size());

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
					const xrect_t& cur_rect = meshp2xbbox< D >::get(*l);
					//const Iso_bbox& cur_ibb = l->iso_bbox();
					const Bbox& cur_bbox = l->bbox();

					std::list< ulong > catched_points;
					for(ulong i = 0; i <= wp_.size(); ++i) {
						// skip already found points
						if(res[i] < m_.size()) continue;
						// check that point lies inside this part
						if(wpi_meshp::point_inside_bbox(cur_bbox, ss_wp(i))) {
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
						ulong cell_id = encode_cell_id(l->lo, m_size_);
						cell_data& cell = m_[cell_id];
						for(std::list< ulong >::iterator pp = catched_points.begin(),
							cp_end = catched_points.end();
							pp != cp_end; ++pp
							)
						{
							if(cell.contains(ss_wp(*pp)))
								res[*pp] = cell_id;
						}

						// check for intersection points
						for(std::list< ulong >::iterator ps = catched_seg.begin(),
							cs_end = catched_seg.end(); ps != cs_end; ++ps
							)
							check_intersection(l->ss_id(0), wp_.begin() + *ps, wseg[*ps]);

						leafs.erase(l++);
					}
					else
						++l;
				}

				// leafs become the new start point for further division
				parts = leafs;
			}

			return res;
		}

		// run it after all dups killed
		void append_wp_nodes(const std::vector< ulong >& hit_idx) {
			if(!wp_.size()) return;

			// walk through the intersection path and add node points
			// of well geometry to the cell with previous intersection
			x_iterator px = x_.begin();
			//ulong facet_id;

			wp_iterator pw = wp_.begin();
			ulong node_idx = 0;
			for(wp_iterator end = wp_.end(); pw != end; ++pw, ++node_idx) {
				//const well_data& wseg = pw->second;
				// upper_bound
				while(px != x_.end() && px->md < pw->md())
					++px;
				// we need prev intersection
				//if(px != x_.begin()) --px;

				px = insert_wp_node(hit_idx[node_idx], pw, px);
			}

			// well path doesn't contain the end-point of trajectory
			// add it manually
			insert_wp_node(hit_idx[node_idx], --pw, x_.end(), true);
		}

		//template< int N >
		void remove_dups2() {
			typedef int dirvec_t[D];
			//typedef intersect_path::iterator x_iterator;

			struct top_surv {
				typedef std::set< x_iterator, spatial_sort > spat_storage_t;
				typedef typename spat_storage_t::iterator spat_iterator;

				top_surv(const dirvec_t& dir, intersect_path& x, const vertex_pos_i& mesh_size)
					: dir_(dir), x_(x), m_size_(mesh_size)
				{}

				x_iterator operator()(x_iterator from, x_iterator to) {
					spat_storage_t sr(spatial_sort(dir_, m_size_));

					// spatially sort iterators
					for(x_iterator px = from; px != to; ++px)
						sr.insert(px);
					// save only frst element
					while(from != to) {
						if(from != *sr.begin())
							x_.erase(from++);
						else ++from;
					}

					return *sr.begin();
				}

				const dirvec_t& dir_;
				intersect_path& x_;
				const vertex_pos_i& m_size_;
			};

			// main processing cycle
			// sanity check
			if(x_.size() < 2) return;

			// position on first intersection
			x_iterator px = x_.begin();
			// walk the nodes and determine direction of trajectory
			double max_md;
			dirvec_t dir;
			for(wp_iterator pw = wp_.begin(), end = wp_.end(); pw != end; ++pw) {
				// identify direction
				const well_data& seg = *pw;
				// calc direction vector
				Point start = seg.start();
				Point finish = seg.finish();
				for(uint i = 0; i < 2; ++i)
					dir[i] = start[i] <= finish[i] ? 0 : 1;
				// judge
				top_surv judge(dir, x_, m_size_);

				// remove dups lying on current well segment
				max_md = seg.md() + seg.len();
				for(; px != x_.end() && px->md <= max_md; ++px) {
					// skeep well node points if any
					//if(px->facet == 4)
					//	continue;

					// find range of cross points with equal MD
					// upper_bound
					ulong cnt = 0;
					x_iterator pn = px;
					for(; pn != x_.end() && abs(px->md - pn->md) < MD_TOL; ++pn, ++cnt)
					{}
					// if we have nonempty range - leave only 1 element
					if(cnt)
						px = judge(px, pn);
				}
			}
		}

		spv_float export_1d() const {
			spv_float res = BS_KERNEL.create_object(v_float::bs_type());
			res->resize(x_.size() * (4 + D));
			vf_iterator pr = res->begin();

			for(typename intersect_path::const_iterator px = x_.begin(), end = x_.end(); px != end; ++px) {
				// cell num
				*pr++ = px->cell;
				// MD
				*pr++ = px->md;
				// intersection point
				for(uint i = 0; i < D; ++i)
					*pr++ = px->where[i];
				// facet id
				*pr++ = px->facet;
				// node flag
				*pr++ = px->is_node;
			}

			return res;
		}

	private:
		x_iterator insert_wp_node(ulong cell_id, wp_iterator pw, x_iterator px, bool end_point = false) {
			// initialization
			const well_data& wseg = *pw;
			Point where = wseg.start();
			t_float wp_md = wseg.md();
			if(end_point) {
				where = wseg.finish();
				wp_md += wseg.len();
				//wp_md = px->md + distance(px->where, where);
			}

			// check if current or prev intersection match with node
			// TODO: refactor this code, looks ugly
			uint facet_id = inner_point_id;
			bool prev_is_node = false;
			if(px != x_.end() && std::abs(px->md - wp_md) < MD_TOL) {
				facet_id = px->facet;
				prev_is_node = px->is_node;
			}
			else if(px != x_.begin()) {
				--px;
				if(std::abs(px->md - wp_md) < MD_TOL) {
					facet_id = px->facet;
					prev_is_node = px->is_node;
				}
			}

			// if node point coinside with existing intersection
			// and that intersection is not a node
			// then just set is_node flag (dont't affect ordering)
			// otherwise insert new intersection point
			if(facet_id == inner_point_id || prev_is_node)
				px = x_.insert(well_hit_cell(
					where, pw, cell_id, wp_md,
					facet_id, true
				));
			else {
				//std::cout << facet_id << ' ' << inner_point_id << " facet_id != inner_point_id" << std::endl;
				well_hit_cell& x = const_cast< well_hit_cell& >(*px);
				x.is_node = true;
			}
			return px;
		}

		void check_intersection(ulong cell_id, wp_iterator pw, const Segment& well_seg) {
			typedef typename strat_t::xpoints_list xpoints_list;

			// find intersection points coord if any
			const xpoints_list& res = strat_t::precise_intersection(m_[cell_id], well_seg);

			// add all points to interscetion path
			for(typename xpoints_list::const_iterator px = res.begin(), end = res.end(); px != end; ++px) {
				x_.insert(
					well_hit_cell(px->first, pw, cell_id, calc_md(pw, px->first), px->second)
				);
			}
		}

		// subscript for accessing well path nodes
		inline Point ss_wp(ulong point_idx) const {
			if(point_idx < wp_.size())
				return wp_[point_idx].start();
			else
				return wp_[wp_.size() - 1].finish();
		}

		// mesh
		trimesh& m_;
		// well path
		well_path& wp_;
		// well-mesh intersection path
		intersect_path& x_;
		// mesh size
		vertex_pos_i m_size_;
	};

};

}} /* blue_sky::wpi */

#endif /* end of include guard: WPI_ALGO_XACTION_PLJRVQ8B */

