/// @file wpi_algo_xaction.h
/// @brief Intersect action class actually track intersections of well & cell boxes
/// @author uentity
/// @version 
/// @date 19.09.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_ALGO_XACTION_PLJRVQ8B
#define WPI_ALGO_XACTION_PLJRVQ8B

#include <CGAL/intersections.h>

#include "wpi_algo_pod.h"
#include "wpi_algo_meshp.h"

// DEBUG
//#include <iostream>

#define MD_TOL 1e-10

namespace blue_sky { namespace wpi {

/*-----------------------------------------------------------------
* holds all data and search actual intersection points
*----------------------------------------------------------------*/
template< class strat_t >
class intersect_base : public helpers< strat_t > {
public:
	// basic types
	typedef typename strat_t::vertex_pos   vertex_pos;
	typedef typename strat_t::vertex_pos_i vertex_pos_i;

	typedef typename strat_t::Point    Point;
	typedef typename strat_t::Segment  Segment;
	typedef typename strat_t::Bbox     Bbox;
	typedef typename strat_t::Iso_bbox Iso_bbox;

	// import global consts
	enum { D = strat_t::D, inner_point_id = strat_t::inner_point_id };
	typedef int dirvec_t[D];

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
	typedef typename intersect_path::iterator x_iterator;

	// import helper functions
	typedef helpers< strat_t > base_t;
	using base_t::decode_cell_id;

	// mesh_part
	typedef mesh_tools< strat_t > mesh_tools_t;
	typedef typename mesh_tools_t::mesh_part mesh_part;

	typedef std::vector< ulong > hit_idx_t;

	//template< int N >
	struct spatial_sort {
		typedef typename intersect_path::iterator x_iterator;

		spatial_sort(const dirvec_t& dir, const vertex_pos_i& mesh_size)
			: dir_(dir), m_size_(mesh_size)
		{}

		spatial_sort(const spatial_sort& rhs)
			: dir_(rhs.dir_), m_size_(rhs.m_size_)
		{}

		bool greater(const well_hit_cell& r1, const well_hit_cell& r2) const {
			// check if r1 or r2 are node points
			if(r1.is_node) {
				if(!r2.is_node)
					return true;
				else
					// merge node points in same position
					return false;
			}
			else if(r2.is_node)
				return false;

			// cell ids
			vertex_pos_i c1, c2;
			decode_cell_id(r1.cell, c1, m_size_);
			decode_cell_id(r2.cell, c2, m_size_);

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

		bool operator()(const x_iterator& r1, const x_iterator& r2) const {
			return greater(*r1, *r2);
		}

		int greater(uint ndim, ulong v1, ulong v2) const {
			if(v1 == v2) return -1;
			return dir_[ndim] == 0 ? v1 > v2 : v2 > v1;
		}

		const dirvec_t& dir_;
		const vertex_pos_i& m_size_;
	};


	// helper to resolve issue with CGAL that can intersect only Bbox_3 in 3D
	// and Iso_rectangle_2 in 2D! holy shit
	template< int dims, class = void >
	struct xbbox {
		typedef Bbox type;

		static type get(const mesh_part& mp) {
			return mp.bbox();
		}
		static type get(const cell_data& c) {
			return c.bbox();
		}
		static type get(const well_data& w) {
			return w.bbox();
		}
	};

	// spec for 2D
	template< class unused >
	struct xbbox< 2, unused > {
		typedef Iso_bbox type;
		static type get(const mesh_part& mp) {
			return mp.iso_bbox();
		}
		static type get(const cell_data& c) {
			return c.iso_bbox();
		}
		static type get(const well_data& w) {
			return w.iso_bbox();
		}
	};

	// ctor
	intersect_base(trimesh& mesh, well_path& wp)
		: m_(mesh), wp_(wp)
	{
		//ca_assign(m_size_, mesh_size);
	}

	static double distance(const Point& p1, const Point& p2) {
		return std::sqrt(Segment(p1, p2).squared_length());
	}

	double calc_md(ulong wseg_id, const Point& target) {
		return wp_[wseg_id].md() + distance(wp_[wseg_id].start(), target);
	}

	hit_idx_t& calc_hit_idx() {
		hit_idx_.clear();
		if(wp_.size()) {
			// prepare vector of well node points
			std::vector< Point > wnodes(wp_.size() + 1);
			ulong i = 0;
			for(; i < wp_.size(); ++i) {
				wnodes[i] = wp_[i].start();
			}
			// last point of well traj
			wnodes[i] = wp_[i - 1].finish();

			// calculate hit_idx
			hit_idx_ = mesh_tools_t::where_is_point(m_, wnodes);
		}

		return hit_idx_;
	}
	// hit_idx getters
	hit_idx_t& hit_idx() {
		return hit_idx_;
	}
	const hit_idx_t& hit_idx() const {
		return hit_idx_;
	}

	// run it after all dups killed
	void append_wp_nodes(const hit_idx_t& hit_idx) {
		if(!wp_.size()) return;

		// walk through the intersection path and add node points
		// of well geometry to the cell with previous intersection
		x_iterator px = x_.begin();
		//ulong facet_id;

		//wp_iterator pw = wp_.begin();
		//ulong node_idx = 0;
		ulong mesh_sz = m_.size_flat();
		for(ulong i = 0; i < wp_.size(); ++i) {
			//const well_data& wseg = pw->second;
			// upper_bound
			while(px != x_.end() && px->md < wp_[i].md())
				++px;
			// we need prev intersection
			//if(px != x_.begin()) --px;

			if(hit_idx[i] < mesh_sz)
				px = insert_wp_node(hit_idx[i], i, px);
		}

		// well path doesn't contain the end-point of trajectory
		// add it manually
		if(hit_idx[hit_idx.size() - 1] < mesh_sz)
			insert_wp_node(hit_idx[hit_idx.size() - 1], wp_.size() - 1, x_.end(), true);
	}

	//template< int N >
	void remove_dups2() {

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
			for(uint i = 0; i < D; ++i)
				dir[i] = start[i] <= finish[i] ? 0 : 1;
			// judge
			top_surv judge(dir, x_, m_.size());

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
				for(; pn != x_.end() && std::abs< double >(px->md - pn->md) < MD_TOL; ++pn, ++cnt)
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

	// directly access intersection path
	intersect_path& path() {
		return x_;
	}
	const intersect_path& path() const {
		return x_;
	}

protected:
	x_iterator insert_wp_node(ulong cell_id, ulong wseg_id, x_iterator px, bool end_point = false) {
		// initialization
		const well_data& wseg = wp_[wseg_id];
		Point where = wseg.start();
		t_float wp_md = wseg.md();
		if(end_point) {
			where = wseg.finish();
			wp_md += wseg.len();
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
				where, wseg_id, cell_id, wp_md,
				facet_id, true
			));
		else {
			//std::cout << facet_id << ' ' << inner_point_id << " facet_id != inner_point_id" << std::endl;
			well_hit_cell& x = const_cast< well_hit_cell& >(*px);
			x.is_node = true;
		}
		return px;
	}

	void check_intersection(ulong cell_id, ulong wseg_id, const Segment& well_seg) {
		typedef typename strat_t::xpoints_list xpoints_list;

		// find intersection points coord if any
		cell_data c = m_[cell_id];
		const xpoints_list& res = strat_t::precise_intersection(c, well_seg);
		// cache modified cell in order to omit polygon or triangulation recalc, used by
		// precise_intersection
		//m_.cache_cell(cell_id, c);

		// add all points to interscetion path
		for(typename xpoints_list::const_iterator px = res.begin(), end = res.end(); px != end; ++px) {
			// prepare intersection to insert
			well_hit_cell new_x(px->first, wseg_id, cell_id, calc_md(wseg_id, px->first), px->second);
			// if intersection for that point exists
			// append new intersection only if it is 'greater' than existing
			// btree: btree_multiset destroy iterators on insertion
			// account that in spatial sorting algo
			// first always insert new point
			x_.insert(new_x);
			// we expect here max 2 points with identical MD
			// check for such case
			std::pair< x_iterator, x_iterator > eqx = x_.equal_range(new_x);
			x_iterator p_secx = eqx.first;
			++p_secx;
			if(p_secx != eqx.second) {
				// calc direction vector
				dirvec_t dir;
				Point start = wp_[wseg_id].start();
				Point finish = wp_[wseg_id].finish();
				for(uint i = 0; i < D; ++i)
					dir[i] = start[i] <= finish[i] ? 0 : 1;

				// because we expect max 2 equal points, compare them
				// end remove one that is lower
				if(spatial_sort(dir, m_.size()).greater(*eqx.first, *p_secx))
					x_.erase(p_secx);
				else
					x_.erase(eqx.first);
			}
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
	intersect_path x_;
	// mesh size
	//vertex_pos_i m_size_;
	// cell IDs of where each well node is located
	hit_idx_t hit_idx_;
};

}} /* blue_sky::wpi */

#endif /* end of include guard: WPI_ALGO_XACTION_PLJRVQ8B */

