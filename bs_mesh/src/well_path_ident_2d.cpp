/// @file well_path_ident.cpp
/// @brief Well path and mesh intersection utilities in 2D
/// @author uentity
/// @version 0.1
/// @date 16.08.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "bs_mesh_stdafx.h"
#include "well_path_ident.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/intersections.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/function_objects.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/algorithm.h>
#include <CGAL/intersections.h>
#include <CGAL/Object.h>

#include "conf.h"
#include "export_python_wrapper.h"
#include "bs_mesh_grdecl.h"

#include <vector>
#include <cmath>
// DEBUG
//#include <iostream>

#define X(n) (3*n)
#define Y(n) (3*n + 1)
#define Z(n) (3*n + 2)
#define C(n, offs) (3*n + offs)
#define MD_TOL 0.000001

namespace bp = boost::python;
using namespace std;

namespace blue_sky {

namespace { 	// hide implementation details

typedef CGAL::Object                                        Object;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_2                                     Point_2;
typedef Kernel::Triangle_2                                  Triangle_2;
typedef Kernel::Segment_2                                   Segment_2;
typedef Kernel::Iso_rectangle_2                             Iso_rectangle_2;
typedef CGAL::Bbox_2                                        Bbox_2;
typedef CGAL::Polygon_2< Kernel >                           Polygon_2;
typedef std::vector<Triangle_2>                             Triangles;
typedef Triangles::iterator                                 tri_iterator;

typedef smart_ptr< bs_mesh_grdecl > sp_grd_mesh;
typedef t_ulong ulong;
typedef t_uint uint;
typedef v_float::iterator vf_iterator;

/*-----------------------------------------------------------------
 * cell description
 *----------------------------------------------------------------*/
typedef t_float vertex_pos[2];
typedef t_float cell_pos[8][3];
typedef ulong cell_pos_i[2];

struct cell_data {
	// vertex coord
	t_float* V;
	// cell facets cover with triangles
	//Triangles cover;

	// empty ctor for map
	cell_data() : V(NULL) {}
	// std ctor
	cell_data(t_float *const cell) : V(cell) {}

	void lo(vertex_pos& b) const {
		bound< std::less >(b);
	}

	void hi(vertex_pos& b) const {
		bound< std::greater >(b);
	}

	cell_pos& cpos() {
		return reinterpret_cast< cell_pos& >(*V);
	}

	const cell_pos& cpos() const {
		return reinterpret_cast< const cell_pos& >(*V);
	}

	Bbox_2 bbox() const {
		vertex_pos p1, p2;
		lo(p1); hi(p2);
		return Bbox_2(p1[0], p1[1], p2[0], p2[1]);
	}

	Polygon_2 polygon() const {
		// 2D upper plane of cell
		Point_2 points[] = {
			Point_2(V[0], V[1]),  Point_2(V[3], V[4]),
			Point_2(V[9], V[10]), Point_2(V[6], V[7])
		};
		return Polygon_2(points, points + 4);
	}


private:
	template< template< class > class pred >
	void bound(vertex_pos& b) const {
		pred< t_float > p = pred< t_float >();
		const cell_pos& cV = cpos();
		for(uint i = 0; i < 2; ++i) {
			t_float c = cV[0][i];
			for(uint j = 1; j < 4; ++j) {
				if(p(cV[j][i], c))
					c = cV[j][i];
			}
			b[i] = c;
		}
	}
};
typedef st_smart_ptr< cell_data > sp_cell_data;

// storage for representing mesh
typedef std::map< t_ulong, cell_data > trimesh;
typedef trimesh::const_iterator trim_iterator;

/*-----------------------------------------------------------------
 * represent rectangular part of mesh with splitting support
 *----------------------------------------------------------------*/
// x_last = last_element + 1 = x_size
// y_last = last_element + 1 = y_size
struct mesh_part {
	typedef set< mesh_part > container_t;

	mesh_part(const trimesh& m, ulong nx, ulong ny)
		: x_first(0), x_last(nx)
		, y_first(0), y_last(ny)
		, m_(m), nx_(nx), ny_(ny)
	{}

	void init(ulong x1, ulong x2, ulong y1, ulong y2) {
		x_first = x1; x_last = x2;
		y_first = y1; y_last = y2;

		// sanity checks
		x_first = min(x_first, nx_ - 1);
		x_last   = min(x_last  , nx_);
		y_first = min(y_first, ny_ - 1);
		y_last   = min(y_last  , ny_);

		x_last = max(x_first + 1, x_last);
		y_last = max(y_first + 1, y_last);
	}

	void init(const vector< ulong >& cell_idx) {
		ulong x1, x2;
		ulong y1, y2;
		// search for bounds
		for(ulong i = 0; i < cell_idx.size(); ++i) {
			ulong y = cell_idx[i] / nx_;
			ulong x = cell_idx[i] - y * nx_;
			if(i == 0) {
				x1 = x2 = x;
				y1 = y2 = y;
			}
			else {
				x1 = min(x1, x); x2 = max(x2, x);
				y1 = min(y1, y); y2 = max(y2, y);
			}
		}
		// +1 for x_last and y_last
		// finally usual init
		init(x1, ++x2, y1, ++y2);
	}

	ulong side_len(int dim) const {
		if(!dim)
			return x_last - x_first;
		else
			return y_last - y_first;
	}

	ulong size() const {
		return side_len(0) * side_len(1);
	}

	trim_iterator ss_iter(ulong x, ulong y) {
		if(x >= side_len(0) || y >= side_len(1))
			return m_.end();
		ulong idx = (y_first + y) * nx_ + x_first + x;
		return m_.find(idx);
	}

	Iso_rectangle_2 bbox() const {
		// calc idx of "upper left" cell
		const ulong start_idx = y_first * nx_ + x_first;
		const ulong end_idx = (y_last - 1) * nx_ + x_last - 1;

		vertex_pos lo, hi;
		mesh_ss(start_idx).lo(lo);
		mesh_ss(end_idx).hi(hi);
		return Iso_rectangle_2(Point_2(lo[0], lo[1]), Point_2(hi[0], hi[1]));
	}

	container_t divide() const {
		// resulting split
		container_t res;
		insert_iterator< container_t > ii(res, res.begin());

		// split points
		ulong x_div[3], y_div[3];
		x_div[0] = x_first; x_div[2] = x_last;
		x_div[1] = x_first + (side_len(0) >> 1);
		y_div[0] = y_first; y_div[2] = y_last;
		y_div[1] = y_first + (side_len(1) >> 1);

		// make splitting only if split containt more than 1 cell
		//ulong x_len, y_len;
		for(ulong i = 0; i < 2; ++i) {
			if(y_div[i + 1] - y_div[i] == 0) continue;
			for(ulong j = 0; j < 2; ++j) {
				if(x_div[j + 1] - x_div[j] == 0) continue;
				*ii++ = mesh_part(m_, nx_, ny_, x_div[j], x_div[j + 1], y_div[i], y_div[i + 1]);
			}
		}
		return res;
	}

	// for sorted containers
	bool operator <(const mesh_part& rhs) const {
		return
		x_first == rhs.x_first ?
			(x_last == rhs.x_last ?
				(y_first == rhs.y_first ?
					y_last < rhs.y_last :
				y_first < rhs.y_first) :
			x_last < rhs.x_last) :
		x_first < rhs.x_first;
	}

	// public members
	ulong x_first, x_last;
	ulong y_first, y_last;

private:
	const trimesh& m_;
	ulong nx_, ny_;

	mesh_part(const trimesh& m, ulong nx, ulong ny,
			ulong x1, ulong x2,
			ulong y1, ulong y2)
		: x_first(x1), x_last(x2)
		, y_first(y1), y_last(y2)
		, m_(m), nx_(nx), ny_(ny)
	{}

	const cell_data& mesh_ss(ulong idx) const {
		// idx SHOULD BE IN MESH!
		return m_.find(idx)->second;
	}
};

/*-----------------------------------------------------------------
 * well description
 *----------------------------------------------------------------*/
struct well_data {
	// segment begin, end and md in raw vector
	t_float* W;
	double md_;

	//empty ctor for map
	well_data() : W(NULL), md_(0) {}
	//std ctor
	well_data(t_float *const segment, double md) : W(segment), md_(md) {}

	t_float md() const {
		return md_;
		//return W[3];
	}

	vertex_pos& cstart() {
		return reinterpret_cast< vertex_pos& >(*W);
	}

	const vertex_pos& cstart() const {
		return reinterpret_cast< const vertex_pos& >(*W);
	}

	Point_2 start() const {
		return Point_2(W[0], W[1]);
	}
	Point_2 finish() const {
		return Point_2(W[4], W[5]);
	}

	//vertex_pos& cend() {
	//	return reinterpret_cast< vertex_pos& >(W + 4);
	//}

	//const vertex_pos& cend() const {
	//	return reinterpret_cast< const vertex_pos& >(W + 4);
	//}

	Segment_2 segment() const {
		return Segment_2(start(), finish());
	}

	Bbox_2 bbox() const {
		return segment().bbox();
	}

	double len() const {
		return std::sqrt(segment().squared_length());
	}
};

typedef std::map< ulong, well_data > well_path;
typedef well_path::iterator wp_iterator;
typedef well_path::const_iterator cwp_iterator;

/*-----------------------------------------------------------------
 * Box description
 *----------------------------------------------------------------*/
// structure to help identify given boxes
class box_handle {
public:
	enum {
		CELL_BOX,
		WELL_BOX
	};

	virtual int type() const = 0;

protected:
	template< class fish_t, class = void >
	struct fish2box_t {
		// default value
		enum { type = CELL_BOX };
	};
	// overload for well path
	template< class unused >
	struct fish2box_t< wp_iterator, unused > {
		enum { type = WELL_BOX };
	};
};
// pointer is really stored as box handle
typedef st_smart_ptr< box_handle > sp_bhandle;

template< class fish >
class box_handle_impl : public box_handle {
public:
	typedef fish fish_t;

	box_handle_impl(const fish_t& f) : f_(f) {}

	int type() const {
		return box_handle::fish2box_t< fish_t >::type;
	}

	fish_t data() const {
		return f_;
	}

private:
	fish_t f_;
};
// handy typedefs
typedef box_handle_impl< trim_iterator > cell_box_handle;
typedef box_handle_impl< wp_iterator > well_box_handle;

// box intersections storage
typedef CGAL::Box_intersection_d::Box_with_handle_d< double, 2, sp_bhandle > Box;

/*-----------------------------------------------------------------
 * intersections description
 *----------------------------------------------------------------*/
struct well_hit_cell {
	// point of intersection
	Point_2 where;
	// what segment of well
	wp_iterator seg;
	// interseect with which cell
	trim_iterator cell;
	// depth along well in point of intersection
	t_float md;
	// cell facet
	uint facet;
	// is that point a node?
	bool is_node;

	well_hit_cell() {}
	// std ctor
	well_hit_cell(const Point_2& where_, const wp_iterator& seg_,
		const trim_iterator& cell_, t_float md_, uint facet_, bool is_node_ = false)
		//double z_ = 0)
		: where(where_), seg(seg_), cell(cell_), md(md_), facet(facet_),
		is_node(is_node_)
	{}
	// for searching
	well_hit_cell(t_float md_) : md(md_) {}

	// hit points ordered by md
	bool operator <(const well_hit_cell& rhs) const {
		return md < rhs.md;
		//if(md < rhs.md)
		//	return true;
		//else if(md == rhs.md) {
		//	if(seg->first < rhs.seg->first)
		//		return true;
		//	else if(seg->first == rhs.seg->first)
		//		return cell->first < rhs.cell->first;
		//}
		//return false;
	}
};

// storage of intersection points
typedef std::multiset< well_hit_cell > intersect_path;

/*-----------------------------------------------------------------
 * callback functor that triggers on intersecting bboxes
 *----------------------------------------------------------------*/
class intersect_action {
public:
	typedef intersect_path::iterator x_iterator;

	// traits for removing duplicates in X and Y direction
	struct dup_traits_x {
		// dimesion to operate onto
		enum { dim_id = 0 };
		// specify IDs of crossing sides of cell in positive direction
		enum { cross_1st = 1, cross_2nd = 3 };
		enum { axe_facet1 = 0, axe_facet2 = 2 };
	};

	struct dup_traits_y {
		// dimesion to operate onto
		enum { dim_id = 1 };
		// specify IDs of crossing sides of cell in positive direction
		enum { cross_1st = 2, cross_2nd = 0 };
		enum { axe_facet1 = 1, axe_facet2 = 3 };
	};

	template< int N >
	struct spatial_sort {
		typedef int dirvec_t[N];
		typedef intersect_path::iterator x_iterator;

		spatial_sort(const dirvec_t& dir, const intersect_action& A)
			: dir_(dir), A_(A)
		{}

		spatial_sort(const spatial_sort& rhs)
			: dir_(rhs.dir_), A_(rhs.A_)
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
			cell_pos_i c1, c2;
			A_.decode_cell_pos(r1->cell->first, c1);
			A_.decode_cell_pos(r2->cell->first, c2);

			//bool res = false;
			for(uint i = 0; i < N; ++i) {
				if(greater(i, c1[i], c2[i]))
					return true;
			}
			return false;
		}

		bool greater(uint ndim, ulong v1, ulong v2) const {
			return dir_[ndim] == 0 ? v1 > v2 : v2 > v1;
		}

		const dirvec_t& dir_;
		const intersect_action& A_;
	};

	// ctor
	intersect_action(trimesh& mesh, well_path& wp, intersect_path& X, ulong nx, ulong ny)
		: m_(mesh), wp_(wp), x_(X), nx_(nx), ny_(ny)
	{}

 /* nodes layout
     *                             X
     *                    0+-------+1
     *                    /|     / |
     *                  /  |   /   |
     *               2/-------+3   |
     *              Y |   4+--|----+5
     *                |   /Z  |   /
     *                | /     | /
     *              6 /-------/7
*/
/*  facets layout (nodes - plane id)
 *  0-1-2-3 - 0
 *  0-1-4-5 - 1
 *  4-5-6-7 - 2
 *  2-3-6-7 - 3
 *  0-2-4-6 - 4
 *  1-3-5-7 - 5
 *  inside cell - 6
*/

	static double distance(const Point_2& p1, const Point_2& p2) {
		return std::sqrt(Segment_2(p1, p2).squared_length());
	}

	static double calc_md(const wp_iterator& fish, const Point_2& target) {
		// walk all segments before the last one;
		double md = fish->second.md();
		// append tail
		const t_float* W = fish->second.W;
		md += distance(Point_2(W[0], W[1]), target);
		return md;
	}

	void operator()(const Box& bc, const Box& bw) {
		//std::cout << "intersection detected!" << std::endl;
		//return;

		trim_iterator cell_fish = static_cast< cell_box_handle* >(bc.handle().get())->data();
		wp_iterator well_fish = static_cast< well_box_handle* >(bw.handle().get())->data();
		const cell_data& c = cell_fish->second;
		const well_data& w = well_fish->second;

		// obtain well segment
		const Segment_2& s = w.segment();
		// and cell polygon
		const Polygon_2& p = c.polygon();

		// intersect well segment with each segment of polygon
		for(ulong i = 0; i < 4; ++i) {
			Object xres = CGAL::intersection(p.edge(i), s);
			// in 99% of cases we should get a point of intersection
			if(const Point_2* xpoint = CGAL::object_cast< Point_2 >(&xres))
				x_.insert(well_hit_cell(
					*xpoint, well_fish, cell_fish,
					calc_md(well_fish, *xpoint),
					i
				));
			else if(const Segment_2* xseg = CGAL::object_cast< Segment_2 >(&xres)) {
				// in rare 1% of segment lying on the facet, add begin and end of segment
				x_.insert(well_hit_cell(
					xseg->source(), well_fish, cell_fish,
					calc_md(well_fish, xseg->source()),
					i
				));
				x_.insert(well_hit_cell(
					xseg->target(), well_fish, cell_fish,
					calc_md(well_fish, xseg->target()),
					i
				));
			}
		}
	}

	// run it after all dups killed
	void append_wp_nodes(const vector< ulong >& hit_idx) {
		if(!wp_.size()) return;

		// walk through the intersection path and add node points
		// of well geometry to the cell with previous intersection
		intersect_path::iterator px = x_.begin();
		//ulong facet_id;

		wp_iterator pw = wp_.begin();
		ulong node_idx = 0;
		for(wp_iterator end = wp_.end(); pw != end; ++pw, ++node_idx) {
			//const well_data& wseg = pw->second;
			// lower_bound
			while(px != x_.end() && px->md < pw->second.md())
				++px;

			px = insert_wp_node(hit_idx[node_idx], pw, px);
		}

		// well path doesn't contain the end-point of trajectory
		// add it manually
		insert_wp_node(hit_idx[node_idx], --pw, x_.end(), true);
	}

	template< class dup_traits >
	void remove_dups(const dup_traits& t) {
		typedef intersect_path::iterator x_iterator;

		// helper to decide which of consequent cross-points whould we keep
		struct who_survive {
			who_survive(int dir, intersect_path& x)
				: surv_id_(dup_traits::cross_2nd), kill_id_(dup_traits::cross_1st), x_(x)
			{
				// in positive direction 2nd point survives
				// in negative - first point
				if(dir)
					swap(surv_id_, kill_id_);
			}

			x_iterator operator()(const x_iterator& p1, const x_iterator& p2) {
				if(p1->facet == surv_id_ && p2->facet == kill_id_) {
					x_.erase(p2);
					return p1;
				}
				else if(p1->facet == kill_id_ && p2->facet == surv_id_) {
					x_.erase(p1);
					return p2;
				}
				return p1;
			}

			uint surv_id_, kill_id_;
			intersect_path& x_;
		};

		// sanity check
		if(x_.size() < 2) return;

		// position on first intersection
		x_iterator px = x_.begin();
		// walk the nodes and determine direction of trajectory
		//int dir;
		const int dim_id = dup_traits::dim_id;
		double max_md;
		for(wp_iterator pw = wp_.begin(), end = wp_.end(); pw != end; ++pw) {
			// identify direction
			const well_data& seg = pw->second;
			//Point_2 seg_end = seg.finish();
			who_survive judge(
				seg.start()[dim_id] < seg.finish()[dim_id] ? 0 : 1,
				x_
			);

			// remove dups lying on current well segment
			max_md = seg.md() + seg.len();
			for(; px != x_.end() && px->md <= max_md; ++px) {
				// skeep well node points if any
				if(px->facet == 4)
					continue;

				// check if next intersection coincides with current one
				// remove dup
				x_iterator pn = px;
				++pn;
				if(pn != x_.end() && abs(px->md - pn->md) < MD_TOL)
					px = judge(px, pn);
			}
		}
	}

	template< int N >
	void remove_dups2() {
		typedef int dirvec_t[N];
		typedef intersect_path::iterator x_iterator;

		struct top_surv {
			typedef set< x_iterator, spatial_sort< N > > spat_storage_t;
			typedef typename spat_storage_t::iterator spat_iterator;

			top_surv(const dirvec_t& dir, intersect_path& x, const intersect_action& A)
				: dir_(dir), x_(x), A_(A)
			{}

			x_iterator operator()(x_iterator from, x_iterator to) {
				spat_storage_t sr(spatial_sort< N >(dir_, A_));

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
			const intersect_action& A_;
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
			const well_data& seg = pw->second;
			// calc direction vector
			Point_2 start = seg.start();
			Point_2 finish = seg.finish();
			for(uint i = 0; i < 2; ++i)
				dir[i] = start[i] <= finish[i] ? 0 : 1;
			// judge
			top_surv judge(dir, x_, *this);

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

	// X-Y order!
	void decode_cell_pos(ulong cell_id, cell_pos_i& res) const {
		res[1] = cell_id / nx_;
		res[0] = cell_id - res[1] * nx_;
	}

	spv_float export_1d() const {
		spv_float res = BS_KERNEL.create_object(v_float::bs_type());
		res->resize(x_.size() * 6);
		vf_iterator pr = res->begin();

		for(intersect_path::const_iterator px = x_.begin(), end = x_.end(); px != end; ++px) {
			// cell num
			*pr++ = px->cell->first;
			// MD
			*pr++ = px->md;
			// intersection point
			*pr++ = px->where.x();
			*pr++ = px->where.y();
			//*pr++ = px->z;
			// facet id
			*pr++ = px->facet;
			// node flag
			*pr++ = px->is_node;
		}

		return res;
	}

	vector< ulong > where_is_point(vector< Point_2 > points) const {
		// start with full mesh
		// and divide it until we come to only one cell

		// mesh partition stored here
		typedef mesh_part::container_t parts_container;
		typedef mesh_part::container_t::iterator part_iterator;
		parts_container parts;
		parts.insert(mesh_part(m_, nx_, ny_));

		// found cell_ids stored here
		vector< ulong > res(points.size(), -1);
		//ulong cell_id;
		while(parts.size()) {
			// split every part
			parts_container leafs;
			for(part_iterator p = parts.begin(), end = parts.end(); p != end; ++p) {
				parts_container kids = p->divide();
				leafs.insert(kids.begin(), kids.end());
			}

			// collection of points inside current partition
			list< ulong > catched_points;
			// process each leaf and find points inside it
			for(part_iterator l = leafs.begin(), end = leafs.end(); l != end; ) {
				const Iso_rectangle_2& cur_rect = l->bbox();
				catched_points.clear();
				for(ulong i = 0; i < points.size(); ++i) {
					// skip already found points
					if(res[i] < m_.size()) continue;
					// check that point lies inside this part
					if(!cur_rect.has_on_unbounded_side(points[i]))
						catched_points.insert(catched_points.begin(), i);
				}

				// if this part don't contain any points - remove it
				// if box contains only 1 cell - test if cell poly contains given points
				if(!catched_points.size())
					leafs.erase(l++);
				else if(l->size() == 1) {
					ulong cell_id = l->y_first * nx_ + l->x_first;
					Polygon_2 cell_poly = m_[cell_id].polygon();
					for(list< ulong >::iterator pp = catched_points.begin(), cp_end = catched_points.end(); pp != cp_end; ++pp) {
						if(!cell_poly.has_on_unbounded_side(points[*pp]))
							res[*pp] = cell_id;
					}
					leafs.erase(l++);
				}
				else ++l;
			}

			// leafs become the new start point for further division
			parts = leafs;
		}
		return res;
	}

	ulong where_is_point(Point_2 point) const {
		return where_is_point(vector< Point_2 >(1, point))[0];
	}

private:

	x_iterator insert_wp_node(ulong cell_id, wp_iterator pw, x_iterator px, bool end_point = false) {
		// initialization
		const well_data& wseg = pw->second;
		Point_2 where = wseg.start();
		t_float wp_md = wseg.md();
		if(end_point) {
			where = wseg.finish();
			wp_md += wseg.len();
			//wp_md = px->md + distance(px->where, where);
		}

		// check if current or prev intersection match with node
		uint facet_id = 4;
		if(px != x_.end() && abs(px->md - wp_md) < MD_TOL)
			facet_id = px->facet;
		else if(px != x_.begin()) {
			--px;
			if(abs(px->md - wp_md) < MD_TOL)
				facet_id = px->facet;
		}

		// if node point coinside with existing intersection
		// then just change is_node flag (dont't affect ordering)
		// otherwise insert new intersection point
		if(facet_id == 4)
			px = x_.insert(well_hit_cell(
				where, pw, m_.find(cell_id), wp_md,
				facet_id, true
			));
		else {
			well_hit_cell& x = const_cast< well_hit_cell& >(*px);
			x.is_node = true;
		}
		return px;
	}

	// mesh
	trimesh& m_;
	// well path
	well_path& wp_;
	// well-mesh intersection path
	intersect_path& x_;
	// mesh dimensions
	ulong nx_, ny_;
};

// helper to create initial cell_data for each cell
spv_float coord_zcorn2trimesh(t_long nx, t_long ny, spv_float coord, spv_float zcorn, trimesh& res) {
	// build mesh_grdecl around given mesh
	sp_grd_mesh grd_src = BS_KERNEL.create_object(bs_mesh_grdecl::bs_type());
	grd_src->init_props(nx, ny, coord, zcorn);
	//t_long nz = (zcorn->size() / nx / ny) >> 3;
	
	// obtain coordinates for all vertices of all cells
	spv_float tops = grd_src->calc_cells_vertices_xyz();
	v_float::iterator pv = tops->begin();

	// fill trimesh with triangles corresponding to each cell
	// use only first plane of cells
	ulong n_cells = ulong(nx * ny);
	for(ulong i = 0; i < n_cells; ++i) {
		res[i] = cell_data(&*pv);
		pv += 3*8;
	}

	return tops;
}

} 	// eof hidden namespace

/*-----------------------------------------------------------------
 * implementation of main routine
 *----------------------------------------------------------------*/
spv_float well_path_ident_2d(t_long nx, t_long ny, spv_float coord, spv_float zcorn,
	spv_float well_info, bool include_well_nodes)
{
	// 1) calculate mesh nodes coordinates and build initial trimesh
	trimesh M;
	spv_float tops;
	tops = coord_zcorn2trimesh(nx, ny, coord, zcorn, M);
	//t_long nz = (zcorn->size() / nx / ny) >> 3;

	// 2) create well path description and
	// bounding boxes for line segments representing well trajectory
	ulong well_node_num = well_info->size() >> 2;
	if(well_node_num < 2) return spv_float();

	// storage
	well_path W;
	vector< Box > well_boxes(well_node_num - 1);
	// build array of well nodes as Point_2
	vector< Point_2 > wnodes(well_node_num);

	// walk along well
	v_float::iterator pw = well_info->begin();
	double md = 0;
	for(ulong i = 0; i < well_node_num - 1; ++i) {
		//well_boxes.push_back(Box(
		//	lo, hi,
		//	new well_box_handle( W.insert(make_pair(i, pw)).first )
		//));

		// calc md
		if(i > 0)
			md += W[i - 1].len();

		well_data wd(pw, md);
		wnodes[i] = wd.start();

		// make bbox
		well_boxes[i] = Box(
			wd.bbox(),
			new well_box_handle(W.insert(make_pair(i, wd)).first)
		);

		pw += 4;
	}
	// put last node to array
	wnodes[well_node_num - 1] = W[well_node_num - 2].finish();

	// 3) find where each node of well is located
	// to restrict search area
	// intersections storage
	intersect_path X;
	intersect_action A(M, W, X, nx, ny);
	const vector< ulong >& hit_idx = A.where_is_point(wnodes);
	// DEBUG
	//cout << "well nodes hit_idx built" << endl;
	// create part of mesh to process based on these cells
	mesh_part hot_mesh(M, nx, ny);
	hot_mesh.init(hit_idx);

	// create bounding box for each cell in given mesh
	std::vector< Box > mesh_boxes(hot_mesh.size());
	ulong cnt = 0;
	trim_iterator pm;
	for(ulong j = 0; j < hot_mesh.side_len(1); ++j) {
		for(ulong i = 0; i < hot_mesh.side_len(0); ++i) {
			pm = hot_mesh.ss_iter(i, j);
			const cell_data& d = pm->second;
			mesh_boxes[cnt++] = Box(d.bbox(), new cell_box_handle(pm));
		}
	}

	// Create the corresponding vector of pointers to cells bounding boxes
	//std::vector< cell_box* > mesh_boxes_ptr;
	//for(std::vector< cell_box >::iterator i = mesh_boxes.begin(); i != mesh_boxes.end(); ++i)
	//	mesh_boxes_ptr.push_back(&*i);

	// Run the intersection algorithm with all defaults on the
	// indirect pointers to cell bounding boxes. Avoids copying the boxes
	CGAL::box_intersection_d(
		mesh_boxes.begin(), mesh_boxes.end(),
		well_boxes.begin(), well_boxes.end(),
		A
	);
	//cout << "facet intersections" << endl;

	// append well path nodes
	if(include_well_nodes)
		A.append_wp_nodes(hit_idx);
	//cout << "wwll path nodes added" << endl;

	// remove duplicates in X and Y directions
	//A.remove_dups(intersect_action::dup_traits_x());
	//A.remove_dups(intersect_action::dup_traits_y());
	A.remove_dups2< 2 >();

	return A.export_1d();
}

}	// eof blue-sky namespace


