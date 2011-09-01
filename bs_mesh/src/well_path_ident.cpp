/// @file well_path_ident.cpp
/// @brief Well path and mesh intersection utilities
/// @author uentity
/// @version 0.1
/// @date 05.07.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "bs_mesh_stdafx.h"
#include "well_path_ident.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/intersections.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Bbox_3.h>
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
#include <boost/array.hpp>
// DEBUG
//#include <iostream>

#define X(n) (3*n)
#define Y(n) (3*n + 1)
#define Z(n) (3*n + 2)
#define C(n, offs) (3*n + offs)

namespace bp = boost::python;
using namespace std;

namespace blue_sky {

namespace { 	// hide implementation details

typedef CGAL::Object                                        Object;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                                     Point_3;
typedef Kernel::Triangle_3                                  Triangle_3;
typedef Kernel::Segment_3                                   Segment_3;
typedef CGAL::Bbox_3                                        Bbox_3;
typedef Kernel::Iso_cuboid_3                                Iso_cuboid_3;
typedef std::vector<Triangle_3>                             Triangles;
typedef Triangles::iterator                                 tri_iterator;

typedef smart_ptr< bs_mesh_grdecl > sp_grd_mesh;
typedef t_ulong ulong;
typedef t_uint uint;
typedef v_float::iterator vf_iterator;

/*-----------------------------------------------------------------
 * strategy
 *----------------------------------------------------------------*/
// dimens num, cell vertex num
enum { D = 3, CVN = 8 };
// boost::array with opertator= and ctor with elements init
template< class T >
class stat_array : public boost::array< T, D > {
public:
	typedef boost::array< t_float, D > base_t;
	// propagate typedefs
	typedef typename base_t::value_type             value_type;
	typedef typename base_t::iterator               iterator;
	typedef typename base_t::const_iterator         const_iterator;
	typedef typename base_t::reverse_iterator       reverse_iterator;
	typedef typename base_t::const_reverse_iterator const_reverse_iterator;
	typedef typename base_t::reference              reference;
	typedef typename base_t::const_reference        const_reference;
	typedef typename base_t::size_type              size_type;
	typedef typename base_t::difference_type        difference_type;

	using base_t::begin;
	using base_t::end;
	using base_t::size;
	using base_t::elems;

	// empty ctor
	stat_array() : base_t() {
		// ensure all elems are filled with zero
		fill(begin(), end(), value_type());
	}

	// ctor accepting C-array
	stat_array(const value_type(&data)[D]) : base_t() {
		copy(&data[0], &data[D], begin());
	}

	// hack-like but useful ctor with per-element init
	stat_array(int N, ...) : base_t() {
		assert(D < N && "vertex_pos overflow!");
		va_list arg_list;
		va_start(arg_list, N);
		for(uint i = 0; i < N; ++i) {
			elems[i] = va_arg(arg_list, value_type);
		}
		va_end(arg_list);
	}

	// assigning arrays
	template< class R >
	stat_array& operator=(const stat_array< R >& rhs) {
		copy(rhs.begin(), rhs.begin() + min(size(), rhs.size()), begin());
		return *this;
	}
};

// strategy base types
//typedef stat_array< t_float > vertex_pos;
//typedef vertex_pos< ulong >   vertex_pos_i;
//typedef vertex_pos            cell_pos[CVN];

typedef t_float vertex_pos[D];
typedef ulong   vertex_pos_i[D];
typedef t_float cell_pos[CVN][D];
// assign for c arrays
// fun with returning reference to array :)
template< class T >
T (&ca_assign(T (&lhs)[D], const T (&rhs)[D]))[D] {
	copy(&rhs[0], &rhs[D], &lhs[0]);
	return lhs;
}

template< class T >
T (&ca_assign(T (&lhs)[D], const T& v))[D] {
	fill(&lhs[0], &lhs[D], v);
	return lhs;
}

// X-Y-Z order!
void decode_cell_id(ulong id, vertex_pos_i& res, const vertex_pos_i& m_size) {
	//vertex_pos_i res;
	res[2] = id / (m_size[0] * m_size[1]);
	res[1] = (id - res[2] * m_size[0] * m_size[1]) / m_size[0];
	res[0] = id - m_size[0]*(res[2] * m_size[1] + res[1]);
}

ulong encode_cell_id(const vertex_pos_i& p, const vertex_pos_i& m_size) {
	return p[0] + m_size[0] * (p[1] + p[2] * m_size[1]);
}

Bbox_3 vertex_pos2bbox(const vertex_pos& lo, const vertex_pos& hi) {
	return Bbox_3(lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]);
}

Iso_cuboid_3 vertex_pos2rect(const vertex_pos& lo, const vertex_pos& hi) {
	return Iso_cuboid_3(Point_3(lo[0], lo[1], lo[2]), Point_3(hi[0], hi[1], hi[2]));
}

/*-----------------------------------------------------------------
 * cell description
 *----------------------------------------------------------------*/
struct cell_data {
	// vertex coord
	t_float* V;
	// cell facets cover with triangles
	Triangles cover;

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

	Bbox_3 bbox() const {
		vertex_pos p1, p2;
		lo(p1); hi(p2);
		return vertex_pos2bbox(p1, p2);
	}

private:
	template< template< class > class pred >
	void bound(vertex_pos& b) const {
		pred< t_float > p = pred< t_float >();
		const cell_pos& cV = cpos();
		for(uint i = 0; i < D; ++i) {
			t_float c = cV[0][i];
			for(uint j = 1; j < CVN; ++j) {
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
typedef trimesh::iterator trim_iterator;
typedef trimesh::const_iterator ctrim_iterator;

/*-----------------------------------------------------------------
 * represent rectangular part of mesh with splitting support
 *----------------------------------------------------------------*/
// x_last = last_element + 1 = x_size
// y_last = last_element + 1 = y_size
struct mesh_part {
	typedef set< mesh_part > container_t;

	mesh_part(const trimesh& m, const vertex_pos_i& mesh_size)
		: m_(m)
	{
		ca_assign(lo, ulong(0));
		ca_assign(hi, mesh_size);
		ca_assign(m_size_, mesh_size);
	}

	void init(const vertex_pos_i& lower, const vertex_pos_i& upper) {
		ca_assign(lo, lower);
		ca_assign(hi, upper);

		// sanity checks
		for(uint i = 0; i < D; ++i) {
			lo[i] = min(lo[i], m_size_[i] - 1);
			hi[i] = min(hi[i], m_size_[i]);
			hi[i] = max(lo[i] + 1, hi[i]);
		}
	}

	void init(const vector< ulong >& cell_idx) {
		vertex_pos_i lower, upper, p;
		// search for bounds
		for(ulong i = 0; i < cell_idx.size(); ++i) {
			decode_cell_id(cell_idx[i], p, m_size_);
			if(i == 0)
				ca_assign(lower, ca_assign(upper, p));
			else {
				for(uint i = 0; i < D; ++i) {
					lower[i] = min(lower[i], p[i]);
					upper[i] = max(upper[i], p[i]);
				}
			}
		}
		// add +1 to upper
		for(uint i = 0; i < D; ++i)
			++upper[i];
		// finally usual init
		init(lower, upper);
	}

	ulong side_len(uint dim) const {
		if(dim < D)
			return hi[dim] - lo[dim];
		else
			return hi[D - 1] - lo[D - 1];
	}

	ulong size() const {
		ulong res = 1;
		for(uint i = 0; i < D; ++i)
			res *= side_len(i);
		return res;
	}

	ctrim_iterator ss_iter(const vertex_pos_i& offset) {
		vertex_pos_i cell;
		ca_assign(cell, lo);
		for(uint i = 0; i < D; ++i)
			cell[i] += offset[i];
		return m_.find(encode_cell_id(cell, m_size_));
	}

	Iso_cuboid_3 bbox() const {
		vertex_pos lo_pos, hi_pos;
		mesh_ss(lo).lo(lo_pos);
		// last = hi - 1
		vertex_pos_i last;
		ca_assign(last, hi);
		transform(&last[0], &last[D], &last[0], bind2nd(minus< ulong >(), 1));
		mesh_ss(last).hi(hi_pos);
		return vertex_pos2rect(lo_pos, hi_pos);
	}

	container_t divide() const {
		// resulting split
		container_t res;
		insert_iterator< container_t > ii(res, res.begin());

		// split points
		vertex_pos_i split_p[3];
		ca_assign(split_p[0], lo);
		ca_assign(split_p[2], hi);
		// middle
		for(uint i = 0; i < D; ++i)
			split_p[1][i] = hi[i] >> 1;

		// TODO: make it dimens-independent
		// make splitting only if split containt more than 1 cell
		vertex_pos_i spl_lo, spl_hi;
		for(ulong z = 0; z < 2; ++z) {
			//if(split_p[z + 1][2] - split_p[z][2] == 0) continue;
			spl_lo[2] = split_p[z][2]; spl_hi[2] = split_p[z + 1][2];
			for(ulong y = 0; y < 2; ++y) {
				spl_lo[1] = split_p[y + 1][1]; spl_hi[1] = split_p[y][1];
					for(ulong x = 0; x < 2; ++x) {
						spl_lo[1] = split_p[y + 1][0]; spl_hi[1] = split_p[y][0];
						// check if any side is zero
						ulong sz = 1;
						for(uint i = 0; i < D; ++i)
							sz *= spl_hi[i] - spl_lo[i];
						// add new child cell
						if(sz)
							*ii++ = mesh_part(m_, m_size_, spl_lo, spl_hi);
					}
			}
		}
		return res;
	}

	// for sorted containers
	bool operator <(const mesh_part& rhs) const {
		for(uint i = 0; i < D; ++i) {
			if(lo[i] < rhs.lo[i])
				return true;
			else if(lo[i] > rhs.lo[i])
				return false;
			else if(hi[i] < rhs.hi[i])
				return true;
			else if(hi[i] > rhs.hi[i])
				return false;
		}
		return false;
	}

	// public members
	vertex_pos_i lo, hi;

private:
	const trimesh& m_;
	vertex_pos_i m_size_;

	mesh_part(const trimesh& m, const vertex_pos_i& mesh_size,
			const vertex_pos_i& first_,
			const vertex_pos_i& last_)
		: m_(m)
	{
		ca_assign(lo, first_);
		ca_assign(hi, last_);
		ca_assign(m_size_, mesh_size);
	}

	const cell_data& mesh_ss(ulong idx) const {
		// idx SHOULD BE IN MESH!
		return m_.find(idx)->second;
	}

	const cell_data& mesh_ss(const vertex_pos_i& idx) const {
		return m_.find(encode_cell_id(idx, m_size_))->second;
	}
};

/*-----------------------------------------------------------------
 * well description
 *----------------------------------------------------------------*/
struct well_data {
	// segment begin, end and md in raw vector
	t_float* W;

	//empty ctor for map
	well_data() : W(NULL) {}
	//std ctor
	well_data(t_float *const segment) : W(segment) {}

	t_float md() const {
		return W[3];
	}

	vertex_pos& cstart() {
		return reinterpret_cast< vertex_pos& >(*W);
	}

	const vertex_pos& cstart() const {
		return reinterpret_cast< const vertex_pos& >(*W);
	}

	//vertex_pos& cend() {
	//	return reinterpret_cast< vertex_pos& >(W + 4);
	//}

	//const vertex_pos& cend() const {
	//	return reinterpret_cast< const vertex_pos& >(W + 4);
	//}

	Segment_3 segment() const {
		return Segment_3(
			Point_3(W[0], W[1], W[2]),
			Point_3(W[4], W[5], W[6])
		);
	}

	Bbox_3 bbox() const {
		return segment().bbox();
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

//class well_box_handle : public box_handle {
//public:
//	well_box_handle(v_float::iterator p) : p_(p) {}
//
//	int type() const {
//		return box_handle::WELL_BOX;
//	}
//
//	well_path::iterator data() const {
//		return p_;
//	}
//
//private:
//	v_float::iterator p_;
//};

typedef CGAL::Box_intersection_d::Box_with_handle_d< double, 3, sp_bhandle > Box;
//typedef CGAL::Box_intersection_d::Box_with_handle_d< double, 3, trim_iterator > cell_box;
//typedef CGAL::Box_intersection_d::Box_with_handle_d< double, 3, typename v_float::iterator > well_box;

/*-----------------------------------------------------------------
 * intersections description
 *----------------------------------------------------------------*/
struct well_hit_cell {
	// point of intersection
	Point_3 where;
	// what segment of well
	wp_iterator seg;
	// interseect with which cell
	trim_iterator cell;
	// depth along well in point of intersection
	t_float md;
	// cell facet
	uint facet;


	well_hit_cell() {}
	well_hit_cell(const Point_3& where_, const wp_iterator& seg_,
		const trim_iterator& cell_, t_float md_, uint facet_)
		: where(where_), seg(seg_), cell(cell_), md(md_), facet(facet_)
	{}

	// hit points ordered first by md
	// next by well segment
	// and at last by cell number
	bool operator <(const well_hit_cell& rhs) const {
		if(md < rhs.md)
			return true;
		else if(md == rhs.md) {
			if(seg->first < rhs.seg->first)
				return true;
			else if(seg->first == rhs.seg->first)
				return cell->first < rhs.cell->first;
		}
		return false;
	}
};

typedef std::set< well_hit_cell > intersect_path;


/*-----------------------------------------------------------------
 * callback functor that triggers on intersecting bboxes
 *----------------------------------------------------------------*/
class intersect_action {
public:
	// ctor
	intersect_action(trimesh& mesh, well_path& wp, intersect_path& X)
		: m_(mesh), wp_(wp), x_(X)
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

	void cell_tri_cover(cell_data& d) {
		const t_float* V = d.V;
		d.cover.resize(12);

		// facet 0-1-2-3
		// 3angle 0-1-3
		d.cover[0] = Triangle_3(
			Point_3(V[X(0)], V[Y(0)], V[Z(0)]),
			Point_3(V[X(1)], V[Y(1)], V[Z(1)]),
			Point_3(V[X(3)], V[Y(3)], V[Z(3)])
		);
		// 3angle 0-2-3
		d.cover[1] = Triangle_3(
			Point_3(V[X(0)], V[Y(0)], V[Z(0)]),
			Point_3(V[X(2)], V[Y(2)], V[Z(2)]),
			Point_3(V[X(3)], V[Y(3)], V[Z(3)])
		);
		// facet 0-1-4-5
		// 3angle 0-1-5
		d.cover[2] = Triangle_3(
			Point_3(V[X(0)], V[Y(0)], V[Z(0)]),
			Point_3(V[X(1)], V[Y(1)], V[Z(1)]),
			Point_3(V[X(5)], V[Y(5)], V[Z(5)])
		);
		// 3angle 0-4-5
		d.cover[3] = Triangle_3(
			Point_3(V[X(0)], V[Y(0)], V[Z(0)]),
			Point_3(V[X(4)], V[Y(4)], V[Z(4)]),
			Point_3(V[X(5)], V[Y(5)], V[Z(5)])
		);
		// facet 4-5-6-7
		// 3angle 4-5-7
		d.cover[4] = Triangle_3(
			Point_3(V[X(4)], V[Y(4)], V[Z(4)]),
			Point_3(V[X(5)], V[Y(5)], V[Z(5)]),
			Point_3(V[X(7)], V[Y(7)], V[Z(7)])
		);
		// 3angle 4-6-7
		d.cover[5] = Triangle_3(
			Point_3(V[X(4)], V[Y(4)], V[Z(4)]),
			Point_3(V[X(6)], V[Y(6)], V[Z(6)]),
			Point_3(V[X(7)], V[Y(7)], V[Z(7)])
		);
		// facet 2-3-6-7
		// 3angle 2-3-7
		d.cover[6] = Triangle_3(
			Point_3(V[X(2)], V[Y(2)], V[Z(2)]),
			Point_3(V[X(3)], V[Y(3)], V[Z(3)]),
			Point_3(V[X(7)], V[Y(7)], V[Z(7)])
		);
		// 3angle 2-6-7
		d.cover[7] = Triangle_3(
			Point_3(V[X(2)], V[Y(2)], V[Z(2)]),
			Point_3(V[X(6)], V[Y(6)], V[Z(6)]),
			Point_3(V[X(7)], V[Y(7)], V[Z(7)])
		);
		// facet 0-2-4-6
		// 3angle 0-2-6
		d.cover[8] = Triangle_3(
			Point_3(V[X(0)], V[Y(0)], V[Z(0)]),
			Point_3(V[X(2)], V[Y(2)], V[Z(2)]),
			Point_3(V[X(6)], V[Y(6)], V[Z(6)])
		);
		// 3angle 0-4-6
		d.cover[9] = Triangle_3(
			Point_3(V[X(0)], V[Y(0)], V[Z(0)]),
			Point_3(V[X(4)], V[Y(4)], V[Z(4)]),
			Point_3(V[X(6)], V[Y(6)], V[Z(6)])
		);
		// facet 1-3-5-7
		// 3angle 1-3-7
		d.cover[10] = Triangle_3(
			Point_3(V[X(1)], V[Y(1)], V[Z(1)]),
			Point_3(V[X(3)], V[Y(3)], V[Z(3)]),
			Point_3(V[X(7)], V[Y(7)], V[Z(7)])
		);
		// 3angle 1-5-7
		d.cover[11] = Triangle_3(
			Point_3(V[X(1)], V[Y(1)], V[Z(1)]),
			Point_3(V[X(5)], V[Y(5)], V[Z(5)]),
			Point_3(V[X(7)], V[Y(7)], V[Z(7)])
		);
	}

	double calc_md(wp_iterator& fish, const Point_3& target) const {
		// walk all segments before the last one;
		double md = fish->second.md();
		//for(cwp_iterator p_seg = wp_.begin(), end = wp_.end(); p_seg != fish && p_seg != end; ++p_seg) {
		//	md += (p_seg->second).md();
		//}
		// append tail
		const t_float* W = fish->second.W;
		md += std::sqrt(Segment_3(Point_3(W[0], W[1], W[2]), target).squared_length());
		return md;
	}

	void operator()(const Box& bc, const Box& bw) {
		//std::cout << "intersection detected!" << std::endl;
		//return;

		trim_iterator cell_fish = static_cast< cell_box_handle* >(bc.handle().get())->data();
		wp_iterator well_fish = static_cast< well_box_handle* >(bw.handle().get())->data();
		cell_data& c = cell_fish->second;
		well_data& w = well_fish->second;

		// cover cell with triangles
		if(!c.cover.size())
			cell_tri_cover(c);

		// check that each triangle really intersects with given well segment
		Segment_3 s = w.segment();
		uint tri_count = 0;
		for(tri_iterator tri = c.cover.begin(), end = c.cover.end(); tri != end; ++tri, ++tri_count) {
			// test intersection
			if(!CGAL::do_intersect(s, *tri)) continue;
			// really do intersection
			Object xres = CGAL::intersection(s, *tri);
			// in 99% of cases we should get a point of intersection
			if(const Point_3* xpoint = CGAL::object_cast< Point_3 >(&xres))
				x_.insert(well_hit_cell(
					*xpoint, well_fish, cell_fish,
					calc_md(well_fish, *xpoint),
					tri_count >> 1
			));
			else if(const Segment_3* xseg = CGAL::object_cast< Segment_3 >(&xres)) {
				// in rare 1% of segment lying on the facet, add begin and end of segment
				x_.insert(well_hit_cell(
					xseg->source(), well_fish, cell_fish,
					calc_md(well_fish, xseg->source()),
					tri_count >> 1
				));
				x_.insert(well_hit_cell(
					xseg->target(), well_fish, cell_fish,
					calc_md(well_fish, xseg->target()),
					tri_count >> 1
				));
			}
		}

		//std::cout << "Box " << (a->handle() - triangles.begin()) << " and "
		//		<< (b->handle() - triangles.begin()) << " intersect";
		//if ( ! a->handle()->is_degenerate() && ! b->handle()->is_degenerate()
		//	&& CGAL::do_intersect( *(a->handle()), *(b->handle()))) {
		//	std::cout << ", and the triangles intersect also";
		//}
		//std::cout << '.' << std::endl;
	}

	void append_wp_nodes() {
		// walk through the intersection path and add node points
		// of well geometry to the cell with previous intersection
		intersect_path::iterator px = x_.begin();
		wp_iterator pw = wp_.begin();
		t_float node_md;
		t_float* W;
		for(wp_iterator end = wp_.end(); pw != end; ++pw) {
			node_md = pw->second.md();
			W = pw->second.W;
			while(px->md < node_md && px != x_.end())
				++px;
			// we need prev intersection
			if(px != x_.begin())
				--px;
			px = x_.insert(well_hit_cell(
				Point_3(W[0], W[1], W[2]),
				pw, px->cell, node_md,
				6
			)).first;
		}
		// well path doen't contain the end-point of trajectory
		// add it manually to the end of intersection
		--pw;
		px = x_.end(); --px;
		W = pw->second.W;
		x_.insert(well_hit_cell(
			Point_3(W[4], W[5], W[6]),
			pw, px->cell, W[7],
			6
		));
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
			*pr++ = px->where.z();
			// facet id
			*pr++ = px->facet;
		}

		return res;
	}

private:
	// mesh
	trimesh& m_;
	// well path
	well_path& wp_;
	// well-mesh intersection path
	intersect_path& x_;
};

// helper to create initial cell_data for each cell
spv_float coord_zcorn2trimesh(t_long nx, t_long ny, spv_float coord, spv_float zcorn, trimesh& res) {
	// build mesh_grdecl around given mesh
	sp_grd_mesh grd_src = BS_KERNEL.create_object(bs_mesh_grdecl::bs_type());
	grd_src->init_props(nx, ny, coord, zcorn);
	t_long nz = (zcorn->size() / nx / ny) >> 3;
	ulong n_cells = ulong(nx * ny * nz);

	// obtain coordinates for all vertices of all cells
	spv_float tops = grd_src->calc_cells_vertices();
	v_float::iterator pv = tops->begin();

	// fill trimesh with triangles corresponding to each cell
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
spv_float well_path_ident(t_long nx, t_long ny, spv_float coord, spv_float zcorn,
	spv_float well_info, bool include_well_nodes)
{
	// calculate mesh nodes coordinates and build initial trimesh
	trimesh M;
	spv_float tops;
	tops = coord_zcorn2trimesh(nx, ny, coord, zcorn, M);
	//t_long nz = (zcorn->size() / nx / ny) >> 3;

	// create bounding box for each cell in given mesh
	std::vector< Box > mesh_boxes(M.size());
	t_float lo[3];
	t_float hi[3];
	ulong cnt = 0;
	for(trim_iterator pm = M.begin(), end = M.end(); pm != end; ++pm, ++cnt) {
		// calc lo and hi for bounding box of cell
		const cell_data& d = pm->second;
		d.lo(lo); d.hi(hi);
		mesh_boxes[cnt] = Box(lo, hi, new cell_box_handle(pm));
	}

	// Create the corresponding vector of pointers to cells bounding boxes
	//std::vector< cell_box* > mesh_boxes_ptr;
	//for(std::vector< cell_box >::iterator i = mesh_boxes.begin(); i != mesh_boxes.end(); ++i)
	//	mesh_boxes_ptr.push_back(&*i);

	// create bounding boxes for line segments representing well trajectory
	ulong well_node_num = well_info->size() >> 2;
	if(well_node_num < 2) return spv_float();

	//vector< well_box > well_boxes;
	well_path W;
	vector< Box > well_boxes(well_node_num - 1);
	v_float::iterator pw = well_info->begin();
	for(ulong i = 0; i < well_node_num - 1; ++i) {
		//well_boxes.push_back(Box(
		//	lo, hi,
		//	new well_box_handle( W.insert(make_pair(i, pw)).first )
		//));

		well_data wd(pw);
		well_boxes[i] = Box(
			wd.bbox(),
			new well_box_handle(W.insert(make_pair(i, wd)).first)
		);

		pw += 4;
	}

	// Run the intersection algorithm with all defaults on the
	// indirect pointers to cell bounding boxes. Avoids copying the boxes
	intersect_path X;
	intersect_action A(M, W, X);
	CGAL::box_intersection_d(
		mesh_boxes.begin(), mesh_boxes.end(),
		well_boxes.begin(), well_boxes.end(),
		A
	);

	// finalize intersection
	if(include_well_nodes)
		A.append_wp_nodes();

	return A.export_1d();
}

/*-----------------------------------------------------------------
 * Python bindings
 *----------------------------------------------------------------*/
BOOST_PYTHON_FUNCTION_OVERLOADS(well_path_ident_overl, well_path_ident, 5, 6)
BOOST_PYTHON_FUNCTION_OVERLOADS(well_path_ident_overl_2d, well_path_ident_2d, 5, 6)

namespace python {

void py_export_wpi() {
	bp::def("well_path_ident", &well_path_ident, well_path_ident_overl());
	bp::def("well_path_ident_2d", &well_path_ident_2d, well_path_ident_overl_2d());
}

}

}	// eof blue-sky namespace

