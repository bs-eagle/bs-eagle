/// @file well_path_ident.cpp
/// @brief Well path and mesh intersection utilities
/// @author uentity
/// @version 0.1
/// @date 05.07.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

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

#include "bs_mesh_stdafx.h"
#include "conf.h"
#include "export_python_wrapper.h"
#include "bs_mesh_grdecl.h"

#include <vector>
#include <cmath>
// DEBUG
#include <iostream>

#define X(n) (3*n)
#define Y(n) (3*n + 1)
#define Z(n) (3*n + 2)
#define C(n, offs) (3*n + offs)

namespace bp = boost::python;
using namespace std;

namespace blue_sky {

typedef CGAL::Object                                        Object;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                                     Point_3;
typedef Kernel::Triangle_3                                  Triangle_3;
typedef Kernel::Segment_3                                   Segment_3;
typedef CGAL::Bbox_3                                        Bbox_3;
typedef std::vector<Triangle_3>                             Triangles;
typedef Triangles::iterator                                 tri_iterator;

typedef smart_ptr< bs_mesh_grdecl > sp_grd_mesh;
typedef t_ulong ulong;
typedef t_uint uint;
typedef v_float::iterator vf_iterator;

/*-----------------------------------------------------------------
 * cell description
 *----------------------------------------------------------------*/
typedef t_float vertex_pos[3];
typedef t_float cell_pos[8][3];

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

private:
	template< template< class > class pred >
	void bound(vertex_pos& b) const {
		pred< t_float > p = pred< t_float >();
		const cell_pos& cV = cpos();
		for(uint i = 0; i < 3; ++i) {
			t_float c = cV[0][i];
			for(uint j = 1; j < 8; ++j) {
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
typedef typename trimesh::iterator trim_iterator;

// well path description
/*-----------------------------------------------------------------
 * well_description
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
typedef typename well_path::iterator wp_iterator;
typedef typename well_path::const_iterator cwp_iterator;

/*-----------------------------------------------------------------
 * Bounding box description
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

	well_hit_cell() {}
	well_hit_cell(const Point_3& where_, const wp_iterator& seg_,
		const trim_iterator& cell_, t_float md_)
		: where(where_), seg(seg_), cell(cell_), md(md_)
	{}

	// hit points ordered first by well segment
	// next by md
	// and at last by cell number
	bool operator <(const well_hit_cell& rhs) const {
		if(seg->first < rhs.seg->first)
			return true;
		else if(seg->first == rhs.seg->first) {
			if(md < rhs.md)
				return true;
			else if(md == rhs.md)
				return cell->first < rhs.cell->first;
		}
		return false;
	}
};

typedef std::set< well_hit_cell > intersect_path;


// callback functor that triggers on intersecting bboxes
class intersect_action {
public:
	// ctor
	intersect_action(trimesh& mesh, well_path& wp, intersect_path& X)
		: m_(mesh), wp_(wp), x_(X)
	{}

	void cell_tri_cover(cell_data& d) {
		const t_float* V = d.V;
		// facet 0-1-2-3
		// 3angle 0-1-3
		d.cover.push_back(Triangle_3(
			Point_3(V[X(0)], V[Y(0)], V[Z(0)]),
			Point_3(V[X(1)], V[Y(1)], V[Z(1)]),
			Point_3(V[X(3)], V[Y(3)], V[Z(3)])
		));
		// 3angle 0-2-3
		d.cover.push_back(Triangle_3(
			Point_3(V[X(0)], V[Y(0)], V[Z(0)]),
			Point_3(V[X(2)], V[Y(2)], V[Z(2)]),
			Point_3(V[X(3)], V[Y(3)], V[Z(3)])
		));
		// facet 0-1-4-5
		// 3angle 0-1-5
		d.cover.push_back(Triangle_3(
			Point_3(V[X(0)], V[Y(0)], V[Z(0)]),
			Point_3(V[X(1)], V[Y(1)], V[Z(1)]),
			Point_3(V[X(5)], V[Y(5)], V[Z(5)])
		));
		// 3angle 0-4-5
		d.cover.push_back(Triangle_3(
			Point_3(V[X(0)], V[Y(0)], V[Z(0)]),
			Point_3(V[X(4)], V[Y(4)], V[Z(4)]),
			Point_3(V[X(5)], V[Y(5)], V[Z(5)])
		));
		// facet 4-5-6-7
		// 3angle 4-5-7
		d.cover.push_back(Triangle_3(
			Point_3(V[X(4)], V[Y(4)], V[Z(4)]),
			Point_3(V[X(5)], V[Y(5)], V[Z(5)]),
			Point_3(V[X(7)], V[Y(7)], V[Z(7)])
		));
		// 3angle 4-6-7
		d.cover.push_back(Triangle_3(
			Point_3(V[X(4)], V[Y(4)], V[Z(4)]),
			Point_3(V[X(6)], V[Y(6)], V[Z(6)]),
			Point_3(V[X(7)], V[Y(7)], V[Z(7)])
		));
		// facet 2-3-6-7
		// 3angle 2-3-7
		d.cover.push_back(Triangle_3(
			Point_3(V[X(2)], V[Y(2)], V[Z(2)]),
			Point_3(V[X(3)], V[Y(3)], V[Z(3)]),
			Point_3(V[X(7)], V[Y(7)], V[Z(7)])
		));
		// 3angle 2-6-7
		d.cover.push_back(Triangle_3(
			Point_3(V[X(2)], V[Y(2)], V[Z(2)]),
			Point_3(V[X(6)], V[Y(6)], V[Z(6)]),
			Point_3(V[X(7)], V[Y(7)], V[Z(7)])
		));
		// facet 0-2-4-6
		// 3angle 0-2-6
		d.cover.push_back(Triangle_3(
			Point_3(V[X(0)], V[Y(0)], V[Z(0)]),
			Point_3(V[X(2)], V[Y(2)], V[Z(2)]),
			Point_3(V[X(6)], V[Y(6)], V[Z(6)])
		));
		// 3angle 0-4-6
		d.cover.push_back(Triangle_3(
			Point_3(V[X(0)], V[Y(0)], V[Z(0)]),
			Point_3(V[X(4)], V[Y(4)], V[Z(4)]),
			Point_3(V[X(6)], V[Y(6)], V[Z(6)])
		));
		// facet 1-3-5-7
		// 3angle 1-3-7
		d.cover.push_back(Triangle_3(
			Point_3(V[X(1)], V[Y(1)], V[Z(1)]),
			Point_3(V[X(3)], V[Y(3)], V[Z(3)]),
			Point_3(V[X(7)], V[Y(7)], V[Z(7)])
		));
		// 3angle 1-5-7
		d.cover.push_back(Triangle_3(
			Point_3(V[X(1)], V[Y(1)], V[Z(1)]),
			Point_3(V[X(5)], V[Y(5)], V[Z(5)]),
			Point_3(V[X(7)], V[Y(7)], V[Z(7)])
		));
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
		for(tri_iterator tri = c.cover.begin(), end = c.cover.end(); tri != end; ++tri) {
			// test intersection
			if(!CGAL::do_intersect(s, *tri)) continue;
			// really do intersection
			Object xres = CGAL::intersection(s, *tri);
			// in 99% of cases we should get a point of intersection
			if(const Point_3* xpoint = CGAL::object_cast< Point_3 >(&xres))
				x_.insert(well_hit_cell(
					*xpoint, well_fish, cell_fish,
					calc_md(well_fish, *xpoint)
			));
			else if(const Segment_3* xseg = CGAL::object_cast< Segment_3 >(&xres)) {
				// in rare 1% of segment lying on the facet, add begin and end of segment
				x_.insert(well_hit_cell(
					xseg->source(), well_fish, cell_fish,
					calc_md(well_fish, xseg->source())
				));
				x_.insert(well_hit_cell(
					xseg->target(), well_fish, cell_fish,
					calc_md(well_fish, xseg->target())
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
				pw, px->cell, node_md
			)).first;
		}
		// well path doen't contain the end-point of trajectory
		// add it manually to the end of intersection
		--pw;
		px = x_.end(); --px;
		W = pw->second.W;
		x_.insert(well_hit_cell(
			Point_3(W[4], W[5], W[6]),
			pw, px->cell, W[7]
		));
	}

	spv_float export_1d() const {
		spv_float res = BS_KERNEL.create_object(v_float::bs_type());
		res->resize(x_.size() * 5);
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

spv_float well_path_ident(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_float well_info) {
	// calculate mesh nodes coordinates and muild initial trimesh
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
		// +1 because we need to skip md field
		//lo[0] = min(pw[X(0)], pw[X(1) + 1]);
		//lo[1] = min(pw[Y(0)], pw[Y(1) + 1]);
		//lo[2] = min(pw[Z(0)], pw[Z(1) + 1]);
		//hi[0] = max(pw[X(0)], pw[X(1) + 1]);
		//hi[1] = max(pw[Y(0)], pw[Y(1) + 1]);
		//hi[2] = max(pw[Z(0)], pw[Z(1) + 1]);
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
	A.append_wp_nodes();

	return A.export_1d();
}

namespace python {

void py_export_wpi() {
	bp::def("well_path_ident", &well_path_ident);
}

}

}	// eof blue-sky namespace

