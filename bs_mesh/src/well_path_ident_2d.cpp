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
typedef trimesh::iterator trim_iterator;

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
	// Alina wants z coord for node points - well, why not
	// Today she don't want Z - terminated
	//t_float z;

	well_hit_cell() {}
	well_hit_cell(const Point_2& where_, const wp_iterator& seg_,
		const trim_iterator& cell_, t_float md_, uint facet_, bool is_node_ = false)
		//double z_ = 0)
		: where(where_), seg(seg_), cell(cell_), md(md_), facet(facet_),
		is_node(is_node_)
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
		cell_data& c = cell_fish->second;
		well_data& w = well_fish->second;

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

	void append_wp_nodes() {
		if(!wp_.size()) return;
		// if intersections list is empty then we have a verticall well
		// that is located in single cell in X-Y plane
		// and we need to find which cell it hurts
		if(!x_.size()) {
			const well_data& wbegin = wp_.begin()->second;
			ulong cell_id = where_is_point(wbegin.start());
			if(cell_id >= m_.size()) return;
			x_.insert(well_hit_cell(
				wbegin.start(), wp_.begin(),
				m_.find(cell_id), 0, 4, true
			));
		}

		// walk through the intersection path and add node points
		// of well geometry to the cell with previous intersection
		intersect_path::iterator px = x_.begin();
		wp_iterator pw = wp_.begin();
		//t_float node_md;
		//t_float* W;
		for(wp_iterator end = wp_.end(); pw != end; ++pw) {
			const well_data& wseg = pw->second;
			//node_md = pw->second.md();
			//W = pw->second.W;
			while(px->md < wseg.md() && px != x_.end())
				++px;
			// we need prev intersection
			if(px != x_.begin())
				--px;
			px = x_.insert(well_hit_cell(
				wseg.start(),
				pw, px->cell, wseg.md(),
				4, true //, W[2]
			)).first;
		}

		// well path doesn't contain the end-point of trajectory
		// add it manually to the end of intersection
		--pw;
		const well_data& wend = pw->second;
		px = x_.end(); --px;
		//W = pw->second.W;
		//Point_2 wend(W[4], W[5]);
		x_.insert(well_hit_cell(
			wend.finish(), pw, px->cell,
			px->md + distance(px->where, wend.finish()),
			4, true //, W[6]
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
			//*pr++ = px->z;
			// facet id
			*pr++ = px->facet;
			// node flag
			*pr++ = px->is_node;
		}

		return res;
	}

	//Bbox_2 mesh_part_bbox(ulong x_start, ulong x_end, ulong y_start, ulong y_end) {
	//	// sanity checks
	//	x_start = min(x_start, nx_ - 1);
	//	x_end   = min(x_end  , nx_ - 1);
	//	y_start = min(y_start, ny_ - 1);
	//	y_end   = min(y_end  , ny_ - 1);

	//	// calc idx of "upper left" cell
	//	const ulong start_idx = y_start * nx_ + x_start;
	//	const ulong end_idx = y_end * nx_ + x_end;

	//	vertex_pos lo, hi;
	//	m_[start_idx].lo(lo);
	//	m_[end_idx].hi(hi);
	//	return Bbox_2(lo[0], lo[1], hi[0], hi[1]);
	//}

	ulong where_is_point(Point_2 point) {
		// start with full mesh
		// and divide it until we come to only one cell

		//mesh_part ceed(m_, nx_, ny_);
		typedef mesh_part::container_t parts_container;
		typedef mesh_part::container_t::iterator part_iterator;
		parts_container parts;
		parts.insert(mesh_part(m_, nx_, ny_));

		//ulong cell_id;
		while(parts.size()) {
			// store all peacies after split here
			parts_container leafs;
			for(part_iterator p = parts.begin(), end = parts.end(); p != end; ++p) {
				parts_container kids = p->divide();
				leafs.insert(kids.begin(), kids.end());
			}

			// now remove leafs that don't contain given point
			for(part_iterator l = leafs.begin(), end = leafs.end(); l != end; ) {
				if(l->bbox().has_on_unbounded_side(point)) {
					leafs.erase(l++);
					continue;
				}
				else if(l->count() == 1) {
					// if box contains only 1 cell - test if cell poly contains given point
					ulong cell_id = l->y_first * nx_ + l->x_first;
					if(!m_[cell_id].polygon().has_on_unbounded_side(point))
						return cell_id;
				}
				++l;
			}

			// leafs become the new start point for further division
			parts = leafs;
		}

		return -1;
	}

private:

	// x_last = last_element + 1 = x_size
	// y_last = last_element + 1 = y_size
	struct mesh_part {
		typedef set< mesh_part > container_t;

		mesh_part(trimesh& m, ulong nx, ulong ny)
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

		ulong size(int dim = 0) const {
			if(!dim)
				return x_last - x_first;
			else
				return y_last - y_first;
		}

		ulong count() const {
			return size() * size(1);
		}

		Iso_rectangle_2 bbox() const {
			// calc idx of "upper left" cell
			const ulong start_idx = y_first * nx_ + x_first;
			const ulong end_idx = (y_last - 1) * nx_ + x_last - 1;

			vertex_pos lo, hi;
			m_[start_idx].lo(lo);
			m_[end_idx].hi(hi);
			return Iso_rectangle_2(Point_2(lo[0], lo[1]), Point_2(hi[0], hi[1]));
		}

		container_t divide() const {
			// resulting split
			container_t res;
			insert_iterator< container_t > ii(res, res.begin());

			// split points
			ulong x_div[3], y_div[3];
			x_div[0] = x_first; x_div[2] = x_last;
			x_div[1] = x_first + (size() >> 1);
			y_div[0] = y_first; y_div[2] = y_last;
			y_div[1] = y_first + (size(1) >> 1);

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
		trimesh& m_;
		ulong nx_, ny_;

		mesh_part(trimesh& m, ulong nx, ulong ny,
				ulong x1, ulong x2,
				ulong y1, ulong y2)
			: x_first(x1), x_last(x2)
			, y_first(y1), y_last(y2)
			, m_(m), nx_(nx), ny_(ny)
		{}
	};

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
	// calculate mesh nodes coordinates and build initial trimesh
	trimesh M;
	spv_float tops;
	tops = coord_zcorn2trimesh(nx, ny, coord, zcorn, M);
	//t_long nz = (zcorn->size() / nx / ny) >> 3;

	// create bounding box for each cell in given mesh
	std::vector< Box > mesh_boxes(M.size());
	//vertex_pos lo, hi;
	ulong cnt = 0;
	for(trim_iterator pm = M.begin(), end = M.end(); pm != end; ++pm, ++cnt) {
		// calc lo and hi for bounding box of cell
		const cell_data& d = pm->second;
		//d.lo(lo); d.hi(hi);
		mesh_boxes[cnt] = Box(d.bbox(), new cell_box_handle(pm));
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
		// make bbox
		well_boxes[i] = Box(
			wd.bbox(),
			new well_box_handle(W.insert(make_pair(i, wd)).first)
		);

		pw += 4;
	}

	// Run the intersection algorithm with all defaults on the
	// indirect pointers to cell bounding boxes. Avoids copying the boxes
	intersect_path X;
	intersect_action A(M, W, X, nx, ny);
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

}	// eof blue-sky namespace


