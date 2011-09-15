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
#include <iterator>
#include <cmath>

#define MD_TOL 0.000001

namespace blue_sky { namespace wpi {

// boost::array with opertator= and ctor with elements init
// not needed right now
template< class T, int D >
class stat_array : public boost::array< T, D > {
public:
	typedef boost::array< T, D > base_t;
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
		std::copy(&data[0], &data[D], begin());
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
	stat_array& operator=(const stat_array< R, D >& rhs) {
		copy(rhs.begin(), rhs.begin() + min(size(), rhs.size()), begin());
		return *this;
	}
};

/*-----------------------------------------------------------------
 * implement well path identification algos depending on strategy
 *----------------------------------------------------------------*/
template< class strat_t >
struct wpi_algo {
	// common typedefs
	typedef t_ulong ulong;
	typedef t_uint uint;
	typedef v_float::iterator vf_iterator;
	typedef smart_ptr< bs_mesh_grdecl > sp_grd_mesh;

	// import strategy typedefs
	typedef typename strat_t::Object   Object;
	typedef typename strat_t::Kernel   Kernel;
	typedef typename strat_t::Point    Point;
	typedef typename strat_t::Segment  Segment;
	typedef typename strat_t::Bbox     Bbox;
	typedef typename strat_t::Iso_bbox Iso_bbox;

	typedef typename strat_t::vertex_pos   vertex_pos;
	typedef typename strat_t::vertex_pos_i vertex_pos_i;
	typedef typename strat_t::cell_pos     cell_pos;

	// import global consts
	// import global consts
	enum { D = strat_t::D, CVN = strat_t::CVN, inner_point_id = strat_t::inner_point_id };
	//using strat_t::D;
	//using strat_t::CVN;
	//using strat_t::inner_point_id;

	// import misc helper functions
	//typedef typename strat_t::decode_cell_id decode_cell_id;
	//typedef typename strat_t::encode_cell_id encode_cell_id;
	//typedef typename strat_t::vertex_pos2bbox vertex_pos2bbox;
	//typedef typename strat_t::vertex_pos2point vertex_pos2point;

	static void decode_cell_id(ulong id, vertex_pos_i& res, const vertex_pos_i& m_size) {
		strat_t::decode_cell_id(id, res, m_size);
	}
	static ulong encode_cell_id(const vertex_pos_i& p, const vertex_pos_i& m_size) {
		return strat_t::encode_cell_id(p, m_size);
	}

	static Bbox vertex_pos2bbox(const vertex_pos& lo, const vertex_pos& hi) {
		return strat_t::vertex_pos2bbox(lo, hi);
	}

	static Point vertex_pos2point(const vertex_pos& p) {
		return strat_t::vertex_pos2point(p);
	}

	// assign for c arrays
	// fun with returning reference to array :)
	template< class T >
	static T (&ca_assign(T (&lhs)[D], const T (&rhs)[D]))[D] {
		std::copy(&rhs[0], &rhs[D], &lhs[0]);
		return lhs;
	}

	template< class T >
	static T (&ca_assign(T (&lhs)[D], const T& v))[D] {
		std::fill(&lhs[0], &lhs[D], v);
		return lhs;
	}

	static Iso_bbox vertex_pos2rect(const vertex_pos& lo, const vertex_pos& hi) {
		return Iso_bbox(vertex_pos2point(lo), vertex_pos2point(hi));
	}

	/*-----------------------------------------------------------------
	* cell description
	*----------------------------------------------------------------*/
	struct cell_data_base {
		// vertex coord
		t_float* V;

		// empty ctor for map
		cell_data_base() : V(NULL) {}
		// std ctor
		cell_data_base(t_float *const cell) : V(cell) {}

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

		Bbox bbox() const {
			vertex_pos p1, p2;
			lo(p1); hi(p2);
			return vertex_pos2bbox(p1, p2);
		}

		Point ss(uint vert_idx) const {
			return vertex_pos2point(cpos()[vert_idx]);
		}
		Point operator[](uint vert_idx) const {
			return ss(vert_idx);
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

	typedef typename strat_t::template cell_data< cell_data_base > cell_data;
	typedef st_smart_ptr< cell_data > sp_cell_data;
	// storage for representing mesh
	typedef std::map< t_ulong, cell_data > trimesh;
	typedef typename trimesh::iterator trim_iterator;
	typedef typename trimesh::const_iterator ctrim_iterator;

	/*-----------------------------------------------------------------
	* represent rectangular part of mesh with splitting support
	*----------------------------------------------------------------*/
	// x_last = last_element + 1 = x_size
	// y_last = last_element + 1 = y_size
	struct mesh_part {
		typedef std::set< mesh_part > container_t;

		mesh_part(trimesh& m, const vertex_pos_i& mesh_size)
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
				lo[i] = std::min(lo[i], m_size_[i] - 1);
				hi[i] = std::min(hi[i], m_size_[i]);
				hi[i] = std::max(lo[i] + 1, hi[i]);
			}
		}

		void init(const std::vector< ulong >& cell_idx) {
			vertex_pos_i lower, upper, p;
			// search for bounds
			for(ulong i = 0; i < cell_idx.size(); ++i) {
				decode_cell_id(cell_idx[i], p, m_size_);
				if(i == 0)
					ca_assign(lower, ca_assign(upper, p));
				else {
					for(uint i = 0; i < D; ++i) {
						lower[i] = std::min(lower[i], p[i]);
						upper[i] = std::max(upper[i], p[i]);
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

		trim_iterator ss_iter(const vertex_pos_i& offset) {
			vertex_pos_i cell;
			ca_assign(cell, lo);
			for(uint i = 0; i < D; ++i)
				cell[i] += offset[i];
			return m_.find(encode_cell_id(cell, m_size_));
		}

		trim_iterator ss_iter(const ulong& offset) {
			// size of this part
			vertex_pos_i part_size;
			for(uint i = 0; i < D; ++i)
				part_size[i] = side_len(i);

			// plain id -> vertex_pos_i
			vertex_pos_i part_pos;
			decode_cell_id(offset, part_pos, part_size);
			// part_pos -> cell
			return ss_iter(part_pos);
		}

		Iso_bbox bbox() const {
			vertex_pos lo_pos, hi_pos;
			mesh_ss(lo).lo(lo_pos);
			// last = hi - 1
			vertex_pos_i last;
			ca_assign(last, hi);
			std::transform(&last[0], &last[D], &last[0], bind2nd(std::minus< ulong >(), 1));
			mesh_ss(last).hi(hi_pos);
			return vertex_pos2rect(lo_pos, hi_pos);
		}

		container_t divide() const {
			// resulting split
			container_t res;
			std::insert_iterator< container_t > ii(res, res.begin());

			// split points
			vertex_pos_i split_p[3];
			ca_assign(split_p[0], lo);
			ca_assign(split_p[2], hi);
			// middle
			for(uint i = 0; i < D; ++i)
				split_p[1][i] = lo[i] + (side_len(i) >> 1);

			// make splitting only if split containt more than 1 cell
			const ulong cube_num = 1 << D;
			vertex_pos_i spl_lo, spl_hi;
			ulong tot_sz = 0;
			for(ulong i = 0; i < cube_num; ++i) {
				ca_assign(spl_lo, split_p[0]);
				ca_assign(spl_hi, split_p[1]);

				ulong mask = 1, sz = 1;
				for(ulong j = 0; j < D; ++j, mask <<= 1) {
					if(i & mask) {
						spl_lo[j] = split_p[1][j];
						spl_hi[j] = split_p[2][j];
					}
					sz *= spl_hi[j] - spl_lo[j];
				}

				// add new child cell
				if(sz) {
					*ii++ = mesh_part(m_, m_size_, spl_lo, spl_hi);
					tot_sz += sz;
				}
			}

			assert(tot_sz == size());
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
		trimesh& m_;
		vertex_pos_i m_size_;

		mesh_part(trimesh& m, const vertex_pos_i& mesh_size,
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
	struct well_data_base {
		// segment begin, end and md in raw vector
		t_float* W;

		//empty ctor for map
		well_data_base() : W(NULL) {}
		//std ctor
		well_data_base(t_float *const segment) : W(segment) {}

		vertex_pos& cstart() {
			return reinterpret_cast< vertex_pos& >(*W);
		}

		const vertex_pos& cstart() const {
			return reinterpret_cast< const vertex_pos& >(*W);
		}

		vertex_pos& cfinish() {
			return reinterpret_cast< vertex_pos& >(*(W + 4));
		}

		const vertex_pos& cfinish() const {
			return reinterpret_cast< const vertex_pos& >(*(W + 4));
		}

		Point start() const {
			return vertex_pos2point(cstart());
		}
		Point finish() const {
			return vertex_pos2point(cfinish());
		}

		Segment segment() const {
			return Segment(start(), finish());
		}

		Bbox bbox() const {
			return segment().bbox();
		}

		double len() const {
			return std::sqrt(segment().squared_length());
		}
	};

	typedef typename strat_t::template well_data< well_data_base > well_data;
	typedef std::map< ulong, well_data > well_path;
	typedef typename well_path::iterator wp_iterator;
	typedef typename well_path::const_iterator cwp_iterator;

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
			typedef cell_data data_t;
		};
		// overload for well path
		template< class unused >
		struct fish2box_t< wp_iterator, unused > {
			enum { type = WELL_BOX };
			typedef well_data data_t;
		};
	};
	// pointer is really stored as box handle
	typedef st_smart_ptr< box_handle > sp_bhandle;

	template< class fish >
	class box_handle_impl : public box_handle {
	public:
		typedef fish fish_t;
		typedef typename box_handle::template fish2box_t< fish_t >::data_t data_t;

		box_handle_impl(const fish_t& f) : f_(f) {}

		int type() const {
			return box_handle::template fish2box_t< fish_t >::type;
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
	typedef CGAL::Box_intersection_d::Box_with_handle_d< double, D, sp_bhandle > Box;

	/*-----------------------------------------------------------------
	* intersections description
	*----------------------------------------------------------------*/
	struct well_hit_cell {
		// point of intersection
		Point where;
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
		well_hit_cell(const Point& where_, const wp_iterator& seg_,
			const trim_iterator& cell_, t_float md_, uint facet_, bool is_node_ = false)
			: where(where_), seg(seg_), cell(cell_), md(md_), facet(facet_), is_node(is_node_)
		{}

		// hit points ordered first by md
		bool operator <(const well_hit_cell& rhs) const {
			return md < rhs.md;
		}
	};

	// storage of intersection points
	typedef std::multiset< well_hit_cell > intersect_path;

	/*-----------------------------------------------------------------
	* callback functor that triggers on intersecting bboxes
	*----------------------------------------------------------------*/
	class intersect_action {
	public:
		typedef typename intersect_path::iterator x_iterator;

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
				decode_cell_id(r1->cell->first, c1, m_size_);
				decode_cell_id(r2->cell->first, c2, m_size_);

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
			const vertex_pos_i& m_size_;
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
			double md = fish->second.md();
			// append tail
			md += distance(fish->second.start(), target);
			return md;
		}

		void operator()(const Box& bc, const Box& bw) {
			typedef typename strat_t::xpoints_list xpoints_list;

			//trim_iterator cell_fish = static_cast< cell_box_handle* >(bc.handle().get())->data();
			//wp_iterator well_fish = static_cast< well_box_handle* >(bw.handle().get())->data();
			//cell_data& c = cell_fish->second;
			//well_data& w = well_fish->second;

			cell_box_handle* cell_h = static_cast< cell_box_handle* >(bc.handle().get());
			well_box_handle* well_h = static_cast< well_box_handle* >(bw.handle().get());

			xpoints_list res = strat_t::on_boxes_intersect(cell_h, well_h, bc, bw);

			// add all points to interscetion path
			for(typename xpoints_list::iterator px = res.begin(), end = res.end(); px != end; ++px) {
				x_.insert(well_hit_cell(
					px->first, well_h->data(), cell_h->data(),
					calc_md(well_h->data(), px->first),
					px->second
				));
			}

			//strat_t::on_boxes_intersect< wpi_impl >(bc, bw, x_);
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
				// lower_bound
				while(px != x_.end() && px->md < pw->second.md())
					++px;

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
				const well_data& seg = pw->second;
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
			res->resize(x_.size() * 7);
			vf_iterator pr = res->begin();

			for(typename intersect_path::const_iterator px = x_.begin(), end = x_.end(); px != end; ++px) {
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
				// node flag
				*pr++ = px->is_node;
			}

			return res;
		}

		std::vector< ulong > where_is_point(std::vector< Point > points) const {
			// start with full mesh
			// and divide it until we come to only one cell

			// mesh partition stored here
			typedef typename mesh_part::container_t parts_container;
			typedef typename mesh_part::container_t::iterator part_iterator;
			parts_container parts;
			parts.insert(mesh_part(m_, m_size_));

			// found cell_ids stored here
			std::vector< ulong > res(points.size(), -1);
			//ulong cell_id;
			while(parts.size()) {
				// split every part
				parts_container leafs;
				for(part_iterator p = parts.begin(), end = parts.end(); p != end; ++p) {
					parts_container kids = p->divide();
					leafs.insert(kids.begin(), kids.end());
				}

				// collection of points inside current partition
				std::list< ulong > catched_points;
				// process each leaf and find points inside it
				for(part_iterator l = leafs.begin(), end = leafs.end(); l != end; ) {
					const Iso_bbox& cur_rect = l->bbox();
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
						ulong cell_id = encode_cell_id(l->lo, m_size_);
						//Polygon_2 cell_poly = m_[cell_id].polygon();
						cell_data& cell = m_[cell_id];
						for(std::list< ulong >::iterator pp = catched_points.begin(),
							cp_end = catched_points.end();
							pp != cp_end; ++pp
							)
						{
							if(cell.contains(points[*pp]))
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

		ulong where_is_point(Point point) const {
			return where_is_point(std::vector< Point >(1, point))[0];
		}

	private:
		x_iterator insert_wp_node(ulong cell_id, wp_iterator pw, x_iterator px, bool end_point = false) {
			// initialization
			const well_data& wseg = pw->second;
			Point where = wseg.start();
			t_float wp_md = wseg.md();
			if(end_point) {
				where = wseg.finish();
				wp_md += wseg.len();
				//wp_md = px->md + distance(px->where, where);
			}

			// check if current or prev intersection match with node
			uint facet_id = inner_point_id;
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
			if(facet_id == inner_point_id)
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
		// mesh size
		vertex_pos_i m_size_;
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
		//typedef wpi_impl< strat_t > impl;
		//typedef typename impl::Point Point;
		//typedef typename impl::Box Box;

		//typedef typename impl::cell_data cell_data;
		//typedef typename impl::well_data well_data;
		//typedef typename impl::mesh_part mesh_part;
		//typedef typename impl::well_path well_path;
		//typedef typename impl::trimesh trimesh;
		//typedef typename impl::intersect_path intersect_path;

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
		for(ulong i = 0; i < well_node_num - 1; ++i) {
			well_data wd(pw);
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
		vertex_pos_i mesh_size = {nx, ny, nz};
		intersect_action A(M, W, X, mesh_size);
		const std::vector< ulong >& hit_idx = A.where_is_point(wnodes);
		// create part of mesh to process based on these cells
		mesh_part hot_mesh(M, mesh_size);
		hot_mesh.init(hit_idx);

		// create bounding box for each cell in given mesh
		std::vector< Box > mesh_boxes(hot_mesh.size());
		ulong cnt = 0;
		trim_iterator pm;
		for(ulong i = 0; i < hot_mesh.size(); ++i) {
			pm = hot_mesh.ss_iter(i);
			const cell_data& d = pm->second;
			mesh_boxes[cnt++] = Box(d.bbox(), new cell_box_handle(pm));
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

