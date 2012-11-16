/// @file wpi_algo_pod.h
/// @brief Basic data structures needed to implement WPI algorithms
/// @author uentity
/// @version 
/// @date 19.09.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_ALGO_POD_BDBOLFWA
#define WPI_ALGO_POD_BDBOLFWA

#include "wpi_common.h"
//#include <boost/pool/pool_alloc.hpp>

namespace blue_sky { namespace wpi {

template< class strat_t >
struct helpers {
	// import strategy typedefs
	typedef typename strat_t::Point    Point;
	typedef typename strat_t::Segment  Segment;
	typedef typename strat_t::Bbox     Bbox;
	typedef typename strat_t::Iso_bbox Iso_bbox;

	typedef typename strat_t::vertex_pos   vertex_pos;
	typedef typename strat_t::vertex_pos_i vertex_pos_i;

	typedef typename strat_t::traits_t strat_traits;
	typedef typename strat_traits::cell_vertex_iterator cell_vertex_iterator;

	// import global consts
	enum { D = strat_t::D };

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

	static Point rawptr2point(const t_float* p) {
		//return strat_t::vertex_pos2point(strat_traits::template iter2pos< const vertex_pos >(p));
		return strat_t::vertex_pos2point(reinterpret_cast< const vertex_pos& >(*p));
	}

	// cell_vertex_iterator -> point
	static Point iter2point(const cell_vertex_iterator& p) {
		return strat_t::vertex_pos2point(strat_traits::template iter2pos< const vertex_pos >(p));
	}

	static Iso_bbox vertex_pos2rect(const vertex_pos& lo, const vertex_pos& hi) {
		return Iso_bbox(vertex_pos2point(lo), vertex_pos2point(hi));
	}
};

template< class strat_t >
struct pods : public helpers< strat_t > {
	// import strategy typedefs
	typedef typename strat_t::Point    Point;
	typedef typename strat_t::Segment  Segment;
	typedef typename strat_t::Bbox     Bbox;
	typedef typename strat_t::Iso_bbox Iso_bbox;

	typedef typename strat_t::vertex_pos   vertex_pos;
	typedef typename strat_t::vertex_pos_i vertex_pos_i;

	typedef typename strat_t::cell_vertex_iterator cell_vertex_iterator;
	typedef typename strat_t::well_traj_iterator   well_traj_iterator;

	typedef typename strat_t::traits_t strat_traits;

	// import base functions
	typedef helpers< strat_t > base_t;
	using base_t::decode_cell_id;
	using base_t::encode_cell_id;
	using base_t::vertex_pos2bbox;
	using base_t::vertex_pos2point;
	using base_t::rawptr2point;
	using base_t::vertex_pos2rect;

	// import global consts
	enum { D = strat_t::D };
	//typedef helpers< strat_t > helper_t;

	/*-----------------------------------------------------------------
	* cell description
	*----------------------------------------------------------------*/
	struct cell_data_base {
		// actual vertex number
		enum { N = (1 << D) };
		// vertex coord
		//t_float* V;
		cell_vertex_iterator V;

		// empty ctor for map
		cell_data_base() : V(NULL) {}
		// std ctor
		cell_data_base(const cell_vertex_iterator& cell) : V(cell) {}

		void lo(vertex_pos& b) const {
			bound< std::less >(b);
		}

		void hi(vertex_pos& b) const {
			bound< std::greater >(b);
		}

		cell_pos& cpos() {
			return strat_traits::template iter2pos< cell_pos& >(V);
			//return reinterpret_cast< cell_pos& >(*V);
		}

		const cell_pos& cpos() const {
			return strat_traits::template iter2pos< const cell_pos& >(V);
			//return reinterpret_cast< const cell_pos& >(*V);
		}

		Bbox bbox() const {
			vertex_pos p1, p2;
			lo(p1); hi(p2);
			return vertex_pos2bbox(p1, p2);
		}

		Iso_bbox iso_bbox() const {
			vertex_pos p1, p2;
			lo(p1); hi(p2);
			return vertex_pos2rect(p1, p2);
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
			t_float c;
			for(uint i = 0; i < D; ++i) {
				c = cV[0][i];
				for(uint j = 1; j < N; ++j) {
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
	//template< class cell_data, class strat_traits >
	class BS_API_PLUGIN trimesh {
	public:
		typedef cell_data value_type;
		typedef value_type& reference;
		typedef const value_type& const_reference;
		typedef value_type* pointer;

		// empty ctor
		trimesh();
		// ctor from given COORD & ZCORN
		trimesh(t_long nx, t_long ny, spv_float coord, spv_float zcorn);

		void init(t_long nx, t_long ny, spv_float coord, spv_float zcorn);

		const vertex_pos_i& size() const {
			return size_;
		}

		ulong size_flat() const {
			ulong sz = 1;
			for(uint i = 0; i < D; ++i)
				sz *= size_[i];
			return sz;
		}

		// direct subscript via backend ignoring cache
		// calls pimpl actually
		value_type ss_backend(ulong idx) const;

		// subscripting operator - returns a new _copy_ of cell every time!
		// if cell was already cached, return a copy from cache
		value_type ss(ulong idx) const {
			//typename cache_t::const_iterator r = cache_.find(idx);
			//if(r != cache_.end())
			//	return r->second;
			return ss_backend(idx);
		}

		// the same in operator form
		value_type operator[](ulong idx) const {
			return ss(idx);
		}

		// if cell was modified, we can return cached version instead of subscripting backend
		//void cache_cell(ulong idx, const value_type& cell) {
		//	cache_[idx] = cell;
		//}

		// obtain iterators on backend
		cell_vertex_iterator begin() const;
		cell_vertex_iterator end() const;

	private:
		struct impl;
		st_smart_ptr< impl > pimpl_;

		vertex_pos_i size_;
		//typedef std::map< ulong, value_type > cache_t;
		//cache_t cache_;
	};

	//typedef std::vector< cell_data > trimesh;
	//typedef typename trimesh::iterator trim_iterator;
	//typedef typename trimesh::const_iterator ctrim_iterator;

	/*-----------------------------------------------------------------
	* well description
	*----------------------------------------------------------------*/
	struct well_data_base {
		// segment begin, end and md in raw vector
		well_traj_iterator W;

		//empty ctor for map
		well_data_base() : W(NULL) {}
		//std ctor
		well_data_base(const well_traj_iterator& segment) : W(segment) {}

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

		Iso_bbox iso_bbox() const {
			return Iso_bbox(start(), finish());
		}

		double len() const {
			return std::sqrt(segment().squared_length());
		}
	};

	typedef typename strat_t::template well_data< well_data_base > well_data;
	typedef std::vector< well_data > well_path;
	typedef typename well_path::iterator wp_iterator;
	typedef typename well_path::const_iterator cwp_iterator;

	/*-----------------------------------------------------------------
	* intersections description
	*----------------------------------------------------------------*/
	struct well_hit_cell {
		// point of intersection
		Point where;
		// what segment of well
		ulong seg;
		// intersect with which cell
		ulong cell;
		// depth along well in point of intersection
		t_float md;
		// cell facet that contains intersection
		uint facet;
		// is that point a node?
		bool is_node;

		well_hit_cell() {}
		well_hit_cell(const Point& where_, ulong seg_,
			ulong cell_, t_float md_, uint facet_, bool is_node_ = false)
			: where(where_), seg(seg_), cell(cell_), md(md_), facet(facet_), is_node(is_node_)
		{}

		// ctor for searching
		well_hit_cell(t_float md_)
			: seg(0), cell(0), md(md_), facet(0), is_node(false)
		{}

		// hit points ordered first by md
		bool operator <(const well_hit_cell& rhs) const {
			return md < rhs.md;
		}
	};

	// storage of intersection points
	typedef std::multiset< well_hit_cell > intersect_path;

	//typedef boost::fast_pool_allocator<
	//	well_hit_cell,
	//	boost::default_user_allocator_new_delete,
	//	boost::details::pool::null_mutex
	//	> whc_allocator;
	//typedef std::multiset< well_hit_cell, std::less< well_hit_cell >, whc_allocator > intersect_path;
};

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
		std::fill(begin(), end(), value_type());
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
		for(int i = 0; i < N; ++i) {
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

}} /* blue_sky::namespace wpi */

#endif /* end of include guard: WPI_ALGO_POD_BDBOLFWA */

