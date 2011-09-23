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

namespace blue_sky { namespace wpi {

template< class strat_t >
struct wpi_algo_helpers {
	// import strategy typedefs
	typedef typename strat_t::Point    Point;
	typedef typename strat_t::Segment  Segment;
	typedef typename strat_t::Bbox     Bbox;
	typedef typename strat_t::Iso_bbox Iso_bbox;

	typedef typename strat_t::vertex_pos   vertex_pos;
	typedef typename strat_t::vertex_pos_i vertex_pos_i;

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
};

template< class strat_t >
struct wpi_algo_pod : public wpi_algo_helpers< strat_t > {
	// import strategy typedefs
	typedef typename strat_t::Point    Point;
	typedef typename strat_t::Segment  Segment;
	typedef typename strat_t::Bbox     Bbox;
	typedef typename strat_t::Iso_bbox Iso_bbox;

	typedef typename strat_t::vertex_pos   vertex_pos;
	typedef typename strat_t::vertex_pos_i vertex_pos_i;

	// import global consts
	enum { D = strat_t::D };
	//typedef wpi_algo_helpers< strat_t > helper_t;

	/*-----------------------------------------------------------------
	* cell description
	*----------------------------------------------------------------*/
	struct cell_data_base {
		// actual vertex number
		enum { N = (1 << D) };
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
			// actual vertex number to search
			//const uint N = (1 << D);
			for(uint i = 0; i < D; ++i) {
				t_float c = cV[0][i];
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
	//typedef std::map< t_ulong, cell_data > trimesh;
	typedef std::vector< cell_data > trimesh;
	typedef typename trimesh::iterator trim_iterator;
	typedef typename trimesh::const_iterator ctrim_iterator;

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
	* intersections description
	*----------------------------------------------------------------*/
	struct well_hit_cell {
		// point of intersection
		Point where;
		// what segment of well
		wp_iterator seg;
		// interseect with which cell
		//trim_iterator cell;
		ulong cell;
		// depth along well in point of intersection
		t_float md;
		// cell facet
		uint facet;
		// is that point a node?
		bool is_node;

		well_hit_cell() {}
		well_hit_cell(const Point& where_, const wp_iterator& seg_,
			const ulong& cell_, t_float md_, uint facet_, bool is_node_ = false)
			: where(where_), seg(seg_), cell(cell_), md(md_), facet(facet_), is_node(is_node_)
		{}

		// hit points ordered first by md
		bool operator <(const well_hit_cell& rhs) const {
			return md < rhs.md;
		}
	};

	// storage of intersection points
	typedef std::multiset< well_hit_cell > intersect_path;
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

}} /* blue_sky::namespace wpi */

#endif /* end of include guard: WPI_ALGO_POD_BDBOLFWA */

