/// @file wpi_algo_meshp.h
/// @brief Mesh part representation and dividing
/// @author uentity
/// @version 
/// @date 19.09.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_ALGO_MESHP_RIC3ZQNS
#define WPI_ALGO_MESHP_RIC3ZQNS

#include "wpi_algo_pod.h"
#include <algorithm>

namespace blue_sky { namespace wpi {

template< class strat_t >
struct mesh_tools : public helpers< strat_t > {
	// basic types
	typedef typename strat_t::vertex_pos   vertex_pos;
	typedef typename strat_t::vertex_pos_i vertex_pos_i;

	typedef typename strat_t::Point    Point;
	typedef typename strat_t::Iso_bbox Iso_bbox;
	typedef typename strat_t::Bbox     Bbox;

	// import global consts
	enum { D = strat_t::D };

	// import pods
	typedef pods< strat_t > pods_t;
	typedef typename pods_t::cell_data cell_data;
	typedef typename pods_t::well_data well_data;
	typedef typename pods_t::trimesh trimesh;
	//typedef typename pods_t::trim_iterator trim_iterator;
	//typedef typename pods_t::ctrim_iterator ctrim_iterator;

	// import helper functions
	typedef helpers< strat_t > base_t;
	using base_t::encode_cell_id;
	using base_t::decode_cell_id;
	using base_t::vertex_pos2bbox;
	using base_t::vertex_pos2rect;

	/*-----------------------------------------------------------------
	* represent rectangular part of mesh with splitting support
	*----------------------------------------------------------------*/
	// x_last = last_element + 1 = x_size
	// y_last = last_element + 1 = y_size
	struct mesh_part : public helpers< strat_t > {
		// CGAL intersection algorithm can find that the same mesh part
		// intersects with different well segments many times,
		// but we need to process only unique mesh_parts,
		// so use std::set instead of std::list as mesh_parts container
		typedef std::set< mesh_part > container_t;

		enum { n_facets = cell_data::n_facets };
		enum { n_edges = cell_data::n_edges };
		typedef ulong cell_neighb_enum[n_facets];
		typedef ulong edge_neighb_enum[n_edges];

		// ctor 1 - mesh part coincide with full mesh
		mesh_part(const trimesh& m, bool dont_calc_bbox = false)
			: m_(m)
		{
			ca_assign(lo, ulong(0));
			ca_assign(hi, m.size());
			//ca_assign(m_size_, mesh_size);
			calc_bounds(dont_calc_bbox);
		}

		void init(const vertex_pos_i& lower, const vertex_pos_i& upper, bool dont_calc_bbox = false) {
			ca_assign(lo, lower);
			ca_assign(hi, upper);

			// sanity checks
			for(uint i = 0; i < D; ++i) {
				lo[i] = std::min(lo[i], m_.size()[i] - 1);
				hi[i] = std::min(hi[i], m_.size()[i]);
				hi[i] = std::max(lo[i] + 1, hi[i]);
			}

			// precalc bounds
			calc_bounds(dont_calc_bbox);
		}

		void init(const std::vector< ulong >& cell_idx, bool dont_calc_bbox = false) {
			vertex_pos_i lower, upper, p;
			// search for bounds
			for(ulong i = 0; i < cell_idx.size(); ++i) {
				decode_cell_id(cell_idx[i], p, m_.size());
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
			init(lower, upper, dont_calc_bbox);
		}

		void init(ulong lower, ulong upper, bool dont_calc_bbox = false) {
			std::vector< ulong > idx(2);
			idx[0] = lower; idx[1] = upper;
			init(idx, dont_calc_bbox);
		}

		ulong side_len(uint dim) const {
			if(dim < D)
				return hi[dim] - lo[dim];
			else
				return hi[D - 1] - lo[D - 1];
		}

		// flat size of this mesh part
		ulong size() const {
			return sz_flat_;
		}

		ulong ss_id(const vertex_pos_i& offset) const {
			vertex_pos_i cell;
			ca_assign(cell, lo);
			for(uint i = 0; i < D; ++i)
				cell[i] += offset[i];
			return encode_cell_id(cell, m_.size());
		}

		ulong ss_id(ulong offset) const {
			// plain id -> vertex_pos_i
			vertex_pos_i part_pos;
			decode_cell_id(offset, part_pos, mp_size_);
			// part_pos -> cell
			return ss_id(part_pos);
		}

		// access to individual cells
		template< class index_t >
		cell_data operator[](index_t idx) const {
			return m_[ss_id(idx)];
		}

		//template< class index_t >
		//const cell_data& operator[](index_t idx) const {
		//	return m_[ss_id(idx)];
		//}

		Iso_bbox iso_bbox() const {
			//vertex_pos lo_pos, hi_pos;
			//bounds(lo_pos, hi_pos);
			return vertex_pos2rect(lo_bbox_, hi_bbox_);
		}

		Bbox bbox() const {
			//vertex_pos lo_pos, hi_pos;
			//bounds(lo_pos, hi_pos);
			return vertex_pos2bbox(lo_bbox_, hi_bbox_);
		}

		container_t divide(bool dont_calc_bbox = false) const {
			// resulting split
			container_t res;
			std::insert_iterator< container_t > ii(res, res.begin());

			// split points
			vertex_pos_i split_p[3];
			ca_assign(split_p[0], lo);
			ca_assign(split_p[2], hi);
			// middle
			for(uint i = 0; i < D; ++i)
				split_p[1][i] = lo[i] + (mp_size_[i] >> 1);

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
					*ii++ = mesh_part(m_, spl_lo, spl_hi, dont_calc_bbox);
					tot_sz += sz;
				}
			}

			assert(tot_sz == size());
			return res;
		}

		// for sorted containers
		bool operator <(const mesh_part& rhs) const {
			// per-element lexicographical compare
			const int r = lexicographical_compare_3way(
				&lo[0], &lo[D], &rhs.lo[0], &rhs.lo[D]
			);
			if(r == 0)
				return std::lexicographical_compare(
					&hi[0], &hi[D], &rhs.hi[0], &rhs.hi[D]
				);
			else
				return (r < 0);

			// if same objects
			//if(&lo[0] == &rhs.lo[0] || &hi[0] == &rhs.hi[0])
			//	return false;
			//for(uint i = 0; i < D; ++i) {
			//	if(lo[i] < rhs.lo[i])
			//		return true;
			//	else if(lo[i] > rhs.lo[i])
			//		return false;
			//	else if(hi[i] < rhs.hi[i])
			//		return true;
			//	else if(hi[i] > rhs.hi[i])
			//		return false;
			//}
			//return false;
		}

		// assignment operator for containers
		// assign only mesh_parts that belong to the same mesh!
		mesh_part& operator=(const mesh_part& rhs) {
			ca_assign(lo, rhs.lo);
			ca_assign(hi, rhs.hi);
			ca_assign(mp_size_, rhs.mp_size_);
			ca_assign(lo_bbox_, rhs.lo_bbox_);
			ca_assign(hi_bbox_, rhs.hi_bbox_);
			return *this;
		}

		// calc size of cell in x-y-z directions
		void cell_size(ulong offset, vertex_pos& res) const {
			if(offset < this->size()) {
				const cell_data& c = (*this)[offset];
				vertex_pos b1, b2;
				c.lo(b1); c.hi(b2);
				std::transform(&b2[0], &b2[D], &b1[0], &res[0], std::minus< t_float >());
			}
			//ctrim_iterator pc = this->ss_iter(offset);
			//if(pc != m_.end()) {
			//	vertex_pos b1, b2;
			//	pc->lo(b1); pc->hi(b2);
			//	std::transform(&b2[0], &b2[D], &b1[0], &res[0], std::minus< t_float >());
			//}
		}

		void cell_size(const vertex_pos_i& offset, vertex_pos& res) {
			cell_size(ss_id(offset), res);
		}

		// returned indexes corresponds to offsets inside this mesh_part!
		// for full mesh correspond to simple index
		void cell_neighbours(ulong idx, cell_neighb_enum& res) const {
			// init resulting array with invalid neighbours
			//cell_neighb_enum res;
			std::fill(&res[0], &res[n_facets], ulong(-1));

			// sanity check
			if(!size() || idx >= size())
				return;

			// obtain D-dim cell id
			vertex_pos_i cell_id, nb;
			decode_cell_id(idx, cell_id, mp_size_);

			// vary different dims and check if cell is inside this mesh_part
			ca_assign(nb, cell_id);
			for(uint i = 0; i < D; ++i) {
				--nb[i];
				//if(nb[i] >= lo[i] && nb[i] < hi[i])
				if(nb[i] < mp_size_[i])
					res[cell_data::facet_id(i, 0)] = encode_cell_id(nb, mp_size_);
				nb[i] += 2;
				//if(nb[i] < hi[i] && nb[i] >= lo[i])
				if(nb[i] < mp_size_[i])
					res[cell_data::facet_id(i, 1)] = encode_cell_id(nb, mp_size_);
				// return nb[i] to initial state
				--nb[i];
			}
		}

		// check that mesh_part consists only of boundaries
		bool is_pure_boundary() const {
			// if any dimesnsion is of size <= 2 then boundary is entirely *this
			for(uint i = 0; i < D; ++i) {
				if(mp_size_[i] < 3)
					return true;
			}
			return false;
		}

		// return array of mesh_parts representing boundary of this mesh
		std::list< mesh_part > boundary() const {
			// if any dimesnsion is of size <= 2 then boundary is entirely *this
			std::list< mesh_part > res;
			if(is_pure_boundary()) {
				res.push_back(*this);
				return res;
			}

			typedef typename strat_t::bbox_bnd_offs bbox_bnd_offs;
			// we have two boundary planes for each dimension
			vertex_pos_i bnd_hi, bnd_lo;
			//uint bnd_idx = 0;
			for(uint i = 0; i < D; ++i) {
				// first boundary - all dims, besides i, from 0 to size
				// i-th dim from 0 to 1
				// bnd_lo = lo;
				ca_assign(bnd_lo, lo);
				ca_assign(bnd_hi, hi);
				bnd_hi[i] = lo[i] + 1;
				// correct dims of 1st boundary using strategy-specific offsets
				// to make them finally non-intersecting
				const bbox_bnd_offs& offs0 = strat_t::bbox_boundary_offs(i, 0);
				std::transform(&bnd_lo[0], &bnd_lo[D], &offs0[0][0], &bnd_lo[0], std::plus< long >());
				std::transform(&bnd_hi[0], &bnd_hi[D], &offs0[1][0], &bnd_hi[0], std::plus< long >());
				// save resulting boundary
				res.push_back(mesh_part(m_, bnd_lo, bnd_hi));

				// first boundary - all dims, besides i, from 0 to size
				// i-th dim from size - 1 to size
				ca_assign(bnd_lo, lo);
				bnd_lo[i] = hi[i] - 1;
				// bnd_hi = hi
				ca_assign(bnd_hi, hi);
				// correct dims of 1st boundary using strategy-specific offsets
				// to make them finally non-intersecting
				const bbox_bnd_offs& offs1 = strat_t::bbox_boundary_offs(i, 1);
				std::transform(&bnd_lo[0], &bnd_lo[D], &offs1[0][0], &bnd_lo[0], std::plus< long >());
				std::transform(&bnd_hi[0], &bnd_hi[D], &offs1[1][0], &bnd_hi[0], std::plus< long >());
				// save resulting boundary
				res.push_back(mesh_part(m_, bnd_lo, bnd_hi));
			}
			return res;
		}

		// access to underlying trimesh object
		const trimesh& backend() const {
			return m_;
		}

		// functions to convert global cell offset to local offset in this mesh_part
		void local_cid(const ulong global_cid, vertex_pos_i& res) const {
			decode_cell_id(global_cid, res, m_.size());
			// res -= lo
			std::transform(&res[0], &res[D], &lo[0], &res[0], std::minus< ulong >());
		}

		ulong local_cid(const ulong global_cid) const {
			vertex_pos_i loc_cid;
			local_cid(global_cid, loc_cid);
			return encode_cell_id(loc_cid, mp_size_);
		}

		// public members
		vertex_pos_i lo, hi;

	private:
		const trimesh& m_;
		vertex_pos_i mp_size_;
		vertex_pos lo_bbox_, hi_bbox_;
		ulong sz_flat_;

		// idx SHOULD BE IN MESH!
		//const cell_data& ss(ulong idx) const {
		//	return m_[idx];
		//}
		cell_data ss(ulong idx) const {
			return m_[idx];
		}

		//const cell_data& ss(const vertex_pos_i& idx) const {
		//	return m_[encode_cell_id(idx, m_.size())];
		//}
		cell_data ss(const vertex_pos_i& idx) const {
			return m_[encode_cell_id(idx, m_.size())];
		}

		mesh_part(const trimesh& m,
				const vertex_pos_i& first_,
				const vertex_pos_i& last_,
				bool dont_calc_bbox = false)
			: m_(m)
		{
			ca_assign(lo, first_);
			ca_assign(hi, last_);
			//ca_assign(m_size_, mesh_size);
			calc_bounds(dont_calc_bbox);
		}

		void calc_bounds(bool dont_calc_bbox) {
			// calc size of this mesh part
			// and check if it is bigger than conjunction of boundaries
			bool is_exact_boundary = false;
			sz_flat_ = 1;
			for(uint i = 0; i < D; ++i) {
				mp_size_[i] = side_len(i);
				if(mp_size_[i] < 3)
					is_exact_boundary = true;
				sz_flat_ *= mp_size_[i];
			}

			if(dont_calc_bbox) return;

			// if mesh consists only from boundaries, then calc bounds
			// by simply iterating over all cells
			if(is_exact_boundary) {
				calc_bbox_raw(lo_bbox_, hi_bbox_);
				return;
			}

			// otherwise extract boundary and calc bounds based on it
			typedef std::list< mesh_part > boundary_cont;
			typedef typename boundary_cont::const_iterator cb_iterator;
			const boundary_cont& B = boundary();
			//vertex_pos lo_part_bbox, hi_part_bbox;
			for(cb_iterator pb = B.begin(), end = B.end(); pb != end; ++pb) {
				//pb->calc_bbox_raw(lo_part_bbox, hi_part_bbox);
				if(pb == B.begin()) {
					ca_assign(lo_bbox_, pb->lo_bbox_);
					ca_assign(hi_bbox_, pb->hi_bbox_);
				}
				else {
					typedef std::pointer_to_binary_function< const t_float&, const t_float&, const t_float& > bin_op_ptr;
					// lo_bbox_ = min(lo_bbox_, lo_part_bbox)
					std::transform(&lo_bbox_[0], &lo_bbox_[D], &pb->lo_bbox_[0], &lo_bbox_[0],
						bin_op_ptr(std::min< t_float >));
					// hi_bbox_ = max(hi_bbox_, hi_part_bbox)
					std::transform(&hi_bbox_[0], &hi_bbox_[D], &pb->hi_bbox_[0], &hi_bbox_[0],
						bin_op_ptr(std::max< t_float >));
				}
			}
		}

		// find bounds by iterating over all mesh_part cells
		void calc_bbox_raw(vertex_pos& lo_bbox, vertex_pos& hi_bbox) const {
			// init bounds from first cell
			ss(lo).lo(lo_bbox); ss(lo).hi(hi_bbox);

			// walk all other cells
			vertex_pos lo_cell, hi_cell;
			ulong sz = size();
			for(ulong i = 1; i < sz; ++i) {
				const cell_data& cell = ss(ss_id(i));
				cell.lo(lo_cell); cell.hi(hi_cell);
				for(uint i = 0; i < D; ++i) {
					if(lo_cell[i] < lo_bbox[i])
						lo_bbox[i] = lo_cell[i];
					if(hi_cell[i] > hi_bbox[i])
						hi_bbox[i] = hi_cell[i];
				}
			}
		}

		void bounds(vertex_pos& lo_pos, vertex_pos& hi_pos) const {
			vertex_pos lo_lo_pos, lo_hi_pos;
			vertex_pos hi_lo_pos, hi_hi_pos;
			ss(lo).lo (lo_lo_pos);
			ss(lo).hi (lo_hi_pos);

			// last = hi - 1
			vertex_pos_i last;
			ca_assign(last, hi);
			std::transform(&last[0], &last[D], &last[0], std::bind2nd(std::minus< ulong >(), 1));
			ss(last).lo (hi_lo_pos);
			ss(last).hi (hi_hi_pos);

			// ensure that lo[i] < hi[i]
			for(uint i = 0; i < D; ++i) {
				if (lo_lo_pos[i] < hi_lo_pos[i]) {
					// global coordinates increasing along cells - normal order
					lo_pos[i] = lo_lo_pos[i]; // std::min (lo_lo_pos[i], lo_hi_pos[i]);
					hi_pos[i] = hi_hi_pos[i]; // std::max (hi_lo_pos[i], hi_hi_pos[i]);
				}
				else {
					// lo_lo_pos[i] > hi_lo_pos[i] // global coordinates decreasing along cells - reverse order 
					lo_pos[i] = hi_lo_pos[i]; // std::min (hi_lo_pos[i], hi_hi_pos[i]);
					hi_pos[i] = lo_hi_pos[i]; // std::max (lo_lo_pos[i], lo_hi_pos[i]);
				}
			}
		}
	};

	static bool point_inside_bbox(const Bbox& b, const Point& p) {
		for(uint i = 0; i < D; ++i) {
			if(p[i] < b.min(i) || p[i] > b.max(i))
				return false;
		}
		return true;
	}

	static std::vector< ulong > where_is_point(
		trimesh& m,
		std::vector< Point > points)
	{
		// start with full mesh
		// and divide it until we come to only one cell

		// mesh partition stored here
		typedef typename mesh_part::container_t parts_container;
		typedef typename mesh_part::container_t::iterator part_iterator;
		parts_container parts;
		parts.insert(mesh_part(m));
		// precalc mesh size
		const ulong m_size = m.size_flat();

		// found cell_ids stored here
		std::vector< ulong > res(points.size(), -1);
		//ulong cell_id;
		//vertex_pos c_lo, c_hi;
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
				//const Iso_bbox& cur_rect = l->iso_bbox();
				const Bbox& cur_rect = l->bbox();

				catched_points.clear();
				for(ulong i = 0; i < points.size(); ++i) {
					// skip already found points
					if(res[i] < m_size) continue;
					// check that point lies inside this part
					if(point_inside_bbox(cur_rect, points[i])) {
						catched_points.push_back(i);
						// parts with > 1 cell anyway undergo splitting
						if(l->size() > 1)
							break;
					}

					//if(!cur_rect.has_on_unbounded_side(points[i]))
					//	catched_points.push_back(i);
				}

				// if this part don't contain any points - remove it
				// if box contains only 1 cell - test if cell poly contains given points
				if(!catched_points.size())
					leafs.erase(l++);
				else if(l->size() == 1) {
					ulong cell_id = encode_cell_id(l->lo, m.size());
					//Polygon_2 cell_poly = m[cell_id].polygon();
					cell_data cell = m[cell_id];
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

	static ulong where_is_point(trimesh& m, Point point) {
		return where_is_point(m, std::vector< Point >(1, point))[0];
	}

};

}} /* blue_sky::wpi */

#endif /* end of include guard: WPI_ALGO_MESHP_RIC3ZQNS */
