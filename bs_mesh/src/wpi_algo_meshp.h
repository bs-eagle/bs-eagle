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

namespace blue_sky { namespace wpi {

template< class strat_t >
struct wpi_algo_meshp : public wpi_algo_helpers< strat_t > {
	// basic types
	typedef typename strat_t::vertex_pos   vertex_pos;
	typedef typename strat_t::vertex_pos_i vertex_pos_i;

	typedef typename strat_t::Point    Point;
	typedef typename strat_t::Iso_bbox Iso_bbox;
	typedef typename strat_t::Bbox     Bbox;

	// import global consts
	enum { D = strat_t::D };

	// import pods
	typedef wpi_algo_pod< strat_t > wpi_pod;
	typedef typename wpi_pod::cell_data cell_data;
	typedef typename wpi_pod::well_data well_data;
	typedef typename wpi_pod::trimesh trimesh;
	typedef typename wpi_pod::trim_iterator trim_iterator;

	/*-----------------------------------------------------------------
	* represent rectangular part of mesh with splitting support
	*----------------------------------------------------------------*/
	// x_last = last_element + 1 = x_size
	// y_last = last_element + 1 = y_size
	struct mesh_part : public wpi_algo_helpers< strat_t > {
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

		void init(ulong lower, ulong upper) {
			std::vector< ulong > idx(2);
			idx[0] = lower; idx[1] = upper;
			init(idx);
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

		ulong ss_id(const vertex_pos_i& offset) const {
			vertex_pos_i cell;
			ca_assign(cell, lo);
			for(uint i = 0; i < D; ++i)
				cell[i] += offset[i];
			return encode_cell_id(cell, m_size_);
		}

		ulong ss_id(ulong offset) const {
			// size of this part
			vertex_pos_i part_size;
			for(uint i = 0; i < D; ++i)
				part_size[i] = side_len(i);

			// plain id -> vertex_pos_i
			vertex_pos_i part_pos;
			decode_cell_id(offset, part_pos, part_size);
			// part_pos -> cell
			return ss_id(part_pos);
		}

		trim_iterator ss_iter(const vertex_pos_i& offset) {
			return m_.begin() + ss_id(offset);
		}

		trim_iterator ss_iter(ulong offset) {
			return m_.begin() + ss_id(offset);
		}

		Iso_bbox iso_bbox() const {
			vertex_pos lo_pos, hi_pos;
			bounds(lo_pos, hi_pos);
			return vertex_pos2rect(lo_pos, hi_pos);
		}

		Bbox bbox() const {
			vertex_pos lo_pos, hi_pos;
			bounds(lo_pos, hi_pos);
			return vertex_pos2bbox(lo_pos, hi_pos);
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

		const cell_data& ss(ulong idx) const {
			// idx SHOULD BE IN MESH!
			return m_[idx];
		}

		const cell_data& ss(const vertex_pos_i& idx) const {
			return m_[encode_cell_id(idx, m_size_)];
		}

		mesh_part(trimesh& m, const vertex_pos_i& mesh_size,
				const vertex_pos_i& first_,
				const vertex_pos_i& last_)
			: m_(m)
		{
			ca_assign(lo, first_);
			ca_assign(hi, last_);
			ca_assign(m_size_, mesh_size);
		}

		void bounds(vertex_pos& lo_pos, vertex_pos& hi_pos) const {
			ss(lo).lo(lo_pos);
			// last = hi - 1
			vertex_pos_i last;
			ca_assign(last, hi);
			std::transform(&last[0], &last[D], &last[0], bind2nd(std::minus< ulong >(), 1));
			ss(last).hi(hi_pos);
		}
	};

	static std::vector< ulong > where_is_point(
		trimesh& m, const vertex_pos_i& m_size,
		std::vector< Point > points)
	{
		// start with full mesh
		// and divide it until we come to only one cell

		// mesh partition stored here
		typedef typename mesh_part::container_t parts_container;
		typedef typename mesh_part::container_t::iterator part_iterator;
		parts_container parts;
		parts.insert(mesh_part(m, m_size));

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
				const Iso_bbox& cur_rect = l->iso_bbox();
				catched_points.clear();
				for(ulong i = 0; i < points.size(); ++i) {
					// skip already found points
					if(res[i] < m.size()) continue;
					// check that point lies inside this part
					if(!cur_rect.has_on_unbounded_side(points[i]))
						catched_points.push_back(i);
				}

				// if this part don't contain any points - remove it
				// if box contains only 1 cell - test if cell poly contains given points
				if(!catched_points.size())
					leafs.erase(l++);
				else if(l->size() == 1) {
					ulong cell_id = encode_cell_id(l->lo, m_size);
					//Polygon_2 cell_poly = m[cell_id].polygon();
					cell_data& cell = m[cell_id];
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

	static ulong where_is_point(trimesh& m, const vertex_pos_i& m_size, Point point) {
		return where_is_point(m, m_size, std::vector< Point >(1, point))[0];
	}

};

}} /* blue_sky::wpi */

#endif /* end of include guard: WPI_ALGO_MESHP_RIC3ZQNS */

