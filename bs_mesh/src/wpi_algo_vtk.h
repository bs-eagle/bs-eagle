/// @file wpi_algo_vtk.h
/// @brief VTK-related algorithms
/// @author uentity
/// @version 1.0
/// @date 14.11.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_ALGO_VTK_ZTTIWKX4
#define WPI_ALGO_VTK_ZTTIWKX4

#include "wpi_common.h"
#include "wpi_algo_pod.h"
#include "wpi_trimesh_impl.h"
#include "wpi_algo_meshp.h"

namespace blue_sky { namespace wpi {

template< class strat_t >
struct algo_vtk : helpers< strat_t > {
	typedef typename strat_t::cell_vertex_iterator cell_vertex_iterator;
	typedef typename strat_t::well_traj_iterator   well_traj_iterator;

	// import global consts
	enum { D = strat_t::D };

	// import pods
	typedef pods< strat_t > pods_t;
	typedef typename pods_t::cell_data cell_data;
	typedef typename pods_t::sp_cell_data sp_cell_data;
	typedef typename pods_t::trimesh trimesh;

	// import mesh_part
	typedef mesh_tools< strat_t > mesh_tools_t;
	typedef typename mesh_tools_t::mesh_part mesh_part;

	/*-----------------------------------------------------------------
	 * vtk index storage backend for facets and edges
	 *----------------------------------------------------------------*/

	// helper structure for filtering vertices
	struct vertex_handle {
		// position of vertex coordinates
		cell_vertex_iterator pv_;
		// offset of this handle in vertex storage
		// filled during exporting step
		ulong offs;

		vertex_handle(const cell_vertex_iterator& V, ulong offs_ = -1)
			: pv_(V), offs(offs_)
		{}

		// sorting criteria
		bool operator<(const vertex_handle& rhs) const {
			return std::lexicographical_compare(
				pv_, pv_ + 3,
				rhs.pv_, rhs.pv_ + 3
			);
		}

		const cell_vertex_iterator& ptr() const {
			return pv_;
		}
	};

	typedef std::set< vertex_handle > vertex_storage;
	typedef typename vertex_storage::iterator vertex_iterator;
	typedef typename vertex_storage::const_iterator vertex_citerator;
	typedef std::pair< vertex_iterator, bool > ins_res;

	static void export_vertex_storage(vertex_storage& vs, spv_float points) {
		// export array of unique points and calc their offsets
		points->resize(vs.size() * 3);
		v_float::iterator p_dst = points->begin();
		ulong cnt = 0;
		for(vertex_iterator p_src = vs.begin(), end = vs.end(); p_src != end; ++p_src, ++cnt) {
			p_dst = std::copy(p_src->ptr(), p_src->ptr() + 3, p_dst);
			const_cast< vertex_handle& >(*p_src).offs = cnt;
		}
	}

	// implementation for drawing facets using vtkQuad
	template< int prim_id, class unused = void >
	struct vtk_index_backend {
		enum { type = prim_id };
		enum { n_fv = cell_data::n_facet_vertex };
		typedef typename cell_data::facet_vid_t facet_vid_t;

		typedef bs_array< t_long, vector_traits > bs_lvector;
		//smart_ptr< bs_lvector > idx_;
		smart_ptr< bs_lvector > cell_idx_;
		cell_vertex_iterator tops_;
		vertex_storage vs_;

		// indexes of points just store vertex_storage iterators
		typedef std::list< vertex_iterator > index_storage;
		index_storage idx_;

		// ctor
		vtk_index_backend(const cell_vertex_iterator& Tops)
			: cell_idx_(BS_KERNEL.create_object(bs_lvector::bs_type())),
			tops_(Tops)
		{}

		// try to insert vertex to storage, return index to be pushed to idx_
		ulong push_vertex(ulong vid) {
			ins_res r = vs_.insert(vertex_handle(vid, tops_));
			return r.first->v_id;
		}

		// return how many primitives from given v were inserted into index
		int operator()(const facet_vid_t& v, ulong cell_id, ulong facet_id) {
			for(uint i = 0; i < n_fv; ++i) {
				// save only unique vertices
				idx_.push_back(vs_.insert(vertex_handle(tops_ + v[i]*3)).first);
			}
			// save cell id
			cell_idx_->push_back(cell_id);
			// we're alwais inserting values
			return 1;
		}

		spv_long get(spv_long cell_ids, spv_float points) {
			// export unque points
			export_vertex_storage(vs_, points);

			// export point ids for vtkQuads
			spv_long res = BS_KERNEL.create_object(v_long::bs_type());
			// each facet reqire a prefix denoting number of facet vertices
			res->resize(idx_.size() + cell_idx_->size());
			v_long::iterator p_res = res->begin();
			typename index_storage::const_iterator p_idx = idx_.begin();
			for(ulong i = 0; i < cell_idx_->size(); ++i) {
				// prefix
				*p_res++ = n_fv;
				for(uint j = 0; j < n_fv; ++j, ++p_res, ++p_idx)
					*p_res = (*p_idx)->offs;
			}

			// cell_ids go inplace
			cell_ids->init_inplace(cell_idx_);
			return res;
		}
	};

	// specialization for edges
	// implementation for drawing facets using vtkQuad
	template< class unused >
	struct vtk_index_backend< 1, unused > {
		enum { type = 1 };
		enum { n_fv = cell_data::n_facet_vertex };
		typedef typename cell_data::facet_vid_t facet_vid_t;

		//typedef bs_array< t_long, vector_traits > bs_lvector;
		//smart_ptr< bs_lvector > idx_;

		cell_vertex_iterator tops_;
		vertex_storage vs_;

		// helper structure for filtering edges
		struct edge_handle {
			ulong cell_id, facet_id;
			// position to the beginning and end of edge
			vertex_iterator beg_, end_;

			edge_handle(ulong cell_id_, ulong facet_id_,
				const vertex_iterator& start, const vertex_iterator& stop)
				: cell_id(cell_id_), facet_id(facet_id_), beg_(start), end_(stop)
			{
				if(!(*beg_ < *end_))
					// v1_ >= v_2
					std::swap(beg_, end_);
			}

			// sorting criteria
			bool operator<(const edge_handle& rhs) const {
				const int r = lexicographical_compare_3way(
					beg_->ptr(), beg_->ptr() + 3, rhs.beg_->ptr(), rhs.beg_->ptr() + 3
				);
				if(r == 0)
					// v1 = rhs.v1, so check v2
					return *end_ < *rhs.end_;
				else
					return (r < 0);
			}
		};

		typedef std::set< edge_handle > edge_storage;
		typedef typename edge_storage::iterator edge_iterator;
		typedef typename edge_storage::const_iterator edge_citerator;
		typedef std::pair< edge_iterator, bool > ins_res;
		edge_storage es_;

		// ctor
		vtk_index_backend(const cell_vertex_iterator& Tops) : tops_(Tops) {}

		inline int push_edge(const edge_handle& e) {
			ins_res r = es_.insert(e);
			// second insertion from DIFFERENT cell kills edge IF (!) we insert the same facet
			if(!r.second && e.cell_id != r.first->cell_id && e.facet_id == r.first->facet_id) {
				es_.erase(r.first);
				return -1;
			}
			return int(r.second);
		}

		// return how many primitives from given v were inserted into index
		int operator()(const facet_vid_t& v, ulong cell_id, ulong facet_id) {
			// try to insert each edge of given facet
			int cnt = 0;
			for(uint i = 0; i < n_fv - 1; ++i) {
				cnt += push_edge(edge_handle(
					cell_id, facet_id,
					vs_.insert(tops_ + v[i]*3).first, vs_.insert(tops_ + v[i + 1]*3).first
				));
			}
			if(n_fv > 2)
				cnt += push_edge(edge_handle(
					cell_id, facet_id,
					vs_.insert(tops_ + v[0]*3).first, vs_.insert(tops_ + v[n_fv - 1]*3).first
				));

			return cnt;
		}

		spv_long get(spv_long cell_ids, spv_float points) {
			// export unque points
			export_vertex_storage(vs_, points);

			spv_long res = BS_KERNEL.create_object(v_long::bs_type());
			// copy collected unique edges as lines to resulting index
			res->resize(es_.size() * 3);
			cell_ids->resize(es_.size());
			v_long::iterator pres = res->begin();
			v_long::iterator pcell = cell_ids->begin();
			for(edge_citerator pe = es_.begin(), end = es_.end(); pe != end; ++pe, ++pcell) {
				*pres++ = 2;
				*pres++ = pe->beg_->offs; *pres++ = pe->end_->offs;
				*pcell = pe->cell_id;
			}
			return res;
		}
	};

	/*-----------------------------------------------------------------
	 * Enumerate border cell facets for drawing mesh in VTK
	 *----------------------------------------------------------------*/
	template< int prim_id >
	static spv_long enum_border_vtk(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_int mask,
		spv_long cell_idx, spv_float points, Loki::Int2Type< prim_id > prim)
	{
		// 1) build trimesh from given tops
		trimesh M(nx, ny, coord, zcorn);
		const ulong n_cells = M.size_flat();
		//vertex_pos_i mesh_size = {ulong(nx), ulong(ny), tops->size() / (nx * ny * 24)};
		//const ulong n_cells = tops->size() / 24;
		//M.resize(n_cells);
		//v_float::iterator pv = tops->begin();
		//for(ulong i = 0; i < n_cells; ++i) {
		//	M[i] = cell_data(&*pv);
		//	pv += 3*8;
		//}

		// make mesh_part containing full mesh
		mesh_part MP(M);

		// 2) loop over all cells
		typedef typename mesh_part::cell_neighb_enum cell_nb_enum;
		typedef typename cell_data::facet_vid_t facet_vid_t;
		enum { n_facets = cell_data::n_facets };
		enum { n_fv = cell_data::n_facet_vertex };
		const ulong mask_sz = mask->size();

		// backend for storing index
		vtk_index_backend< prim_id > vib(M.begin());

		cell_nb_enum cell_nb;
		facet_vid_t cell_fvid;

		// loop over all cells
		for(ulong i = 0; i < n_cells; ++i) {
			// skip masked cells
			if(i < mask_sz && mask->ss(i) == 0)
				continue;

			// 2.1) if some facet has no neighbors - include it in results
			MP.cell_neighbours(i, cell_nb);
			for(ulong j = 0; j < n_facets; ++j) {
				// skip facet if it has non-masked neighbour
				if(cell_nb[j] < n_cells && (cell_nb[j] >= mask_sz || mask->ss(cell_nb[j]) != 0))
					continue;

				cell_data::facet_vid(j, cell_fvid);
				// vertex_id[i] = vertex_id[i] + i*8
				std::transform(
					&cell_fvid[0], &cell_fvid[n_fv], &cell_fvid[0],
					std::bind2nd(std::plus< ulong >(), i * 8)
				);
				vib(cell_fvid, i, j);
			}
		}

		//cell_idx->clear();
		return vib.get(cell_idx, points);
	}

	// specialization for facets
	static spv_long enum_border_facets_vtk(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_int mask,
		spv_long cell_idx, spv_float points) {
		return enum_border_vtk(nx, ny, coord, zcorn, mask, cell_idx, points, Loki::Int2Type< 0 >());
	}

	// specialization for edges
	static spv_long enum_border_edges_vtk(t_long nx, t_long ny, spv_float coord, spv_float zcorn, spv_int mask,
		spv_long cell_idx, spv_float points) {
		return enum_border_vtk(nx, ny, coord, zcorn, mask, cell_idx, points, Loki::Int2Type< 1 >());
	}

}; // algo_vtk

}} /* blue_sky::wpi */

#endif /* end of include guard: WPI_ALGO_VTK_ZTTIWKX4 */

