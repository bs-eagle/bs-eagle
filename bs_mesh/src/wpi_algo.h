/// @file wpi_algo.h
/// @brief Algorithms implementing well path identification using given strategy
/// @author uentity
/// @version 
/// @date 15.09.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_ALGO_I21Y0RBS
#define WPI_ALGO_I21Y0RBS

#include <iterator>
#include <cmath>
#include <algorithm>

#include "wpi_common.h"
#include "wpi_algo_pod.h"
#include "wpi_trimesh_impl.h"
#include "wpi_algo_meshp.h"
#include "wpi_algo_xaction.h"
#include "wpi_algo_xaction_build.h"
#include "wpi_algo_xaction_build2.h"
#include "wpi_algo_xaction_build3.h"

#include "conf.h"
#include "i_cant_link_2_mesh.h"

// DEBUG
//#include <iostream>

namespace blue_sky { namespace wpi {

/*-----------------------------------------------------------------
 * implement well path identification algos depending on strategy
 *----------------------------------------------------------------*/
template< class strat_t >
struct algo : public helpers< strat_t > {
	// import strategy typedefs
	typedef typename strat_t::Point    Point;
	typedef typename strat_t::Segment  Segment;
	typedef typename strat_t::Bbox     Bbox;
	typedef typename strat_t::Iso_bbox Iso_bbox;

	typedef typename strat_t::vertex_pos   vertex_pos;
	typedef typename strat_t::vertex_pos_i vertex_pos_i;

	typedef typename strat_t::cell_vertex_iterator cell_vertex_iterator;
	typedef typename strat_t::well_traj_iterator   well_traj_iterator;

	// import global consts
	enum { D = strat_t::D };

	// import pods
	typedef pods< strat_t > pods_t;
	typedef typename pods_t::cell_data cell_data;
	typedef typename pods_t::sp_cell_data sp_cell_data;
	typedef typename pods_t::trimesh trimesh;
	//typedef typename pods_t::trim_iterator trim_iterator;
	//typedef typename pods_t::ctrim_iterator ctrim_iterator;

	typedef typename pods_t::well_data well_data;
	typedef typename pods_t::well_path well_path;
	typedef typename pods_t::wp_iterator wp_iterator;
	typedef typename pods_t::cwp_iterator cwp_iterator;

	typedef typename pods_t::well_hit_cell well_hit_cell;
	typedef typename pods_t::intersect_path intersect_path;

	// import mesh_part
	typedef mesh_tools< strat_t > mesh_tools_t;
	typedef typename mesh_tools_t::mesh_part mesh_part;

	// import intersect_action
	typedef intersect_base< strat_t > xbase;
	typedef typename xbase::hit_idx_t hit_idx_t;
	typedef intersect_builder2< strat_t > xbuilder;

	// helper to create initial cell_data for each cell
	//static spv_float coord_zcorn2trimesh(t_long nx, t_long ny, spv_float coord, spv_float zcorn,
	//		trimesh& res, vertex_pos_i& mesh_size, bool free_cz_mem = false)
	//{
	//	//typedef smart_ptr< bs_mesh_grdecl, true > sp_grd_mesh;
	//	// build mesh_grdecl around given mesh
	//	//sp_grd_mesh grd_src = BS_KERNEL.create_object(bs_mesh_grdecl::bs_type());
	//	//grd_src->init_props(nx, ny, coord, zcorn);

	//	// init mesh size
	//	const ulong full_sz[] = {ulong(nx), ulong(ny), (zcorn->size() / (nx * ny)) >> 3};
	//	const ulong n_cells = ulong(full_sz[0] * full_sz[1] * full_sz[2]);
	//	std::copy(full_sz, full_sz + D, mesh_size);

	//	// obtain coordinates for all vertices of all cells
	//	sp_himesh handy = BS_KERNEL.create_object("handy_mesh_iface");
	//	spv_float tops = handy->calc_cells_vertices_xyz(nx, ny, coord, zcorn);
	//	// clear COORD & ZCORN arrays
	//	if(free_cz_mem) {
	//		spv_float t = BS_KERNEL.create_object(v_float::bs_type());
	//		t->swap(*coord);
	//		t = BS_KERNEL.create_object(v_float::bs_type());
	//		t->swap(*zcorn);
	//	}

	//	// fill trimesh with triangles corresponding to each cell
	//	res.resize(n_cells);
	//	v_float::iterator pv = tops->begin();
	//	for(ulong i = 0; i < n_cells; ++i) {
	//		// DEBUG
	//		//if(i < 100) {
	//		//	std::cout << std::fixed << std::setprecision(2);
	//		//	for(uint j = 0; j < 24; ++j)
	//		//		std::cout << *(pv + j) << ' ';
	//		//	std::cout << std::endl;
	//		//}
	//		res[i] = cell_data(&*pv);
	//		pv += 3*8;
	//	}

	//	return tops;
	//}

	static ulong fill_well_path(spv_float well_info, well_path& W) {
		ulong well_node_num = well_info->size() >> 2;
		if(well_node_num < 2) return 0;

		// storage
		W.resize(well_node_num - 1);

		// walk along well
		v_float::iterator pw = well_info->begin();
		W[0] = well_data(pw);
		//well_data wd;
		for(ulong i = 1; i < well_node_num - 1; ++i) {
			pw += 4;
			//if(i)
			//	wd = well_data(pw, &W[i - 1]);
			//else
			//	wd = well_data(pw);
	
			// insert well segment
			W[i] = well_data(pw, &W[i - 1]);
			//pw += 4;
		}

		return well_node_num;
	}

	/*-----------------------------------------------------------------
	* implementation of main routine
	*----------------------------------------------------------------*/
	template< bool pythonish, class = void >
	struct wpi_return {
		typedef spv_float type;

		static type make(xbase& A) {
			return A.export_1d();
		}
	};

	template< class unused >
	struct wpi_return< false, unused > {
		typedef std::vector< well_hit_cell > type;

		static type make(xbase& A) {
			type res(A.path().size());
			ulong i = 0;
			for(typename intersect_path::const_iterator px = A.path().begin(), end = A.path().end(); px != end; ++px)
				res[i++] = *px;

			return res;
		}
	};

	template< bool pythonish >
	static typename wpi_return< pythonish >::type well_path_ident_d(
		t_long nx, t_long ny, spv_float coord, spv_float zcorn,
		spv_float well_info, bool include_well_nodes)
	{
		typedef typename wpi_return< pythonish >::type ret_t;

		// 1) calculate mesh nodes coordinates and build initial trimesh
		trimesh M(nx, ny, coord, zcorn);
		//vertex_pos_i mesh_size;
		//spv_float tops = coord_zcorn2trimesh(nx, ny, coord, zcorn, M, mesh_size);
		// DEBUG
		//std::cout << "trimesh built" << std::endl;

		// 2) create well path description
		well_path W;
		if(!fill_well_path(well_info, W)) return ret_t();
		// DEBUG
		//std::cout << "well_path created" << std::endl;

		// 3) construct main object
		xbuilder A(M, W);
		// DEBUG
		//std::cout << "hit_idx found" << std::endl;

		// 4) narrow search space via branch & bound algo
		hit_idx_t& hit_idx = A.build();
		// DEBUG
		//std::cout << "build() done" << std::endl;

		// 5) remove duplicates in X,Y,Z directions
		//A.remove_dups2();
		// DEBUG
		//std::cout << "remove_dups2 done" << std::endl;

		// 6) finalize intersection
		if(include_well_nodes)
			A.append_wp_nodes(hit_idx);
		// DEBUG
		//std::cout << "well nodes inserted" << std::endl;

		return wpi_return< pythonish >::make(A);
	}

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
};

}} /* blue_sky::wpi */

#endif /* end of include guard: WPI_ALGO_I21Y0RBS */

