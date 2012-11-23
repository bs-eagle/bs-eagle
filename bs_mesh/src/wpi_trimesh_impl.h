/// @file wpi_trimesh_impl.h
/// @brief Full implementation of wpi::pods::trimesh class
/// @author uentity
/// @version 1.0
/// @date 06.11.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_TRIMESH_IMPL_1I00E960
#define WPI_TRIMESH_IMPL_1I00E960

#include "wpi_algo_pod.h"
#include "wpi_online_tops_traits.h"
// we depend on mesh_grdecl
#include "i_cant_link_2_mesh.h"
#include "rs_smesh_iface.h"
#include "coord_zcorn_tools.h"

namespace blue_sky { namespace wpi {

// implementation of simple approach using continous array of cells tops
template< class strat_t >
struct pods< strat_t >::trimesh::impl {
	typedef typename strat_t::traits_t strat_traits;
	typedef pods< strat_t > pods_t;
	typedef typename pods_t::trimesh trimesh_t;
	typedef typename pods_t::vertex_pos_i vertex_pos_i;
	typedef typename pods_t::cell_vertex_iterator iterator_t;
	typedef typename trimesh_t::value_type value_type;
	//typedef typename trimesh_t::const_reference cref_t;

	enum { D = strat_t::D };

	/*-----------------------------------------------------------------
	 * internals specializations for different traits
	 *----------------------------------------------------------------*/

	// most simple = for raw tops array
	template< class traits_t, class = void >
	struct iimpl {
		typedef typename traits_t::cell_vertex_iterator iterator_t;

		// holder for tops array
		spv_float tops_;

		void init(t_long nx, t_long ny, spv_float coord, spv_float zcorn) {
			using namespace blue_sky;

			// obtain coordinates for all vertices of all cells
			sp_himesh handy = BS_KERNEL.create_object("handy_mesh_iface");
			tops_ = handy->calc_cells_vertices_xyz(nx, ny, coord, zcorn);
		}

		iterator_t begin() const {
			return tops_->begin();
		}
		iterator_t end() const {
			return tops_->end();
		}
		typename v_float::const_iterator data() const {
			return tops_->begin();
		}
	};

	// specialization for online tops traits
	// holds pointer to rs_smesh_iface
	template< class unused >
	struct iimpl< online_tops_traits, unused > {
		typedef typename online_tops_traits::cell_vertex_iterator iterator_t;
		typedef smart_ptr< rs_smesh_iface, true > sp_smesh;

		void init(t_long nx, t_long ny, spv_float coord, spv_float zcorn) {
			using namespace blue_sky;
			// build mesh_grdecl around given mesh
			sp_himesh handy = BS_KERNEL.create_object("handy_mesh_iface");
			mesh_ = handy->make_mesh_grdecl(nx, ny, coord, zcorn);
			assert(mesh_);
		}

		iterator_t begin() const {
			return iterator_t(mesh_);
		}
		iterator_t end() const {
			return iterator_t(mesh_, mesh_->get_n_elements());
		}
		rs_smesh_iface* data() const {
			return mesh_.lock();
		}

		sp_smesh mesh_;
	};

	// specialization for sgrid_ti_traits
	// holds structured grid array of size (nx + 1)*(ny + 1)*(nz + 1)
	template< class traits_t >
	struct iimpl_sgrid {
		typedef typename traits_t::cell_vertex_iterator iterator_t;
		typedef typename iterator_t::strat_ctor_param_t sgrid_handle;

		void init(t_long nx, t_long ny, spv_float coord, spv_float zcorn) {
			using namespace blue_sky::coord_zcorn_tools;
			// extract structured grid from mesh
			sgrid_ = tops2struct_grid(nx, ny, coord, zcorn);
			assert(sgrid_);

			// save handle
			sgrid_handle h = { sgrid_->begin(), ulong(nx), ulong(ny), zcorn->size() >> 3 };
			h_ = h;
		}

		iterator_t begin() const {
			return iterator_t(h_);
		}
		iterator_t end() const {
			return iterator_t(h_, h_.size);
		}
		sgrid_handle data() const {
			return h_;
		}

		spv_float sgrid_;
		sgrid_handle h_;
	};

	template< class unused >
	struct iimpl< sgrid_traits, unused > : public iimpl_sgrid< sgrid_traits > {};

	// bufpool-related strategies can utilize the same code :)
	// can't stop pushing the template limits :)
	// implementation is base on simple ti_traits without _bofpool postfix
	// so pass such strategy as template parameter
	// the purpose of this class is to add cell_buf_storage to original uncached traits
	template< class traits_t >
	struct iimpl_bufpool : public iimpl< typename traits_t::uncached_traits > {
		typedef typename traits_t::cell_vertex_iterator iterator_t;
		typedef typename iterator_t::buf_ctor_param_t cell_buf_storage;
		iimpl< typename traits_t::uncached_traits > base_t;

		using base_t::init;

		iterator_t begin() {
			return iterator_t(std::make_pair(base_t::data(), &store_));
		}
		iterator_t end() {
			return iterator_t(std::make_pair(base_t::data(), &store_),
				base_t::end() - base_t::begin());
		}

		~iimpl_bufpool() {
			// free memory allocated by cell buffers
			typedef typename cell_buf_storage::allocator_type alloc;
			typedef typename alloc::value_type V;
			boost::singleton_pool< boost::fast_pool_allocator_tag, sizeof(V) >::release_memory();
		}

		cell_buf_storage store_;
	};

	// specializations of the above for *_bufpool ti_traits
	template< class unused >
	struct iimpl< online_tops_traits_bufpool, unused >
		: public iimpl_bufpool< online_tops_traits_bufpool >
	{};

	/*-----------------------------------------------------------------
	 * higher-level trimesh::impl definition
	 *----------------------------------------------------------------*/
	// empty ctor
	impl() {}

	void init(t_long nx, t_long ny, spv_float coord, spv_float zcorn) {
		ii_.init(nx, ny, coord, zcorn);
	}

	value_type ss(ulong idx) {
		return value_type(ii_.begin() + 24 * idx);
	}

	iterator_t begin() {
		return ii_.begin();
	}
	iterator_t end() {
		return ii_.end();
	}

	iimpl< strat_traits > ii_;
};

// common implementation of trimesh fuunctions

template< class strat_t >
pods< strat_t >::trimesh::trimesh()
	: pimpl_(new impl)
{
	ca_assign(size_, ulong(0));
}

template< class strat_t >
pods< strat_t >::trimesh::trimesh(t_long nx, t_long ny, spv_float coord, spv_float zcorn)
	: pimpl_(new impl)
{
	this->init(nx, ny, coord, zcorn);
}

template< class strat_t >
void pods< strat_t >::trimesh::init(t_long nx, t_long ny, spv_float coord, spv_float zcorn) {
	pimpl_->init(nx, ny, coord, zcorn);
	// set size
	const ulong full_sz[] = {ulong(nx), ulong(ny), (zcorn->size() / (nx * ny)) >> 3};
	std::copy(full_sz, full_sz + D, size_);
}

template< class strat_t >
typename pods< strat_t >::trimesh::value_type
pods< strat_t >::trimesh::ss_backend(ulong idx) const {
	return pimpl_->ss(idx);
}

template< class strat_t >
typename pods< strat_t >::cell_vertex_iterator
pods< strat_t >::trimesh::begin() const {
	return pimpl_->begin();
}

template< class strat_t >
typename pods< strat_t >::cell_vertex_iterator
pods< strat_t >::trimesh::end() const {
	return pimpl_->end();
}

// instantiate trimesh for known strategies
//template class pods< strategy_3d >::trimesh;
//template class pods< strategy_2d >::trimesh;

}} /* blue_sky::wpi */

#endif /* end of include guard: WPI_TRIMESH_IMPL_1I00E960 */

