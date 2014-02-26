/// @file wpi_algo_trimesh.h
/// @brief Mesh abstraction layer for WPI algorithms
/// @author uentity
/// @version 1.0
/// @date 04.12.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_ALGO_TRIMESH_TD8OAAW1
#define WPI_ALGO_TRIMESH_TD8OAAW1

#include "wpi_common.h"
#include "wpi_strategy_traits.h"
#include "i_cant_link_2_mesh.h"
#include "rs_smesh_iface.h"
#include "coord_zcorn_tools.h"

namespace blue_sky { namespace wpi {
// implementation details
namespace detail {
/*-----------------------------------------------------------------
 * internals specializations for different traits
 *----------------------------------------------------------------*/
// most simple = for raw tops array
template< class traits_t >
struct trimpl {
	typedef typename traits_t::cell_vertex_iterator iterator_t;

	// holder for tops array
	spv_float tops_;

	void init(t_long nx, t_long ny, spv_float coord, spv_float zcorn) {
		using namespace blue_sky;

		// obtain coordinates for all vertices of all cells
		sp_himesh handy = BS_KERNEL.create_object("handy_mesh_iface");
		tops_ = handy->calc_cells_vertices_xyz(nx, ny, coord, zcorn);
	}

	void init(t_long, t_long, const sp_obj& backend) {
		tops_ = backend;
		if(!tops_)
			throw bs_exception("trimesh::impl", "Incompatible backend passed");
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
	sp_obj backend() const {
		return tops_;
	}
};

// specialization for online tops traits
// holds pointer to rs_smesh_iface
template<  >
struct trimpl< online_tops_traits > {
	typedef typename online_tops_traits::cell_vertex_iterator iterator_t;
	typedef smart_ptr< rs_smesh_iface, true > sp_smesh;

	void init(t_long nx, t_long ny, spv_float coord, spv_float zcorn) {
		using namespace blue_sky;
		// build mesh_grdecl around given mesh
		sp_himesh handy = BS_KERNEL.create_object("handy_mesh_iface");
		mesh_ = handy->make_mesh_grdecl(nx, ny, coord, zcorn);
		assert(mesh_);
	}

	void init(t_long, t_long, const sp_obj& backend) {
		mesh_ = backend;
		if(!mesh_)
			throw bs_exception("trimesh::impl", "Incompatible backend passed");
	}

	iterator_t begin() const {
		return iterator_t(mesh_);
	}
	iterator_t end() const {
		return iterator_t(mesh_, mesh_->get_n_elements());
	}
	rs_smesh_iface* data() const {
#ifdef BS_DISABLE_MT_LOCKS
		return mesh_.get();
#else
		return mesh_.lock();
#endif
	}
	sp_obj backend() const {
		return mesh_;
	}

	sp_smesh mesh_;
};

// specialization for sgrid_ti_traits
// holds structured grid array of size (nx + 1)*(ny + 1)*(nz + 1)
template< class traits_t >
struct trimpl_sgrid {
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

	void init(t_long nx, t_long ny, const sp_obj& backend) {
		sgrid_ = backend;
		if(!sgrid_)
			throw bs_exception("trimesh::impl", "Incompatible backend passed");

		// save handle
		sgrid_handle h = {
			sgrid_->begin(), ulong(nx), ulong(ny), sgrid_->size() / ((nx + 1)*(ny + 1)) - 1
		};
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
	sp_obj backend() const {
		return sgrid_;
	}

	spv_float sgrid_;
	sgrid_handle h_;
};

template<  >
struct trimpl< sgrid_traits > : public trimpl_sgrid< sgrid_traits > {};

// bufpool-related strategies can utilize the same code :)
// can't stop pushing the template limits :)
// implementation is base on simple ti_traits without _bofpool postfix
// so pass such strategy as template parameter
// the purpose of this class is to add cell_buf_storage to original uncached traits
template< class traits_t >
struct trimpl_bufpool : public trimpl< typename traits_t::uncached_traits > {
	typedef typename traits_t::cell_vertex_iterator iterator_t;
	typedef typename iterator_t::buf_ctor_param_t cell_buf_storage;
	typedef trimpl< typename traits_t::uncached_traits > base_t;

	using base_t::init;
	using base_t::data;

	iterator_t begin() {
		return iterator_t(std::make_pair(data(), &store_));
	}
	iterator_t end() {
		return iterator_t(std::make_pair(data(), &store_),
			base_t::end() - base_t::begin());
	}

	~trimpl_bufpool() {
		// free memory allocated by cell buffers
		typedef typename cell_buf_storage::allocator_type alloc;
		typedef typename alloc::value_type V;
		boost::singleton_pool< boost::fast_pool_allocator_tag, sizeof(V) >::release_memory();
	}

	cell_buf_storage store_;
};

// specializations of the above for *_bufpool ti_traits
template<  >
struct trimpl< online_tops_traits_bufpool >
	: public trimpl_bufpool< online_tops_traits_bufpool > {};

} // eof detail namespace

/*-----------------------------------------------------------------
 * Main trimesh class declaration
 *----------------------------------------------------------------*/
template< class strat_traits >
class Trimesh : public detail::trimpl< strat_traits > {
public:
	typedef strat_traits traits_t;
	typedef detail::trimpl< traits_t > base_t;
	typedef ulong vertex_pos_i[3];
	typedef typename traits_t::cell_vertex_iterator cell_vertex_iterator;

	enum { n_cell_coords = sizeof(cell_pos) / sizeof(t_float) };

	using base_t::init;
	using base_t::begin;

	// empty ctor
	Trimesh() {
		ca_assign(size_, ulong(0));
	}

	// ctor from given COORD & ZCORN
	Trimesh(t_long nx, t_long ny, spv_float coord, spv_float zcorn) {
		init(nx, ny, coord, zcorn);
	}

	const vertex_pos_i& size() const {
		return size_;
	}

	ulong size_flat() const {
		ulong sz = 1;
		for(uint i = 0; i < 3; ++i)
			sz *= size_[i];
		return sz;
	}

	// subscripting - returns a new _copy_ of cell every time!
	template< class cell_data_t >
	cell_data_t ss(ulong idx) const {
		return cell_data_t(begin() + n_cell_coords * idx);
	}

	// the same in operator[] form
	template< class cell_data_t >
	cell_data_t operator[](ulong idx) const {
		return ss(idx);
	}

	static sp_obj create_backend(t_long nx, t_long ny, spv_float coord, spv_float zcorn) {
		Trimesh M(nx, ny, coord, zcorn);
		return M.backend();
	}

private:
	//struct impl;
	//st_smart_ptr< impl > pimpl_;

	vertex_pos_i size_;
};

/*-----------------------------------------------------------------
 * Definition of trimesh functions
 *----------------------------------------------------------------*/
//template< class strat_traits >
//trimesh< strat_traits >::trimesh()
//	: pimpl_(new impl)
//{
//	ca_assign(size_, ulong(0));
//}
//
//template< class strat_traits >
//trimesh< strat_traits >::trimesh(t_long nx, t_long ny, spv_float coord, spv_float zcorn)
//	: pimpl_(new impl)
//{
//	this->init(nx, ny, coord, zcorn);
//}
//
//template< class strat_traits >
//void trimesh< strat_traits >::init(t_long nx, t_long ny, spv_float coord, spv_float zcorn) {
//	pimpl_->init(nx, ny, coord, zcorn);
//	// set size
//	const ulong full_sz[] = {ulong(nx), ulong(ny), (zcorn->size() / (nx * ny)) >> 3};
//	std::copy(full_sz, full_sz + 3, size_);
//}
//
//template< class strat_traits >
//void trimesh< strat_traits >::init(t_long nx, t_long ny, sp_obj backend) {
//	return pimpl_->init(nx, ny, backend);
//}
//
//template< class strat_traits >
//typename trimesh< strat_traits >::cell_vertex_iterator
//trimesh< strat_traits >::begin() const {
//	return pimpl_->begin();
//}
//
//template< class strat_traits >
//typename trimesh< strat_traits >::cell_vertex_iterator
//trimesh< strat_traits >::end() const {
//	return pimpl_->end();
//}
//
//template< class strat_traits >
//sp_obj trimesh< strat_traits >::backend() const {
//	return pimpl_->backend();
//}

}}  /* blue_sky::wpi */

#endif /* end of include guard: WPI_ALGO_TRIMESH_TD8OAAW1 */

