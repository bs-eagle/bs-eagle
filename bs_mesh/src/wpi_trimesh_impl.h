/// @file wpi_trimesh_impl.h
/// @brief Full implementation of wpi::pods::trimesh class
/// @author uentity
/// @version 1.0
/// @date 06.11.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_TRIMESH_IMPL_1I00E960
#define WPI_TRIMESH_IMPL_1I00E960

#include "wpi_common.h"
#include "wpi_strategy_traits.h"
#include "i_cant_link_2_mesh.h"
#include "rs_smesh_iface.h"
#include "coord_zcorn_tools.h"

namespace blue_sky { namespace detail { namespace wpi {
using namespace blue_sky::wpi;
using blue_sky::wpi::ulong;
/*-----------------------------------------------------------------
 * internals trimesh implementation for different traits
 *----------------------------------------------------------------*/
// most simple = for raw tops array
template< class traits_t >
struct trimpl {
	typedef typename traits_t::cell_vertex_iterator iterator_t;

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

protected:
	// holder for tops array
	spv_float tops_;

	void init(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn) {
		using namespace blue_sky;

		// obtain coordinates for all vertices of all cells
		sp_himesh handy = BS_KERNEL.create_object("handy_mesh_iface");
		tops_ = handy->calc_cells_vertices_xyz(nx, ny, coord, zcorn);
		assert(tops_);
	}

	// init returns nz
	ulong init(t_ulong nx, t_ulong ny, const sp_obj& backend) {
		// backend can come from Python and be array with different traits
		// so try to steal buffer instead of direct assignment
		//sgrid_ = backend;
		smart_ptr< bs_arrbase< t_float > > cont(backend, bs_dynamic_cast());
		if(!cont)
			throw bs_exception("trimesh::impl", "Incompatible backend passed");
		tops_ = BS_KERNEL.create_object(v_float::bs_type());
		tops_->init_inplace(cont);

		// for 2D strategy we'll return 0 here, that is OK
		// for 3D returned value is valid nz
		return tops_->size() / ulong(nx * ny * 24);
	}
};

// specialization for online tops traits
// holds pointer to rs_smesh_iface
// NOTE: this implementation most probably is INCOMPATIBLE with 2D strategies!
// because we're using mesh_grdecl inside, which works only for 3D
template< uint D >
struct trimpl< online_tops_traits< D > > {
	typedef typename online_tops_traits< D >::cell_vertex_iterator iterator_t;
	typedef smart_ptr< rs_smesh_iface, true > sp_smesh;

	iterator_t begin() const {
		return iterator_t(mesh_);
	}
	iterator_t end() const {
		return iterator_t(mesh_, mesh_->get_n_elements());
	}
	rs_smesh_iface* data() const {
		return mesh_.lock();
	}
	sp_obj backend() const {
		return mesh_;
	}

protected:
	sp_smesh mesh_;

	void init(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn) {
		using namespace blue_sky;
		// build mesh_grdecl around given mesh
		sp_himesh handy = BS_KERNEL.create_object("handy_mesh_iface");
		mesh_ = handy->make_mesh_grdecl(nx, ny, coord, zcorn);
		assert(mesh_);
	}

	ulong init(t_ulong nx, t_ulong ny, const sp_obj& backend) {
		mesh_ = backend;
		if(!mesh_)
			throw bs_exception("trimesh::impl", "Incompatible backend passed");
		return ulong(mesh_->get_n_elements() / (nx * ny));
	}
};

// specialization for sgrid_ti_traits
// holds structured grid array of size (nx + 1)*(ny + 1)*(nz + 1)
template< class traits_t >
struct trimpl_sgrid {
	typedef typename traits_t::cell_vertex_iterator iterator_t;
	typedef typename iterator_t::strat_t::sgrid_handle sgrid_handle;

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

protected:
	spv_float sgrid_;
	sgrid_handle h_;

	enum { D = traits_t::D };

	void init(t_long nx, t_long ny, spv_float coord, spv_float zcorn) {
		using namespace blue_sky::coord_zcorn_tools;
		// extract structured grid from mesh
		sgrid_ = tops2struct_grid(nx, ny, coord, zcorn);
		assert(sgrid_);

		// save handle
		h_.init(sgrid_->begin(), ulong(nx), ulong(ny), zcorn->size() >> 3);
	}

	ulong init(t_long nx, t_long ny, const sp_obj& backend) {
		// backend can come from Python and be array with different traits
		// so try to steal buffer instead of direct assignment
		//sgrid_ = backend;
		smart_ptr< bs_arrbase< t_float > > cont(backend, bs_dynamic_cast());
		if(!cont)
			throw bs_exception("trimesh::impl", "Incompatible backend passed");
		sgrid_ = BS_KERNEL.create_object(v_float::bs_type());
		sgrid_->init_inplace(cont);

		// sanity
		const ulong nz = sgrid_->size() / ((nx + 1)*(ny + 1)*3) - 1;
		if(nz >= sgrid_->size())
			throw bs_exception("trimesh::impl", "Sgrid backend with wrong dimensions passed");

		// save handle
		if(D == 2) {
			h_.init(sgrid_->begin(), ulong(nx), ulong(ny), ulong(nx * ny));
		}
		else {
			h_.init(sgrid_->begin(), ulong(nx), ulong(ny), ulong(nx * ny * nz));
		}
		return nz;
	}
};

template< uint D >
struct trimpl< sgrid_traits< D > > : public trimpl_sgrid< sgrid_traits< D > > {};

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

	iterator_t begin() const {
		return iterator_t(std::make_pair(data(), &store_));
	}
	iterator_t end() const {
		return iterator_t(std::make_pair(data(), &store_),
			base_t::end() - base_t::begin());
	}

	~trimpl_bufpool() {
		// free memory allocated by cell buffers
		typedef typename cell_buf_storage::allocator_type alloc;
		typedef typename alloc::value_type V;
		boost::singleton_pool< boost::fast_pool_allocator_tag, sizeof(V) >::release_memory();
	}

	mutable cell_buf_storage store_;
};

// specializations of the above for *_bufpool ti_traits
template< uint D >
struct trimpl< online_tops_traits_bufpool< D > >
	: public trimpl_bufpool< online_tops_traits_bufpool< D > > {};

}}} /* blue_sky::wpi::detail */

#endif /* end of include guard: WPI_TRIMESH_IMPL_1I00E960 */

