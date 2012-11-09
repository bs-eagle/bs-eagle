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
//#include "wpi_strategy_3d.h"
//#include "wpi_strategy_2d.h"
// we depend on mesh_grdecl
#include "i_cant_link_2_mesh.h"

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

	// internals specializations for different traits
	template< class traits_t, class = void >
	struct iimpl {
		// holder for tops array
		spv_float tops_;

		void init(t_long nx, t_long ny, spv_float coord, spv_float zcorn, vertex_pos_i& mesh_size) {
			using namespace blue_sky;

			// obtain coordinates for all vertices of all cells
			sp_himesh handy = BS_KERNEL.create_object("handy_mesh_iface");
			tops_ = handy->calc_cells_vertices_xyz(nx, ny, coord, zcorn);
			// set size
			const ulong full_sz[] = {ulong(nx), ulong(ny), (zcorn->size() / (nx * ny)) >> 3};
			std::copy(full_sz, full_sz + D, mesh_size);
		}

		iterator_t begin() const {
			return tops_->begin();
		}
		iterator_t end() const {
			return tops_->end();
		}
	};

	// empty ctor
	impl() {}

	void init(t_long nx, t_long ny, spv_float coord, spv_float zcorn, vertex_pos_i& mesh_size) {
		ii_.init(nx, ny, coord, zcorn, mesh_size);
	}

	value_type ss(ulong idx) const {
		return value_type(ii_.begin() + 24 * idx);
	}

	iterator_t begin() const {
		return ii_.begin();
	}
	iterator_t end() const {
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
	pimpl_->init(nx, ny, coord, zcorn, size_);
}

template< class strat_t >
void pods< strat_t >::trimesh::init(t_long nx, t_long ny, spv_float coord, spv_float zcorn) {
	pimpl_->init(nx, ny, coord, zcorn, size_);
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

