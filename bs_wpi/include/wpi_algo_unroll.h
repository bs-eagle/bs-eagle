/// @file wpi_algo_unroll.h
/// @brief Unroll well trajectory into 2D grid
/// @author uentity
/// @version 1.0
/// @date 18.04.2013
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_ALGO_UNROLL_YAUH9FX2
#define WPI_ALGO_UNROLL_YAUH9FX2

#include "wpi_common.h"
#include "wpi_algo_pod.h"
#include "wpi_trimesh_impl.h"
#include "wpi_algo_meshp.h"

namespace blue_sky { namespace wpi {

template< class strat_t >
struct algo_unroll : helpers< strat_t > {
	typedef typename strat_t::cell_vertex_iterator cell_vertex_iterator;
	typedef typename strat_t::well_traj_iterator   well_traj_iterator;

	// import global consts
	enum { D = strat_t::D };

	// import pods
	typedef pods< strat_t > pods_t;
	typedef typename pods_t::cell_data cell_data;
	typedef typename pods_t::sp_cell_data sp_cell_data;
	typedef typename pods_t::trimesh trimesh;
	typedef typename pods_t::vertex_pos_i vertex_pos_i;

	// import mesh_part
	typedef mesh_tools< strat_t > mesh_tools_t;
	typedef typename mesh_tools_t::mesh_part mesh_part;

	typedef helpers< strat_t > base_t;
	using base_t::decode_cell_id;

	/*-----------------------------------------------------------------
	 * main well path unrolling function
	 *----------------------------------------------------------------*/
	static well_path_unroll()
};

}} /* blue_sky::wpi */

#endif /* end of include guard: WPI_ALGO_UNROLL_YAUH9FX2 */

