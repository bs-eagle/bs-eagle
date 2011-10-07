/// @file frac_comp_ient.cpp
/// @brief Implementation of finding fractions/completions intersections with mesh
/// @author uentity
/// @version 
/// @date 06.10.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "frac_comp_ident.h"
#include "bcsr_matrix_iface.h"
#include "wpi_algo.h"
#include "wpi_strategy_3d.h"

#include <boost/format.hpp>

using namespace wpi;

namespace blue_sky { namespace fci {
// impl details
namespace {
using namespace std;

// assign for c arrays
// fun with returning reference to array :)
template< class T, unsigned int D >
static T (&ca_assign(T (&lhs)[D], const T (&rhs)[D]))[D] {
	std::copy(&rhs[0], &rhs[D], &lhs[0]);
	return lhs;
}

template< class T, unsigned int D >
static T (&ca_assign(T (&lhs)[D], const T& v))[D] {
	std::fill(&lhs[0], &lhs[D], v);
	return lhs;
}

// POD that holds all info needed by COMPDAT
struct compfrac {
	typedef stat_array< ulong, 4 > cell_info;
	typedef strategy_3d::vertex_pos_i vertex_pos_i;

	std::string well_name;
	std::string branch_name;
	cell_info cell_id;
	char dir;

	compfrac(const string& well_name_, const string& branch_name_, const cell_info& cell_info, char dir_)
		: well_name(well_name_), branch_name(branch_name_), cell_id(cell_info), dir(dir_)
	{}

	compfrac(const string& well_name_, const string& branch_name_, const vertex_pos_i& id)
		: well_name(well_name_), branch_name(branch_name_)
	{
		copy(&id[0], &id[strategy_3d::D], cell_id.begin());
		cell_id[3] = cell_id[2];
	}
};

typedef std::list< compfrac > cf_storage;

} /* hidden namespace*/

spv_float completions_ident(smart_ptr< sql_well > sw, double date,
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn)
{
	typedef algo< strategy_3d > wpi_algo;
	typedef typename wpi_algo::trimesh trimesh;
	typedef typename wpi_algo::well_path well_path;
	typedef typename wpi_algo::well_hit_cell whc;
	typedef typename wpi_algo::intersect_path xpath;
	typedef typename wpi_algo::xbuilder xbuilder;
	typedef typename wpi_algo::hit_idx_t hit_idx_t;
	typedef typename wpi_algo::vertex_pos_i vertex_pos_i;

	typedef typename sql_well::sp_traj_t sp_traj_t;
	typedef typename sql_well::sp_table_t sp_table_t;

	typedef multimap< string, string > wb_storage;
	typedef typename xpath::iterator x_iterator;

	// 1 find all unique well+branch names on given date
	wb_storage wb;
	string q = (boost::format(
		"SELECT DISTINCT well_name, branch_name FROM completions WHERE d=%f") % date).str();
	sw->prepare_sql(q);
	// fill storage with all unique well+branch
	while(sw->step_sql() == 0) {
		wb.insert(make_pair(sw->get_sql_str(0), sw->get_sql_str(1)));
	}
	sw->finalize_sql();

	// 2 build trimesh for given coord+zcorn
	trimesh M;
	vertex_pos_i mesh_size;
	spv_float tops = wpi_algo::coord_zcorn2trimesh(nx, ny, coord, zcorn, M, mesh_size, true);
	const ulong plane_sz = mesh_size[0] * mesh_size[1];

	// 3 for each well+branch combo do
	cf_storage cfs;
	for(wb_storage::iterator pwb = wb.begin(), wb_end = wb.end(); pwb != wb_end; ++pwb) {
		// 3.1 from 'branches' select well+branch_i trajectory (sql_well::get_branch_traj)
		sp_traj_t traj = sw->get_branch_traj(pwb->first, pwb->second);
		if(!traj) return NULL;
		sp_table_t traj_t = traj->get_table();
		if(!traj_t) return NULL;

		// 3.2 find intersections of given branch with mesh (well_path_ident)
		// fill array with branch trajectory
		spv_float traj_v = BS_KERNEL.create_object(v_float::bs_type());
		traj_v->resize(traj_t->get_n_rows() * 4);
		v_float::iterator ptv = traj_v->begin();
		for(ulong i = 0, trows = traj_t->get_n_rows(); i < trows; ++i) {
			for(ulong j = 0; j < 4; ++j)
				*ptv++ = traj_t->get_value(i, 0);
		}
		// make well_path
		well_path W;
		if(!wpi_algo::fill_well_path(traj_v, W)) return NULL;
		// find intersections
		xbuilder A(M, W, mesh_size);
		A.build();
		A.remove_dups2();
		//A.append_wp_nodes(hit_idx);

		// 3.3 select all completions that belong to well+branch_i
		q = (boost::format(
			"SELECT md, length FROM completions WHERE well_name='%1%' and well_branch='%2%'")
			% pwb->first % pwb->second).str();
		sw->prepare_sql(q);

		// 3.4 for all completions do
		while(sw->step_sql() == 0) {
			// 3.4.1 search for completion_j begin_j and end_j using md_j and lentgh_j
			xpath& xp = A.path();
			x_iterator px = xp.upper_bound(whc(sw->get_sql_real(0)));
			x_iterator xend = xp.upper_bound(whc(sw->get_sql_real(0) + sw->get_sql_real(1)));

			// 3.4.2 consider all intersections between begin_j and end_j
			for(; px != xend; ++px) {
				// prepare compfrac
				vertex_pos_i cell_id;
				wpi_algo::decode_cell_id(px->cell, cell_id, mesh_size);
				compfrac cf(pwb->first, pwb->second, cell_id);

				// 3.4.3.1 calc delta between consequent xpoint_k and xpoint_(k + 1)
				// position to prev point
				x_iterator pprev_x = px;
				if(pprev_x != xp.begin())
					--pprev_x;
				ulong delta = std::abs(px->cell - pprev_x->cell);

				// 3.4.3.2 if delta == 1 mark direction as 'X'
				//         else if delta == dx direction = 'Y'
				//         else if delta = dx*dy direction = 'Z'

				if(delta == 1)
					cf.dir = 'X';
				else if(delta >= mesh_size[0] && delta < plane_sz)
					cf.dir = 'Y';
				else
					cf.dir = 'Z';
				// add new COMPDAT record
				cfs.push_back(cf);
			} // 3.4.4 end of intersections loop
		} // 3.5 end of completions loop
	} // 4 end of well+branch loop

	return NULL;
}

}} /* blue_sky::fci */
