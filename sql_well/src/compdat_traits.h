/// @file compdat_traits.h
/// @brief COMPDAT traits implementation
/// @author uentity
/// @version 
/// @date 12.12.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef COMPDAT_TRAITS_RXNBKXX0
#define COMPDAT_TRAITS_RXNBKXX0

#include "frac_comp_builder.h"

namespace blue_sky { namespace fci {

struct compl_traits {
	typedef strategy_3d strat_t;
	typedef strat_t::vertex_pos_i vertex_pos_i;
	typedef strat_t::vertex_pos vertex_pos;
	typedef wpi_algo::intersect_path xpath;
	typedef cd_storage storage_t;
	typedef xpath::iterator x_iterator;
	typedef wpi_algo::well_hit_cell whc;
	typedef wpi_algo::mesh_part mesh_part;

	typedef frac_comp_builder< compl_traits > fcb_t;

	static boost::format select_unique_well_branch() {
		return boost::format("SELECT DISTINCT well_name, branch_name FROM completions WHERE d=%f");
	}

	static boost::format select_segment() {
		return boost::format(
			//"SELECT md, length FROM completions WHERE d=%f and well_name='%s' and branch_name='%s'"
			"SELECT md, length, status, rw, kh, skin, kh_mult FROM\
			completions WHERE d=%f and well_name='%s' and branch_name='%s' ORDER BY md"
		);
	}

	static void proc_length(
			const std::string& well_name, const std::string& branch_name,
			xpath& xp, fcb_t& fcb
	){
		typedef fcb_t::sp_srcwell sp_srcwell;
		sp_srcwell sw = fcb.sw_;

		// prepare mesh_part representing full mesh
		mesh_part fullmesh(fcb.m_, fcb.m_size_);

		// 3.4.1 search for completion_j begin_j and end_j using md_j and lentgh_j
		t_double md      = sw->get_sql_real (0);
		t_double len     = sw->get_sql_real (1);
		t_uint status    = sw->get_sql_int (2);
		t_double rw      = sw->get_sql_real (3);
		t_double kh      = sw->get_sql_real (4);
		t_double skin    = sw->get_sql_real (5);
		t_double kh_mult = sw->get_sql_real (6);

		x_iterator px      = xp.upper_bound(whc(md));
		x_iterator xend    = xp.upper_bound(whc(md + len));
		// sanity check if completion inside well segment
		if(px == xp.end()) return;

		// always start with prev intersection
		x_iterator pnext_x = px;
		if(px != xp.begin())
			--px;
		else
			++pnext_x;

		// 3.4.2 consider all intersections between begin_j and end_j
		for(; px != xend; ++px) {
			// prepare compdat
			compdat cf(well_name, branch_name, px->cell, fcb.m_size_);

			// 3.4.3.1 calc delta between consequent xpoint_k and xpoint_(k - 1)
			// position to previous point
			double delta_l = 0;
			//if(px != xp.end()) {
			// completion fully inside well segment
			if (md >= px->md && (md + len) <= pnext_x->md) {
				delta_l = len;
				cf.md = md;
				cf.len = delta_l;
				for (t_uint j = 0; j < strat_t::D; ++j) {
					cf.x1[j] = px->where[j] + (md - px->md) /
						(pnext_x->md - px->md) * (pnext_x->where[j] - px->where[j]);
					cf.x2[j] = px->where[j] + (md + len - px->md) /
						(pnext_x->md - px->md) * (pnext_x->where[j] - px->where[j]);
				}
			}
			// well segment fully inside completion
			else if (md <= px->md && (md + len) >= pnext_x->md) {
				delta_l = pnext_x->md - px->md;
				cf.md = px->md;
				cf.len = delta_l;
				for (t_uint j = 0; j < strat_t::D; ++j) {
					cf.x1[j] = px->where[j];
					cf.x2[j] = pnext_x->where[j];
				}
			}
			// start of completion is inside well segment
			// end of completion is out of well segment
			else if (md >= px->md && (md + len) >= pnext_x->md) {
				delta_l = pnext_x->md - md;
				cf.md = md;
				cf.len = delta_l;
				for (t_uint j = 0; j < strat_t::D; ++j) {
					cf.x1[j] = px->where[j] + (md - px->md) /
						(pnext_x->md - px->md) * (pnext_x->where[j] - px->where[j]);
					cf.x2[j] = pnext_x->where[j];
				}
			}
			// start of completion is outside well segment
			// end of completion is inside of well segment
			else if (md <= px->md && (md + len) <= pnext_x->md) {
				delta_l = md + len - px->md;
				cf.md = px->md;
				cf.len = delta_l;
				for (t_uint j = 0; j < strat_t::D; ++j) {
					cf.x1[j] = px->where[j];
					cf.x2[j] = px->where[j] + (md + len - px->md) /
						(pnext_x->md - px->md) * (pnext_x->where[j] - px->where[j]);
				}
			}
			// set next element
			if(++pnext_x == xp.end()) {
				pnext_x = px;
				++pnext_x;
			}

			// 3.4.3.2 if delta == 1 mark direction as 'X'
			//         else if delta == dx direction = 'Y'
			//         else if delta = dx*dy direction = 'Z'
			//         also calc kh_mult assuming that cells are rectangular (!)
			// update: use better recognition of direction using max coord change
			//         between two intersections
			if(delta_l) {
				// obtain cell size
				vertex_pos cell_sz = {1., 1., 1.};
				fullmesh.cell_size(px->cell, cell_sz);

				// select direction, calc kh increase
				const char dirs[] = {'X', 'Y', 'Z'};
				double max_step = 0;
				for(uint i = 0; i < strat_t::D; ++i) {
					double dir_step = std::fabs(cf.x2[i] - cf.x1[i]);
					if (max_step < dir_step) {
						max_step = dir_step;
						cf.dir = dirs[i];
						cf.kh_mult = delta_l / cell_sz[i];
					}
				}
				cf.kh_mult = std::min(cf.kh_mult, 1.);

				// set data
				cf.diam      = 2 * rw;
				cf.kh        = kh;
				cf.status    = status;
				cf.skin      = skin;
				cf.kh_mult  *= kh_mult;


				// if compdat for this cell is already added
				// then just update kh_mult
				// otherwise add new COMPDAT record
				// TODO: handle case of different directions inside one cell
				cd_storage::iterator pcd = fcb.s_.find(compdat(px->cell));
				if(pcd != fcb.s_.end()) {
					compdat& cur_cd = const_cast< compdat& >(*pcd);
					cur_cd.kh_mult = std::min(cur_cd.kh_mult + cf.kh_mult, 1.);
				}
				else
					fcb.s_.insert(cf);
			}
		} // 3.4.4 end of intersections loop
	}
};

}}

#endif /* end of include guard: COMPDAT_BUILDER_RXNBKXX0 */

