/// @file fracture_traits.h
/// @brief FRACTURE traits implementation
/// @author uentity
/// @version 
/// @date 12.12.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef FRACTURE_TRAITS_RAF5LK1T
#define FRACTURE_TRAITS_RAF5LK1T

#include "frac_comp_builder.h"

namespace blue_sky { namespace fci {

struct fract_traits  {
	typedef strategy_3d strat_t;
	typedef strat_t::vertex_pos_i vertex_pos_i;
	typedef strat_t::vertex_pos vertex_pos;
	typedef frac_storage storage_t;
	typedef wpi_algo::intersect_path xpath;
	typedef xpath::iterator x_iterator;
	typedef wpi_algo::well_hit_cell whc;

	typedef frac_comp_builder< fract_traits > fcb_t;

	static boost::format select_unique_well_branch() {
		return boost::format("SELECT DISTINCT well_name, branch_name FROM fractures WHERE d=%f");
	}

	static boost::format select_segment() {
		return boost::format(
			"SELECT status, md, half_length_1, half_length_2, angle, half_up, half_down, perm,\
			half_thin, skin FROM fractures WHERE d=%f and well_name='%s' and branch_name='%s'\
			ORDER BY md"
			//"SELECT status, md, half_length_1, half_length_2, angle, half_up, half_down, perm,\
			//half_thin, skin FROM fractures WHERE d=%f and well_name='%s' and branch_name='%s'"
		);
	}

	static void proc_length(
			const std::string& well_name, const std::string& branch_name,
			xpath& xp, fcb_t& fcb
	) {
		typedef fcb_t::sp_srcwell sp_srcwell;
		sp_srcwell sw = fcb.sw_;

		// 3.4.1 search for completion_j begin_j and end_j using md_j and lentgh_j
		t_int    status        = sw->get_sql_int  (0);
		t_double md            = sw->get_sql_real (1);
		t_double half_length_1 = sw->get_sql_real (2);
		t_double half_length_2 = sw->get_sql_real (3);
		t_double angle         = sw->get_sql_real (4);
		t_double half_up       = sw->get_sql_real (5);
		t_double half_down     = sw->get_sql_real (6);
		t_double perm          = sw->get_sql_real (7);
		t_double half_thin     = sw->get_sql_real (8);
		t_double frac_skin     = sw->get_sql_real (9);
		x_iterator px = xp.upper_bound(whc(md));

		if (px == xp.end ())
			return;
		// position to next point
		x_iterator pprev_x = px;
		// always start with prev intersection
		if(px != xp.begin())
			--pprev_x;

		// prepare compdat
		fracture frac(well_name, branch_name, pprev_x->cell, fcb.m_size_);

		vertex_pos   frac_coords;
		vertex_pos_i cell_pos;
		wpi_algo::decode_cell_id (pprev_x->cell, cell_pos, fcb.m_size_);

		frac.md_cell_pos[0] = frac.cell_pos[0] = cell_pos[0];   // I - position
		frac.md_cell_pos[1] = frac.cell_pos[1] = cell_pos[1];   // J - position
		if (strat_t::D == 3)
			frac.md_cell_pos[2] = cell_pos[2];
		frac.frac_half_length_1 = half_length_1;
		frac.frac_half_length_2 = half_length_2;
		frac.frac_status = status;
		frac.frac_angle = angle;
		frac.frac_perm  = perm;
		frac.frac_half_thin = half_thin;
		frac.frac_skin = frac_skin;

		for (t_uint j = 0; j < strat_t::D; ++j) {
			frac_coords[j] = pprev_x->where[j] + (md - pprev_x->md) /
				(px->md - pprev_x->md) * (px->where[j] - pprev_x->where[j]);
		}

		//std::cout<<"up "<<half_up<<" down "<<half_down<<"\n";


		t_uint kw1_flag = 0;
		// find kw1 position of fracture
		for (t_long kw = cell_pos[2]; kw >= 0; kw--) {
			//std::cout<<"kw "<<kw<<" 0";

			vertex_pos_i cell_up;
			copy (&cell_pos[0], &cell_pos[strat_t::D], &cell_up[0]);
			cell_up[2] = kw;

			t_long k_cell_id_ = wpi_algo::encode_cell_id (cell_up, fcb.m_size_);
			t_double z_top = 0;
			t_float* tops_data = &(*fcb.tops_)[0];

			// get mean of Z top plane
			for (t_uint j = 0; j < 4; ++j) {
				z_top += 0.25 * tops_data[24 * k_cell_id_ + 3 * j + 2];
			}
			//std::cout<<" frac "<<frac_coords[2] - half_up<<" z_down "<<z_top<<"\n";

			// fracture inside this layer
			if (frac_coords[2] - half_up >= z_top) {
				kw1_flag = 1;
				frac.cell_pos[2] = kw;  // Kw1 - position
				//std::cout<<"    break\n";
				break;
			}
		}

		if (kw1_flag == 0)  {
			// cell_pos[2] hasn't changed
			// reach top of reservoir
			frac.cell_pos[2] = 0;
		}

		t_uint kw2_flag = 0;
		// find kw2 position of fracture
		for (t_ulong kw = cell_pos[2]; kw < fcb.m_size_[2]; kw++) {
			//std::cout<<"kw "<<kw<<"m_size "<<fcb.m_size_[2];

			vertex_pos_i cell_down;
			copy (&cell_pos[0], &cell_pos[strat_t::D], &cell_down[0]);
			cell_down[2] = kw;

			t_ulong k_cell_id_ = wpi_algo::encode_cell_id (cell_down, fcb.m_size_);
			t_double z_down = 0;
			t_float* tops_data = &(*fcb.tops_)[0];

			// get mean of Z bottom plane
			for (t_uint j = 4; j < 8; ++j) {
				z_down += 0.25 * tops_data[24 * k_cell_id_ + 3 * j + 2];
			}

			//std::cout<<" frac "<<frac_coords[2] + half_down<<" z_down "<<z_down<<"\n";

			// fracture inside this layer
			if (frac_coords[2] + half_down <= z_down) {
				kw2_flag = 1;
				frac.cell_pos[3] = kw; // Kw2 position
				//std::cout<<"    break\n";
				break;
			}

		}
		if (kw2_flag == 0) {
			// reach bottom of reservoir
			frac.cell_pos[3] = fcb.m_size_[2] - 1;
		}

		if (!half_down)
			frac.cell_pos[3] = frac.cell_pos[2];

		if (!half_up && !half_down)
			return;

		//std::cout<<"FRACTURE K1 "<<frac.cell_pos[2]<<" K2 "<<frac.cell_pos[3]<<"\n";

		// if compdat for this cell is already added
		// then just update kh_mult
		// otherwise add new COMPDAT record
		// TODO: handle case of different directions inside one cell
		frac_storage::iterator pcd = fcb.s_.find(fracture(pprev_x->cell));
		if(pcd != fcb.s_.end()) {
			//fracture& cur_cd = const_cast< fracture& >(*pcd);
			// TODO: cur_cd.kh_mult = std::min(cur_cd.kh_mult + cf.kh_mult, 1.);
		}
		else
			fcb.s_.insert(frac);
	}
};

}}

#endif /* end of include guard: FRACTURE_BUILDER_IMPL_RAF5LK1T */

