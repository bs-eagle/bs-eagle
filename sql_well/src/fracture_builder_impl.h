/// @file fracture_builder_impl.h
/// @brief FRACTURE builder implementation
/// @author uentity
/// @version 
/// @date 12.12.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef FRACTURE_BUILDER_IMPL_RAF5LK1T
#define FRACTURE_BUILDER_IMPL_RAF5LK1T

#include "frac_comp_ident.h"

namespace blue_sky { namespace fci {

/*-----------------------------------------------------------------
 * fracture_builder::impl
 *----------------------------------------------------------------*/
class fracture_builder::impl {
public:
	typedef strategy_3d strat_t;
	typedef wpi_algo::trimesh trimesh;
	typedef wpi_algo::well_path well_path;
	typedef wpi_algo::well_hit_cell whc;
	typedef wpi_algo::intersect_path xpath;
	typedef wpi_algo::hit_idx_t hit_idx_t;
	typedef wpi_algo::xbuilder xbuilder;
	//typedef intersect_builder2< strategy_3d > xbuilder;

	typedef strategy_3d::vertex_pos_i vertex_pos_i;
	typedef strategy_3d::vertex_pos vertex_pos;
	typedef mesh_tools< strat_t >::mesh_part mesh_part;

	typedef sql_well::sp_traj_t sp_traj_t;
	typedef sql_well::sp_table_t sp_table_t;

	typedef multimap< string, string > wb_storage;
	typedef xpath::iterator x_iterator;

	impl(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn) {
		init(nx, ny, coord, zcorn);
	}

	impl(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn,
		smart_ptr< well_pool_iface, true > src_well)
	{
		init(nx, ny, coord, zcorn);
		init(src_well);
	}

	void init(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn) {
		// build trimesh for given coord+zcorn
		tops_ = wpi_algo::coord_zcorn2trimesh(nx, ny, coord, zcorn, m_, m_size_);
	}

	void init(smart_ptr< well_pool_iface, true > src_well) {
		sw_ = BS_KERNEL.create_object_copy(src_well);
	}

	template< class cd_traits >
	void build(double date, const cd_traits& t = cd_traits()) {
		//cfs_.clear();

		// 1 fill storage with all unique well+branch
		wb_storage wb;
		std::string q = (cd_traits::select_unique_well_branch() % date).str();
		sw_->prepare_sql(q);
		while(sw_->step_sql() == 0) {
			wb.insert(make_pair(sw_->get_sql_str(0), sw_->get_sql_str(1)));
		}
		sw_->finalize_sql();

		// 2 precalc plane size
		//const ulong plane_sz = m_size_[0] * m_size_[1];
		// prepare mesh_part representin full mesh
		mesh_part fullmesh(m_, m_size_);

		// 3 for each well+branch combo do
		for(wb_storage::iterator pwb = wb.begin(), wb_end = wb.end(); pwb != wb_end; ++pwb) {
			// 3.1 from 'branches' select well+branch_i trajectory (sql_well::get_branch_traj)
			sp_traj_t traj = sw_->get_branch_traj(pwb->first, pwb->second);
			if(!traj) return;
			sp_table_t traj_t = traj->get_table();
			if(!traj_t) return;

			// find column id's for X, Y, Z and MD
			t_long col_ids[4] = {1, 2, 3, 0};
			//for(t_long i = 0; i < traj_t->get_n_cols(); ++i) {
			//	std::string cur_col = traj_t->get_col_name(i);
			//	if(cur_col == "X")
			//		col_ids[0] = i;
			//	else if(cur_col == "Y")
			//		col_ids[1] = i;
			//	else if(cur_col == "Z")
			//		col_ids[2] = i;
			//	else if(cur_col == "MD")
			//		col_ids[3] = i;
			//}

			// fill vector with traj data
			spv_float traj_v = BS_KERNEL.create_object(v_float::bs_type());
			traj_v->resize(traj_t->get_n_rows() * 4);
			v_float::iterator ptv = traj_v->begin();
			for(ulong i = 0, trows = traj_t->get_n_rows(); i < trows; ++i) {
				for(ulong j = 0; j < 4; ++j)
					*ptv++ = traj_t->get_value(i, col_ids[j]);
			}

			// make well_path
			well_path W;
			if(!wpi_algo::fill_well_path(traj_v, W)) return;

			// 3.2 find intersections of given branch with mesh (well_path_ident)
			xbuilder A(m_, W, m_size_);
			A.build();
			A.remove_dups2();
			//A.append_wp_nodes(hi);
			xpath& xp = A.path();
			// check for empty intersections
			if (xp.empty ())
				continue;

			// 3.3 select all completions that belong to well+branch_i
			q = (cd_traits::select_segment() % date % pwb->first % pwb->second).str();
			sw_->prepare_sql(q);

			// 3.4 for all completions do
			while(sw_->step_sql() == 0) {
				// 3.4.1 search for completion_j begin_j and end_j using md_j and lentgh_j
				t_int    status        = sw_->get_sql_int  (0);
				t_double md            = sw_->get_sql_real (1);
				t_double half_length_1 = sw_->get_sql_real (2);
				t_double half_length_2 = sw_->get_sql_real (3);
				t_double angle         = sw_->get_sql_real (4);
				t_double half_up       = sw_->get_sql_real (5);
				t_double half_down     = sw_->get_sql_real (6);
				t_double perm          = sw_->get_sql_real (7);
				t_double half_thin     = sw_->get_sql_real (8);

				x_iterator px = xp.upper_bound(whc(md));

				if (px == xp.end ())
					continue;
				// position to next point
				x_iterator pprev_x = px;
				// always start with prev intersection
				if(px != xp.begin())
					--pprev_x;

				// prepare compdat
				fracture frac(pwb->first, pwb->second, px->cell, m_size_);

				vertex_pos   frac_coords;
				vertex_pos_i cell_pos;
				wpi_algo::decode_cell_id (px->cell, cell_pos, m_size_);

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

				for (t_uint j = 0; j < strat_t::D; ++j) {
					frac_coords[j] = pprev_x->where[j] + (md - pprev_x->md) / (px->md - pprev_x->md) * (px->where[j] - pprev_x->where[j]);
				}

				t_uint kw1_flag = 0;
				// find kw1 position of fracture
				for (t_long kw = cell_pos[2]; kw >= 0; kw--) {
					vertex_pos_i cell_up;
					copy (&cell_pos[0], &cell_pos[strat_t::D], &cell_up[0]);
					cell_up[2] = kw;

					t_ulong k_cell_id_ = wpi_algo::encode_cell_id (cell_up, m_size_);
					t_double z_top = 0;
					t_float* tops_data = &(*tops_)[0];

					// get mean of Z top plane
					for (t_uint j = 0; j < 4; ++j) {
						z_top += 0.25 * tops_data[24 * k_cell_id_ + 3 * j + 2];
					}
					// fracture inside this layer
					if (frac_coords[2] - half_up > z_top) {
						kw1_flag = 1;
						frac.cell_pos[2] = kw;  // Kw1 - position
						break;
					}
				}

				if (kw1_flag == 0) {
					frac.md_cell_pos[2] = 0;
				}

				t_uint kw2_flag = 0;
				// find kw2 position of fracture
				for (t_ulong kw = cell_pos[2]; kw < m_size_[2]; kw++) {
					vertex_pos_i cell_down;
					copy (&cell_pos[0], &cell_pos[strat_t::D], &cell_down[0]);
					cell_down[2] = kw;

					t_ulong k_cell_id_ = wpi_algo::encode_cell_id (cell_down, m_size_);
					t_double z_down = 0;
					t_float* tops_data = &(*tops_)[0];

					// get mean of Z bottom plane
					for (t_uint j = 4; j < 8; ++j) {
						z_down += 0.25 * tops_data[24 * k_cell_id_ + 3 * j + 2];
					}
					// fracture inside this layer
					if (frac_coords[2] + half_down < z_down) {
						kw2_flag = 1;
						frac.cell_pos[3] = kw; // Kw2 position
						break;
					}
				}

				if (kw2_flag == 0) {
					frac.md_cell_pos[3] = m_size_[2] - 1;
				}

				// if compdat for this cell is already added
				// then just update kh_mult
				// otherwise add new COMPDAT record
				// TODO: handle case of different directions inside one cell
				frac_storage::iterator pcd = cfs_.find(fracture(px->cell));
				if(pcd != cfs_.end()) {
					//fracture& cur_cd = const_cast< fracture& >(*pcd);
					// TODO: cur_cd.kh_mult = std::min(cur_cd.kh_mult + cf.kh_mult, 1.);
				}
				else
					cfs_.insert(frac);
			} // 3.5 end of completions loop

			sw_->finalize_sql();
		} // 4 end of well+branch loop
	}

	//virtual std::string select_unique_well_branch(double d) = 0;
	//virtual std::string select_segment(double d, const std::string well_name,
	//	const std::string& branch_name) = 0;

	trimesh m_;
	vertex_pos_i m_size_;
	// tops should live as long as mesh lives
	spv_float tops_;
	// copy of source sql_well
	smart_ptr< well_pool_iface, true > sw_;
	// storage for compdat records
	frac_storage cfs_;
};

}}

#endif /* end of include guard: FRACTURE_BUILDER_IMPL_RAF5LK1T */

