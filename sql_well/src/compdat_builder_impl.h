/// @file compdat_builder.h
/// @brief COMPDAT builder implementation
/// @author uentity
/// @version 
/// @date 12.12.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef COMPDAT_BUILDER_IMPL_RXNBKXX0
#define COMPDAT_BUILDER_IMPL_RXNBKXX0

#include "frac_comp_ident.h"

namespace blue_sky { 
namespace fci {

/*-----------------------------------------------------------------
 * compdat_builder::impl
 *----------------------------------------------------------------*/
class compdat_builder::impl {
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
			if (xp.empty ()) // no intersections
				continue;

			// 3.3 select all completions that belong to well+branch_i
			q = (cd_traits::select_segment() % date % pwb->first % pwb->second).str();
			sw_->prepare_sql(q);

			// 3.4 for all completions do
			while(sw_->step_sql() == 0) {
				// 3.4.1 search for completion_j begin_j and end_j using md_j and lentgh_j
				t_double md = sw_->get_sql_real (0);
				t_double len = sw_->get_sql_real (1);
				//t_uint status = sw_->get_sql_int (2);
				//t_double rw = sw_->get_sql_real (3);
				//t_double kh = sw_->get_sql_real (4);
				//t_double skin = sw_->get_sql_real (5);
				x_iterator px = xp.upper_bound(whc(md));
				x_iterator xend = xp.upper_bound(whc(md + len));
				x_iterator pprev_x = px;
				// always start with prev intersection
				if(px != xp.begin())
					--pprev_x;

				// add last segment to intersection loop
				if (xend != xp.end ())
					++xend;

				// 3.4.2 consider all intersections between begin_j and end_j
				for(; px != xend; ++px) {
					// prepare compdat
					compdat cf(pwb->first, pwb->second, px->cell, m_size_);

					// 3.4.3.1 calc delta between consequent xpoint_k and xpoint_(k - 1)
					// position to previous point
					double delta_l = 0;
					if(px != xp.end()) {
						// completion fully inside well segment
						if (md >= pprev_x->md && (md + len) <= px->md) {
							delta_l = len;
							cf.md = md;
							cf.len = delta_l;
							for (t_uint j = 0; j < strat_t::D; ++j) {
								cf.x1[j] = pprev_x->where[j] + (md - pprev_x->md) / (px->md - pprev_x->md) * (px->where[j] - pprev_x->where[j]);
								cf.x2[j] = pprev_x->where[j] + (md + len - pprev_x->md) / (px->md - pprev_x->md) * (px->where[j] - pprev_x->where[j]);
							}
						}
						// well segment fully inside completion
						else if (md <= pprev_x->md && (md + len) >= px->md) {
							delta_l = px->md - pprev_x->md;
							cf.md = pprev_x->md;
							cf.len = delta_l;
							for (t_uint j = 0; j < strat_t::D; ++j) {
								cf.x1[j] = pprev_x->where[j];
								cf.x2[j] = px->where[j];
							}
						}
						// start of completion is inside well segment
						// end of completion is out of well segment
						else if (md >= pprev_x->md && (md + len) >= px->md) {
							delta_l = px->md - md;
							cf.md = md;
							cf.len = delta_l;
							for (t_uint j = 0; j < strat_t::D; ++j) {
								cf.x1[j] = pprev_x->where[j] + (md - pprev_x->md) / (px->md - pprev_x->md) * (px->where[j] - pprev_x->where[j]);
								cf.x2[j] = px->where[j];
							}
						}
						// start of completion is outside well segment
						// end of completion is inside of well segment
						else if (md <= pprev_x->md && (md + len) <= px->md) {
							delta_l = md + len - pprev_x->md;
							cf.md = pprev_x->md;
							cf.len = delta_l;
							for (t_uint j = 0; j < strat_t::D; ++j) {
								cf.x1[j] = pprev_x->where[j];
								cf.x2[j] = pprev_x->where[j] + (md + len - pprev_x->md) / (px->md - pprev_x->md) * (px->where[j] - pprev_x->where[j]);
							}
						}
					}
					else {
						// TODO: identify coordinates of last cell intersection???
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

						//if(delta == 1) {
						//	cf.dir = 'X';
						//	cf.kh_mult = delta_l / cell_sz[0];
						//}
						//else if(delta >= m_size_[0] && delta < plane_sz) {
						//	cf.dir = 'Y';
						//	cf.kh_mult = delta_l / cell_sz[1];
						//}
						//else {
						//	cf.dir = 'Z';
						//	cf.kh_mult = delta_l / cell_sz[2];
						//}
						cf.kh_mult = std::min(cf.kh_mult, 1.);

						// set data
						//cf.diam = 2 * rw;
						//cf.kh = kh;
						//cf.status = status;
						//cf.skin = skin;

						// if compdat for this cell is already added
						// then just update kh_mult
						// otherwise add new COMPDAT record
						// TODO: handle case of different directions inside one cell
						cd_storage::iterator pcd = cfs_.find(compdat(px->cell));
						if(pcd != cfs_.end()) {
							compdat& cur_cd = const_cast< compdat& >(*pcd);
							cur_cd.kh_mult = std::min(cur_cd.kh_mult + cf.kh_mult, 1.);
						}
						else
							cfs_.insert(cf);
					}
          // set next element
          pprev_x = px;
				} // 3.4.4 end of intersections loop
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
	cd_storage cfs_;
};

}}

#endif /* end of include guard: COMPDAT_BUILDER_RXNBKXX0 */

