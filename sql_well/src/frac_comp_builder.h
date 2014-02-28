/// @file frac_comp_builder.h
/// @brief Implementation of common COMPDAT/FRACTURE building algorithm
/// @author uentity
/// @version 
/// @date 12.12.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef FRAC_COMP_BUILDER_NIGNH6R8
#define FRAC_COMP_BUILDER_NIGNH6R8

#include "frac_comp_ident.h"

namespace blue_sky { namespace fci {

/*-----------------------------------------------------------------
 * compdat_builder::impl
 *----------------------------------------------------------------*/
template< class cd_traits >
class frac_comp_builder {
public:
	typedef wpi_algo::trimesh         trimesh;
	typedef wpi_algo::well_path       well_path;
	typedef wpi_algo::well_paths      well_paths;
	typedef wpi_algo::well_hit_cell   whc;
	typedef wpi_algo::intersect_path  xpath;
	typedef wpi_algo::intersect_paths xpaths;
	typedef wpi_algo::hit_idx_t       hit_idx_t;
	typedef wpi_algo::xbuilder        xbuilder;
	typedef wpi_algo::xbuilder_mp     xbuilder_mp;

	typedef strat_t::vertex_pos_i            vertex_pos_i;
	typedef strat_t::vertex_pos              vertex_pos;
	typedef mesh_tools< strat_t >::mesh_part mesh_part;

	typedef sql_well::sp_traj_t sp_traj_t;
	typedef sql_well::sp_table_t sp_table_t;

	typedef multimap< string, string > wb_storage;
	typedef xpath::iterator x_iterator;

	typedef typename cd_traits::storage_t storage_t;
	typedef smart_ptr< well_pool_iface, true > sp_srcwell;

	typedef std::map< std::string, xpath > xp_cache_t;
	typedef st_smart_ptr< xp_cache_t >     spxp_cache_t;
	typedef typename xp_cache_t::iterator  xpc_iterator;

	frac_comp_builder() {}

	frac_comp_builder(t_ulong nx, t_ulong ny, const spv_float& coord, const spv_float& zcorn) {
		init(nx, ny, coord, zcorn);
	}

	frac_comp_builder(t_ulong nx, t_ulong ny, const spv_float& coord, const spv_float& zcorn,
		smart_ptr< well_pool_iface, true > src_well)
	{
		init(nx, ny, coord, zcorn);
		init(src_well);
	}

	void init(t_ulong nx, t_ulong ny, const spv_float& coord, const spv_float& zcorn) {
		// build trimesh for given coord+zcorn
		m_.init(nx, ny, coord, zcorn);
	}

	void init(t_long nx, t_long ny, const sp_obj& trim_backend) {
		m_.init(nx, ny, trim_backend);
	}

	void init(smart_ptr< well_pool_iface, true > src_well) {
		sw_ = src_well;
		//sw_ = BS_KERNEL.create_object_copy(src_well);
	}

	void init_cache(const spxp_cache_t& xp_cache = NULL, const ulong cache_limit = 0) {
		if(xp_cache)
			xp_cache_ = xp_cache;
		else
			xp_cache_ = new xp_cache_t;
		cache_limit_ = cache_limit;
	}

	// access to stored cache
	spxp_cache_t cache() const {
		return xp_cache_;
	}

	bool build_cache() {
		// 1 fill storage with all unique well+branch
		if(!sw_) return false;
		wb_storage wb;
		// select all wells+branch in completions and fractures tables
		std::string q =
			"SELECT DISTINCT well_name, branch_name FROM branches";
			//"SELECT well_name, branch_name FROM completions UNION SELECT well_name, branch_name FROM fractures";
		//std::string q = (cd_traits::select_unique_well_branch() % date).str();
		sw_->prepare_sql(q);
		while(sw_->step_sql() == 0) {
			wb.insert(make_pair(sw_->get_sql_str(0), sw_->get_sql_str(1)));
		}
		sw_->finalize_sql();

		// column id's for X, Y, Z and MD
		const t_long col_ids[4] = {1, 2, 3, 0};

		// 2. fill well_paths vector
		well_paths W(wb.size());
		// well_paths reference source traj data, so keep it while algo works
		std::vector< spv_float > trajes(wb.size());
		ulong k = 0;
		// for each well+branch combo build well_infos
		for(wb_storage::iterator pwb = wb.begin(), wb_end = wb.end(); pwb != wb_end; ++pwb, ++k) {
			// from 'branches' select well+branch_i trajectory (sql_well::get_branch_traj)
			sp_traj_t traj = sw_->get_branch_traj(pwb->first, pwb->second);
			if(!traj) return false;
			sp_table_t traj_t = traj->get_table();
			if(!traj_t) return false;

			// fill vector with traj data
			spv_float traj_v = BS_KERNEL.create_object(v_float::bs_type());
			traj_v->resize(traj_t->get_n_rows() * 4);
			v_float::iterator ptv = traj_v->begin();
			for(ulong i = 0, trows = traj_t->get_n_rows(); i < trows; ++i) {
				for(ulong j = 0; j < 4; ++j)
					*ptv++ = traj_t->get_value(i, col_ids[j]);
			}
			// prevent deleting source traj data
			trajes[k] = traj_v;

			// convert traj_v to k-th well_path
			if(!wpi_algo::fill_well_path(traj_v, W[k])) return false;
		}

		// 3 find intersections of all branches with mesh (well_paths_ident)
		xbuilder_mp A(m_, W);
		A.build(true);

		// 4. store results in cache
		// cache limit is ignored
		xp_cache_->clear();
		k = 0;
		for(
			wb_storage::iterator pwb = wb.begin(), wb_end = wb.end();
			pwb != wb_end, k < A.xbricks().size(); ++pwb, ++k
		) {
			(*xp_cache_)[pwb->first + pwb->second] = A.xbricks()[k].path();
		}

		return true;
	}

	void build(double date) {
		//s_.clear();

		// 1 fill storage with all unique well+branch
		if(!sw_) return;
		wb_storage wb;
		std::string q = (cd_traits::select_unique_well_branch() % date).str();
		sw_->prepare_sql(q);
		while(sw_->step_sql() == 0) {
			wb.insert(make_pair(sw_->get_sql_str(0), sw_->get_sql_str(1)));
		}
		sw_->finalize_sql();

		// 2 precalc plane size
		//const ulong plane_sz = m_size_[0] * m_size_[1];

		// 3 for each well+branch combo do
		for(wb_storage::iterator pwb = wb.begin(), wb_end = wb.end(); pwb != wb_end; ++pwb) {
			bool do_build_xp = true;
			xpath* xp;
			st_smart_ptr< xbuilder > sp_A;

			// check if this branch is already in cache
			if(xp_cache_) {
				xpc_iterator p_xpc = xp_cache_->find(pwb->first + pwb->second);
				if(p_xpc != xp_cache_->end()) {
					//cache hit
					xp = &p_xpc->second;
					do_build_xp = false;
				}
			}
			if(do_build_xp) {
				// 3.1 from 'branches' select well+branch_i trajectory (sql_well::get_branch_traj)
				sp_traj_t traj = sw_->get_branch_traj(pwb->first, pwb->second);
				if(!traj) return;
				sp_table_t traj_t = traj->get_table();
				if(!traj_t) return;

				// find column id's for X, Y, Z and MD
				const t_long col_ids[4] = {1, 2, 3, 0};
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
				sp_A = new xbuilder(m_, W);
				const hit_idx_t& hi = sp_A->build();
				//sp_A->remove_dups2();
				sp_A->append_wp_nodes(hi);
				xp = &sp_A->path();

				// if cache enabled - then store new xpath in cache
				if(xp_cache_ && (cache_limit_ == 0 || xp_cache_->size() < cache_limit_))
					(*xp_cache_)[pwb->first + pwb->second] = *xp;
			}

			if (xp->empty()) // no intersections
				continue;

			// 3.3 select all completions that belong to well+branch_i
			q = (cd_traits::select_segment() % date % pwb->first % pwb->second).str();
			sw_->prepare_sql(q);

			// 3.4 for all completions do
			while(sw_->step_sql() == 0) {
				cd_traits::proc_length(pwb->first, pwb->second, *xp, *this);
			} // 3.5 end of completions loop

			sw_->finalize_sql();
		} // 4 end of well+branch loop
	}

	trimesh m_;
	//vertex_pos_i m_size_;
	// tops should live as long as mesh lives
	//spv_float tops_;
	// copy of source sql_well
	sp_srcwell sw_;
	// storage for compdat records
	storage_t s_;
	// xpath cache
	spxp_cache_t xp_cache_;
	// max elems in cache
	ulong cache_limit_;
};

}} // eof blue_sky::fci


#endif /* end of include guard: FRAC_COMP_BUILDER_NIGNH6R8 */

