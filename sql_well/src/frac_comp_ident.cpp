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
#include "export_python_wrapper.h"
#include "py_list_converter.h"
#include "sql_well.h"


#include <boost/format.hpp>
#include <iosfwd>
#include <iostream>

namespace blue_sky { namespace fci {
// impl details
using namespace std;
using namespace wpi;

typedef algo< strategy_3d > wpi_algo;

// hidden details
namespace {

/*-----------------------------------------------------------------
 * search different tables
 *----------------------------------------------------------------*/
struct fract_traits  {
	static boost::format select_unique_well_branch() {
		return boost::format("SELECT DISTINCT well_name, branch_name FROM fractures WHERE d=%f");
	}

	static boost::format select_segment() {
		return boost::format(
			"SELECT status, md, half_length_1, half_length_2, angle, half_up, half_down, perm, half_thin FROM fractures WHERE d=%f and well_name='%s' and branch_name='%s'"
		);
	}
};

struct compl_traits {
	static boost::format select_unique_well_branch() {
		return boost::format("SELECT DISTINCT well_name, branch_name FROM completions WHERE d=%f");
	}

	static boost::format select_segment() {
		return boost::format(
			"SELECT md, length FROM completions WHERE d=%f and well_name='%s' and branch_name='%s'"
		);
	}
};

// simple dump
void dump(std::ostream& os, const compdat& cd) {
	os << setw(15) << cd.well_name << ' ' << setw(15) << cd.branch_name << ' ';
	for(int i = 0; i < 4; ++i)
		os << setw(4) << cd.cell_pos[i];
	os << setw(2) << cd.dir;
	os << setw(10) << std::fixed << std::setprecision(3) << cd.kh_mult << std::endl;
}

void dump(std::ostream& os, const cd_storage& cfs) {
	for(cd_storage::const_iterator pc = cfs.begin(), end = cfs.end(); pc != end; ++pc)
		dump(os, *pc);
}

}  // oef hidden namespace


/*-----------------------------------------------------------------
 * well_path_builder
 *----------------------------------------------------------------*/
well_path_builder::well_path_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn)
{
  init (nx, ny, coord, zcorn);
}

well_path_builder::well_path_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn,
	smart_ptr< well_pool_iface, true > src_well)
{
	init(nx, ny, coord, zcorn);
	init(src_well);
}

void 
well_path_builder::init(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn) 
{
  tops_ = wpi_algo::coord_zcorn2trimesh(nx, ny, coord, zcorn, m_, m_size_);
}

void 
well_path_builder::init(smart_ptr< well_pool_iface, true > src_well) 
{
	sw_ = BS_KERNEL.create_object_copy(src_well);
}

template <class cd_traits> 
xpath_storage well_path_builder::build (const cd_traits& t = cd_traits ()) 
{
	// 1 fill storage with all unique well+branch
  xpath_storage xpath_storage_;
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
    xpath_storage_.insert (xp);
  }
  return xpath_storage_;
}


/*-----------------------------------------------------------------
 * compdat
 *----------------------------------------------------------------*/
compdat::compdat(const string& well_name_, const string& branch_name_, const pos_i& cell_pos_, const pos_i& mesh_size)
	: well_name(well_name_), branch_name(branch_name_), dir(' '), kh_mult(0)
{
	copy(&cell_pos_[0], &cell_pos_[strategy_3d::D], &cell_pos[0]);
	cell_pos[3] = cell_pos[2];
	cell_id_ = wpi_algo::encode_cell_id(cell_pos_, mesh_size);
  
  x1[0] = x1[1] = x1[2] = 0;
  x2[0] = x2[1] = x2[2] = 0;
}

compdat::compdat(const string& well_name_, const string& branch_name_, ulong cell_id, const pos_i& mesh_size)
	: well_name(well_name_), branch_name(branch_name_), dir(' '), kh_mult(0)
{
	init(cell_id, mesh_size);
}

compdat::compdat(ulong cell_id)
	: dir(' '), kh_mult(0), cell_id_(cell_id)
{
  x1[0] = x1[1] = x1[2] = 0;
  x2[0] = x2[1] = x2[2] = 0;
}

void compdat::init(ulong cell_id, const pos_i& mesh_size) {
	// convert cell id -> cell pos
	pos_i cell_pos_;
	wpi_algo::decode_cell_id(cell_id, cell_pos_, mesh_size);
	copy(&cell_pos_[0], &cell_pos_[strategy_3d::D], &cell_pos[0]);
	cell_pos[3] = cell_pos[2];
	// store cell_id
	cell_id_ = cell_id;
  x1[0] = x1[1] = x1[2] = 0;
  x2[0] = x2[1] = x2[2] = 0;
}

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

			// 3.3 select all completions that belong to well+branch_i
			q = (cd_traits::select_segment() % date % pwb->first % pwb->second).str();
			sw_->prepare_sql(q);

			// 3.4 for all completions do
			while(sw_->step_sql() == 0) {
				// 3.4.1 search for completion_j begin_j and end_j using md_j and lentgh_j
        t_double cur_md = sw_->get_sql_real (0);
        t_double cur_len = sw_->get_sql_real (1);
				x_iterator px = xp.upper_bound(whc(cur_md));
				x_iterator xend = xp.upper_bound(whc(cur_md + cur_len));
				// always start with prev intersection
				if(px != xp.begin())
					--px;

				// 3.4.2 consider all intersections between begin_j and end_j
				for(; px != xend; ++px) {
					// prepare compdat
					compdat cf(pwb->first, pwb->second, px->cell, m_size_);

					// 3.4.3.1 calc delta between consequent xpoint_k and xpoint_(k + 1)
					// position to next point
					x_iterator pnext_x = px;
					++pnext_x;
					double delta_l = 0;
					if(pnext_x != xp.end())
						{
              // completion fully inside well segment
              if (cur_md >= px->md && (cur_md + cur_len) <= pnext_x->md)
                {
                  delta_l = cur_len;
                  for (t_uint j = 0; j < strat_t::D; ++j)
                    {
                      cf.x1[j] = px->where[j] + (cur_md - px->md) / (pnext_x->md - px->md) * (pnext_x->where[j] - px->where[j]);
                      cf.x2[j] = px->where[j] + (cur_md + cur_len - px->md) / (pnext_x->md - px->md) * (pnext_x->where[j] - px->where[j]);
                    }
                }
              // well segment fully inside completion
              else if (cur_md <= px->md && (cur_md + cur_len) >= pnext_x->md)
                {
                  delta_l = pnext_x->md - px->md;
                  for (t_uint j = 0; j < strat_t::D; ++j)
                    {
                      cf.x1[j] = px->where[j];
                      cf.x2[j] = pnext_x->where[j];
                    }
                } 
              // start of completion is inside well segment 
              // end of completion is out of well segment
              else if (cur_md >= px->md && (cur_md + cur_len) >= pnext_x->md)
                {
                  delta_l = pnext_x->md - cur_md;
                  for (t_uint j = 0; j < strat_t::D; ++j)
                    {
                      cf.x1[j] = px->where[j] + (cur_md - px->md) / (pnext_x->md - px->md) * (pnext_x->where[j] - px->where[j]);
                      cf.x2[j] = pnext_x->where[j];
                    }
                }
              // start of completion is outside well segment 
              // end of completion is inside of well segment
              else if (cur_md <= px->md && (cur_md + cur_len) <= pnext_x->md)
                {
                  delta_l = cur_md + cur_len - px->md;
                  for (t_uint j = 0; j < strat_t::D; ++j)
                    {
                      cf.x1[j] = px->where[j];
                      cf.x2[j] = px->where[j] + (cur_md + cur_len - px->md) / (pnext_x->md - px->md) * (pnext_x->where[j] - px->where[j]);
                    }
                }
            } 
          else
            {
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
							double cur_step = std::abs(cf.x2[i] - cf.x1[i]);
							if(max_step < cur_step) {
								max_step = cur_step;
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

/*-----------------------------------------------------------------
 * compdat_builder implementation
 *----------------------------------------------------------------*/
compdat_builder::compdat_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn)
	: pimpl_(new impl(nx, ny, coord, zcorn))
{}

compdat_builder::compdat_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn,
	smart_ptr< well_pool_iface, true > src_well)
	: pimpl_(new impl(nx, ny, coord, zcorn, src_well))
{}

void compdat_builder::init(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn) {
	pimpl_->init(nx, ny, coord, zcorn);
}

void compdat_builder::init(smart_ptr< well_pool_iface, true > src_well) {
	pimpl_->init(src_well);
}

const cd_storage& compdat_builder::build(double date, int mode) {
	if(mode == 0)
		pimpl_->build(date, compl_traits());
	else
		pimpl_->build(date, fract_traits());
	return storage();
}

const cd_storage& compdat_builder::storage() const {
	return pimpl_->cfs_;
}

void compdat_builder::clear() {
	pimpl_->cfs_.clear();
}


/*-----------------------------------------------------------------
 * fracture
 *----------------------------------------------------------------*/
fracture::fracture(const string& well_name_, const string& branch_name_, const pos_i& cell_pos_, const pos_i& mesh_size)
	: well_name(well_name_), branch_name(branch_name_), frac_half_length_1 (0), frac_half_length_2 (0),  frac_angle (0),  
    frac_half_thin (0), frac_status("SHUT"), frac_perm (0)
{
	copy(&cell_pos_[0], &cell_pos_[strategy_3d::D], &cell_pos[0]);
	cell_pos[3] = cell_pos[2];
	cell_id_ = wpi_algo::encode_cell_id(cell_pos_, mesh_size);
}

fracture::fracture(const string& well_name_, const string& branch_name_, ulong cell_id, const pos_i& mesh_size)
	: well_name(well_name_), branch_name(branch_name_), frac_half_length_1 (0), frac_half_length_2 (0),  frac_angle (0), 
    frac_half_thin (0), frac_status("SHUT"), frac_perm (0)
{
	init(cell_id, mesh_size);
}

fracture::fracture(ulong cell_id)
	: frac_half_length_1 (0), frac_half_length_2 (0),  frac_angle (0),  
    frac_half_thin (0), frac_status("SHUT"), frac_perm (0), cell_id_(cell_id)
{
}

void fracture::init(ulong cell_id, const pos_i& mesh_size) {
	// convert cell id -> cell pos
	pos_i cell_pos_;
	wpi_algo::decode_cell_id(cell_id, cell_pos_, mesh_size);
	copy(&cell_pos_[0], &cell_pos_[strategy_3d::D], &cell_pos[0]);
	cell_pos[3] = cell_pos[2];
	// store cell_id
	cell_id_ = cell_id;
}

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

			// 3.3 select all completions that belong to well+branch_i
			q = (cd_traits::select_segment() % date % pwb->first % pwb->second).str();
			sw_->prepare_sql(q);

			// 3.4 for all completions do
			while(sw_->step_sql() == 0) {
				// 3.4.1 search for completion_j begin_j and end_j using md_j and lentgh_j
        std::string status     = sw_->get_sql_str  (0);
        t_double md            = sw_->get_sql_real (1);
        t_double half_length_1 = sw_->get_sql_real (2);
        t_double half_length_2 = sw_->get_sql_real (3);
        t_double angle         = sw_->get_sql_real (4);
        t_double half_up       = sw_->get_sql_real (5);
        t_double half_down     = sw_->get_sql_real (6);
        t_double perm          = sw_->get_sql_real (7);
        t_double half_thin     = sw_->get_sql_real (8);

				x_iterator px = xp.upper_bound(whc(md));
				// position to next point
				x_iterator pnext_x = px;
				// always start with prev intersection
				if(px != xp.begin())
					--px;

  			// prepare compdat
	  		fracture frac(pwb->first, pwb->second, px->cell, m_size_);

        vertex_pos   frac_coords;
        vertex_pos_i cell_pos;
        wpi_algo::decode_cell_id (px->cell, cell_pos, m_size_);

        frac.cell_pos[0] = cell_pos[0];   // I - position
        frac.cell_pos[1] = cell_pos[1];   // J - position
        frac.frac_half_length_1 = half_length_1;
        frac.frac_half_length_2 = half_length_2;
        frac.frac_status = status;
        frac.frac_angle = angle;
        frac.frac_perm  = perm;
        frac.frac_half_thin = half_thin;

        for (t_uint j = 0; j < strat_t::D; ++j)
          {
            frac_coords[j] = px->where[j] + (md - px->md) / (pnext_x->md - px->md) * (pnext_x->where[j] - px->where[j]);
          }

        // find kw1 position of fracture
        for (t_ulong kw = cell_pos[2]; kw >= 0; kw--)
          {
            vertex_pos_i cell_up;
            copy (&cell_pos[0], &cell_pos[strat_t::D], &cell_up[0]);
            cell_up[2] = kw; 

            t_ulong k_cell_id_ = wpi_algo::encode_cell_id (cell_up, m_size_); 
            t_double z_top = 0;
            t_float* tops_data = &(*tops_)[0];

            // get mean of Z top plane 
            for (t_uint j = 0; j < 4; ++j)
              {
                z_top += 0.25 * tops_data[24 * k_cell_id_ + 3 * j + 2]; 
              }
            // fracture inside this layer
            if (frac_coords[2] - half_up > z_top) 
              {
                frac.cell_pos[2] = kw;  // Kw1 - position
                break;
              }
          }
        // find kw2 position of fracture
        for (t_ulong kw = cell_pos[2]; kw < m_size_[2]; kw++)
          {
            vertex_pos_i cell_down;
            copy (&cell_pos[0], &cell_pos[strat_t::D], &cell_down[0]);
            cell_down[2] = kw; 

            t_ulong k_cell_id_ = wpi_algo::encode_cell_id (cell_down, m_size_); 
            t_double z_down = 0;
            t_float* tops_data = &(*tops_)[0];

            // get mean of Z bottom plane 
            for (t_uint j = 4; j < 8; ++j)
              {
                z_down += 0.25 * tops_data[24 * k_cell_id_ + 3 * j + 2]; 
              }
            // fracture inside this layer
            if (frac_coords[2] + half_down < z_down) 
              {
                frac.cell_pos[2] = kw; // Kw2 position
                break;
              }
          }

						// if compdat for this cell is already added
						// then just update kh_mult
						// otherwise add new COMPDAT record
						// TODO: handle case of different directions inside one cell
						frac_storage::iterator pcd = cfs_.find(fracture(px->cell));
						if(pcd != cfs_.end()) {
							fracture& cur_cd = const_cast< fracture& >(*pcd);
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

/*-----------------------------------------------------------------
 * fracture_builder implementation
 *----------------------------------------------------------------*/
fracture_builder::fracture_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn)
	: pimpl_(new impl(nx, ny, coord, zcorn))
{}

fracture_builder::fracture_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn,
	smart_ptr< well_pool_iface, true > src_well)
	: pimpl_(new impl(nx, ny, coord, zcorn, src_well))
{}

void fracture_builder::init(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn) {
	pimpl_->init(nx, ny, coord, zcorn);
}

void fracture_builder::init(smart_ptr< well_pool_iface, true > src_well) {
	pimpl_->init(src_well);
}

const frac_storage& fracture_builder::build(double date, int mode) {
	if(mode == 0)
		pimpl_->build(date, compl_traits());
	else
		pimpl_->build(date, fract_traits());
	return storage();
}

const frac_storage& fracture_builder::storage() const {
	return pimpl_->cfs_;
}

void fracture_builder::clear() {
	pimpl_->cfs_.clear();
}

/*-----------------------------------------------------------------
 * global functions impl
 *----------------------------------------------------------------*/
spv_float completions_ident(smart_ptr< well_pool_iface, true > src_well, double date,
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn)
{
	compdat_builder b(nx, ny, coord, zcorn, src_well);
	b.build(date, 0);
	dump(std::cout, b.storage());
	//return BS_KERNEL.create_object(v_float::bs_type());
	return NULL;
}

spv_float fractures_ident(smart_ptr< well_pool_iface, true > src_well, double date,
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn)
{
	compdat_builder b(nx, ny, coord, zcorn, src_well);
	b.build(date, 1);
	dump(std::cout, b.storage());
	//return BS_KERNEL.create_object(v_float::bs_type());
	return NULL;
}

} // eof blue_sky::fci

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(cdb_build_overl, blue_sky::fci::compdat_builder::build, 1, 2)

namespace python {
namespace bp = boost::python;

static void print_compdat(const fci::compdat& cd) {
	fci::dump(std::cout, cd);
}

void py_export_compdat_ident() {
	using namespace fci;

	// export global functions
	bp::def("completions_ident", &fci::completions_ident);
	bp::def("fractures_ident", &fci::fractures_ident);

	// export compdat
	bp::class_< compdat >("compdat", bp::no_init)
		.def_readwrite("well_name", &compdat::well_name)
		.def_readwrite("branch_name", &compdat::branch_name)
		.def_readwrite("dir", &compdat::dir)
		.def("dump", &print_compdat)
	;

	// export cd_storage as list
	typedef bspy_converter< list_traits< cd_storage, 2 > > py_cds_conv;
	py_cds_conv::register_to_py();
	py_cds_conv::register_from_py();

	void (compdat_builder::*init1)(t_ulong, t_ulong, spv_float, spv_float) = &compdat_builder::init;
	void (compdat_builder::*init2)(smart_ptr< well_pool_iface, true >) = &compdat_builder::init;
	bp::class_< compdat_builder >("compdat_builder",
		bp::init< t_ulong, t_ulong, spv_float, spv_float >())
		.def(bp::init< t_ulong, t_ulong, spv_float, spv_float, smart_ptr< well_pool_iface, true> >())
		.def("init", init1)
		.def("init", init2)
		.def("build", &compdat_builder::build,
			bp::return_value_policy< bp::copy_const_reference >(),
			cdb_build_overl())
		.def("clear", &compdat_builder::clear)
		.def("storage", &compdat_builder::storage,
			bp::return_value_policy< bp::copy_const_reference >())
	;

}

} /* blue_sky::python */

} /* blue_sky */