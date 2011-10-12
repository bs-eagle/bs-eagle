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

#include <boost/format.hpp>
#include <iosfwd>
#include <iostream>

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

/*-----------------------------------------------------------------
 * implementation of COMPDAT building algo
 *----------------------------------------------------------------*/
template< class cd_traits >
class compdat_builder {
public:
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

	compdat_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn) {
		init(nx, ny, coord, zcorn);
	}

	compdat_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn,
		smart_ptr< sql_well > src_well)
	{
		init(nx, ny, coord, zcorn);
		init(src_well);
	}

	void init(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn) {
		// build trimesh for given coord+zcorn
		tops_ = wpi_algo::coord_zcorn2trimesh(nx, ny, coord, zcorn, m_, m_size_, true);
	}

	void init(smart_ptr< sql_well > src_well) {
		sw_ = BS_KERNEL.create_object_copy(src_well);
	}

	void go(double date) {
		cfs_.clear();

		// 1 fill storage with all unique well+branch
		wb_storage wb;
		std::string q = (cd_traits::select_unique_well_branch() % date).str();
		sw_->prepare_sql(q);
		while(sw_->step_sql() == 0) {
			wb.insert(make_pair(sw_->get_sql_str(0), sw_->get_sql_str(1)));
		}
		sw_->finalize_sql();

		// 2 precalc plane size
		const ulong plane_sz = m_size_[0] * m_size_[1];

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
			xpath& xp = A.path();

			// 3.3 select all completions that belong to well+branch_i
			q = (cd_traits::select_segment() % date % pwb->first % pwb->second).str();
			sw_->prepare_sql(q);

			// 3.4 for all completions do
			while(sw_->step_sql() == 0) {
				// 3.4.1 search for completion_j begin_j and end_j using md_j and lentgh_j
				x_iterator px = xp.upper_bound(whc(sw_->get_sql_real(0)));
				x_iterator xend = xp.upper_bound(whc(sw_->get_sql_real(0) + sw_->get_sql_real(1)));
				// always start with prev intersection
				if(px != xp.begin())
					--px;

				// 3.4.2 consider all intersections between begin_j and end_j
				for(; px != xend; ++px) {
					// prepare compfrac
					vertex_pos_i cell_id;
					wpi_algo::decode_cell_id(px->cell, cell_id, m_size_);
					compfrac cf(pwb->first, pwb->second, cell_id);

					// 3.4.3.1 calc delta between consequent xpoint_k and xpoint_(k + 1)
					// position to next point
					x_iterator pnext_x = px;
					++pnext_x;
					ulong delta = 0;
					if(pnext_x != xp.end())
						delta = std::abs(pnext_x->cell - px->cell);

					// 3.4.3.2 if delta == 1 mark direction as 'X'
					//         else if delta == dx direction = 'Y'
					//         else if delta = dx*dy direction = 'Z'
					if(delta) {
						if(delta == 1)
							cf.dir = 'X';
						else if(delta >= m_size_[0] && delta < plane_sz)
							cf.dir = 'Y';
						else
							cf.dir = 'Z';
						// add new COMPDAT record
						cfs_.push_back(cf);
					}
				} // 3.4.4 end of intersections loop
			} // 3.5 end of completions loop

			sw_->finalize_sql();
		} // 4 end of well+branch loop
	}

	void dump(std::ostream& os) {
		for(cf_storage::iterator pc = cfs_.begin(), end = cfs_.end(); pc != end; ++pc) {
			os << setw(15) << pc->well_name << ' ' << setw(15) << pc->branch_name << ' ';
			for(int i = 0; i < 4; ++i)
				os << setw(3) << pc->cell_id[i] << ' ';
			os << pc->dir << std::endl;
		}
	}

	//virtual std::string select_unique_well_branch(double d) = 0;
	//virtual std::string select_segment(double d, const std::string well_name,
	//	const std::string& branch_name) = 0;

private:
	trimesh m_;
	vertex_pos_i m_size_;
	// tops should live as long as mesh lives
	spv_float tops_;
	// copy of source sql_well
	smart_ptr< sql_well > sw_;
	// storage for compdat records
	cf_storage cfs_;
};

/*-----------------------------------------------------------------
 * traits for completions
 *----------------------------------------------------------------*/
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

struct fract_traits {
	static boost::format select_unique_well_branch() {
		return boost::format("SELECT DISTINCT well_name, branch_name FROM fractures WHERE d=%f");
	}

	static boost::format select_segment() {
		return boost::format(
			"SELECT md, half_length1 + half_length_2 FROM fractures WHERE d=%f and well_name='%s' and branch_name='%s'"
		);
	}
};

} /* hidden namespace*/

spv_float completions_ident(smart_ptr< sql_well > src_well, double date,
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn)
{
	typedef compdat_builder< compl_traits > builder_t;

	builder_t b(nx, ny, coord, zcorn, src_well);
	b.go(date);
	b.dump(std::cout);
	//return BS_KERNEL.create_object(v_float::bs_type());
	return NULL;
}

spv_float fractures_ident(smart_ptr< sql_well > src_well, double date,
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn)
{
	typedef compdat_builder< fract_traits > builder_t;

	builder_t b(nx, ny, coord, zcorn, src_well);
	b.go(date);
	b.dump(std::cout);
	//return BS_KERNEL.create_object(v_float::bs_type());
	return NULL;
}

} // eof blue_sky::fci

namespace python {
namespace bp = boost::python;

void py_export_compdat_ident() {
	bp::def("completions_ident", &fci::completions_ident);
	bp::def("fractions_ident", &fci::fractures_ident);
}

} /* blue_sky::python */

} /* blue_sky */
