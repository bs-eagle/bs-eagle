/// @file frac_comp_ident.h
/// @brief Fractures and completions identification inside mesh
/// @author uentity
/// @version 
/// @date 06.10.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef FRAC_COMP_IDENT_3WM7KTWH
#define FRAC_COMP_IDENT_3WM7KTWH

#include "conf.h"
#include "well_pool_iface.h"
#include "wpi_algo.h"
#include "wpi_strategy_3d.h"

namespace blue_sky { 
namespace fci {

// impl details
using namespace std;
using namespace wpi;

typedef algo< strategy_3d > wpi_algo;
typedef wpi_algo::intersect_path xpath;
typedef std::vector <xpath>  xpath_storage;

/*-----------------------------------------------------------------
 * interface of well_path_builder building algo
 *----------------------------------------------------------------*/
class BS_API_PLUGIN well_path_builder {

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

	typedef well_pool_iface::sp_traj_t sp_traj_t;
	typedef well_pool_iface::sp_table_t sp_table_t;

  typedef xpath::iterator x_iterator;
	typedef multimap< string, string > wb_storage;

	// ctors
	well_path_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn);

	well_path_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn,
	                 	smart_ptr< well_pool_iface, true > src_well);

	void init(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn);
	void init(smart_ptr< well_pool_iface, true > src_well);

  template <class cd_traits> 
  xpath_storage build (const cd_traits& t);
	// clear storage
	void clear();
	
  trimesh m_;
	vertex_pos_i m_size_;
	// tops should live as long as mesh lives
	spv_float tops_;
	// copy of source sql_well
	smart_ptr< well_pool_iface, true > sw_;

private:
};

/*-----------------------------------------------------------------
 *  POD that holds info needed by COMPDAT
 *----------------------------------------------------------------*/
struct BS_API_PLUGIN compdat {
	typedef unsigned long ulong;
	typedef ulong cell_info[4];
	typedef ulong pos_i[3];
  typedef t_float coord[3];

	std::string well_name;
	std::string branch_name;
	cell_info cell_pos;
	char dir;
	t_double kh_mult;
  coord x1;       //!< start point coordinates (in 3D) of completion part
  coord x2;       //!< end point coordinates  (in 3D) of completion part

	// ctors
	compdat(const std::string& well_name_, const std::string& branch_name_,
		const pos_i& cell_pos_, const pos_i& mesh_size);
	// cell_pos decoded from cell_id
	compdat(const std::string& well_name_, const std::string& branch_name_,
		ulong cell_id, const pos_i& mesh_size);
	// for searching only
	compdat(ulong cell_id);

	// decode cell id
	void init(ulong cell_id, const pos_i& mesh_size);

	bool operator<(const compdat& rhs) const {
		return cell_id_ < rhs.cell_id_;
	}

private:
	ulong cell_id_;
};

typedef std::set< compdat > cd_storage;

/*-----------------------------------------------------------------
 * interface of COMPDAT building algo
 *----------------------------------------------------------------*/
class BS_API_PLUGIN compdat_builder {

public:
	// ctors
	compdat_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn);

	compdat_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn,
		smart_ptr< well_pool_iface, true > src_well);

	void init(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn);
	void init(smart_ptr< well_pool_iface, true > src_well);

	// mode == 0 - search completions, otherwise - fractures
	const cd_storage& build(double date, int mode = 0); 

	// clear storage
	void clear();

	// storage getter
	const cd_storage& storage() const;

private:
	class impl;
	st_smart_ptr< impl > pimpl_;
};

/*-----------------------------------------------------------------
 *  POD that holds info needed by FRACTURE
 *----------------------------------------------------------------*/
struct BS_API_PLUGIN fracture {
	typedef unsigned long ulong;
	typedef ulong cell_info[4];
	typedef ulong pos_i[3];
  typedef t_float coord[3];

	std::string well_name;
	std::string branch_name;
  std::string frac_status;
	cell_info cell_pos;
	t_double frac_half_length_1;
  t_double frac_half_length_2;
	t_double frac_angle;
  t_double frac_half_thin;
  t_double frac_perm;


	// ctors
	fracture(const std::string& well_name_, const std::string& branch_name_,
		const pos_i& cell_pos_, const pos_i& mesh_size);
	// cell_pos decoded from cell_id
	fracture(const std::string& well_name_, const std::string& branch_name_,
		ulong cell_id, const pos_i& mesh_size);
	// for searching only
	fracture(ulong cell_id);

	// decode cell id
	void init(ulong cell_id, const pos_i& mesh_size);

	bool operator<(const fracture& rhs) const {
		return cell_id_ < rhs.cell_id_;
	}

private:
	ulong cell_id_;
};

typedef std::set< fracture > frac_storage;

/*-----------------------------------------------------------------
 * interface of COMPDAT building algo
 *----------------------------------------------------------------*/
class BS_API_PLUGIN fracture_builder {

public:
	// ctors
	fracture_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn);

	fracture_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn,
		smart_ptr< well_pool_iface, true > src_well);

	void init(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn);
	void init(smart_ptr< well_pool_iface, true > src_well);

	// mode == 0 - search completions, otherwise - fractures
	const frac_storage& build(double date, int mode = 0);

	// clear storage
	void clear();

	// storage getter
	const frac_storage& storage() const;

private:
	class impl;
	st_smart_ptr< impl > pimpl_;
};


/*-----------------------------------------------------------------
 * global functions
 *----------------------------------------------------------------*/
spv_float completions_ident(smart_ptr< well_pool_iface, true > src_well, double date,
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn);

spv_float fractures_ident(smart_ptr< well_pool_iface, true > src_well, double date,
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn);

}}  // eof blue_sky::fci

#endif /* end of include guard: FRAC_COMP_IDENT_3WM7KTWH */

