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
#include "wpi_strategies.h"

namespace blue_sky { namespace fci {

// impl details
using namespace std;
using namespace wpi;

typedef strategy_3d strat_t;
typedef algo< strat_t > wpi_algo;
typedef wpi_algo::intersect_path xpath;
typedef std::vector <xpath>  xpath_storage;

/*-----------------------------------------------------------------
 *  POD that holds info needed by COMPDAT
 *----------------------------------------------------------------*/
struct BS_API_PLUGIN compdat {
	typedef std::set< compdat > storage_t;
	typedef wpi::ulong ulong;
	typedef ulong cell_info[4];
	typedef ulong pos_i[3];
	typedef t_float coord[3];

	std::string well_name;
	std::string branch_name;
	cell_info   cell_pos;
	char        dir;
	t_double    kh_mult;
	t_double    md;
	t_double    len;
	t_double    skin;
	t_double    diam;
	t_double    kh;
	t_uint      status;
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
		if(well_name == rhs.well_name) {
			if(md == rhs.md)
				return cell_id_ < rhs.cell_id_;
			return md < rhs.md;
		}
		else
			return well_name < rhs.well_name;
	}

	static storage_t::const_iterator find_first_cd(const storage_t& cds, const std::string& well_name) {
		compdat t(0);
		t.well_name = well_name;
		const storage_t::const_iterator p_cd = cds.lower_bound(t);
		if(p_cd != cds.end() && p_cd->well_name == well_name)
			return p_cd;
		else
			return cds.end();
	}

private:
	ulong cell_id_;
};

typedef compdat::storage_t cd_storage;

/*-----------------------------------------------------------------
 *  POD that holds info needed by FRACTURE
 *----------------------------------------------------------------*/
struct BS_API_PLUGIN fracture {
	typedef std::set< fracture > storage_t;
	typedef wpi::ulong ulong;
	typedef ulong cell_info[4];
	typedef ulong pos_i[3];
	typedef t_float coord[3];

	std::string well_name;
	std::string branch_name;
	t_int       frac_status;
	cell_info   cell_pos;
	pos_i       md_cell_pos;
	t_double    frac_half_length_1;
	t_double    frac_half_length_2;
	t_double    frac_angle;
	t_double    frac_half_thin;
	t_double    frac_perm;
	t_double    frac_skin;
	t_int       frac_main_k;

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

typedef fracture::storage_t frac_storage;

/*-----------------------------------------------------------------
 * interface of COMPDAT/FRACTURE building algo
 *----------------------------------------------------------------*/
// brick - compdat or fracture
template< class brick >
class BS_API_PLUGIN builder {
public:
	template< class B > friend class builder;

	typedef typename brick::storage_t storage_t;

	builder();

	void init(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn);
	void init(t_ulong nx, t_ulong ny, sp_obj trim_backend);
	void init(smart_ptr< well_pool_iface, true > src_well);

	// invoke main algo
	const storage_t& build(double date);

	// clear storage
	void clear();

	// storage getter
	const storage_t& storage() const;

	template< class B >
	void share_cache_with(const builder< B >& rhs);

private:
	class impl;
	st_smart_ptr< impl > pimpl_;
};


/*-----------------------------------------------------------------
 * interface of COMPDAT building algo
 *----------------------------------------------------------------*/
class BS_API_PLUGIN compdat_builder : public builder< compdat > {
public:
	// ctors
	compdat_builder();

	compdat_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn);

	compdat_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn,
		smart_ptr< well_pool_iface, true > src_well);
};

/*-----------------------------------------------------------------
 * interface of FRACTURE building algo
 *----------------------------------------------------------------*/
class BS_API_PLUGIN fracture_builder : public builder< fracture > {
public:
	// ctors
	fracture_builder();

	fracture_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn);

	fracture_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn,
		smart_ptr< well_pool_iface, true > src_well);
};


/*-----------------------------------------------------------------
 * interface of COMPDAT and FRACTURE building algo
 *----------------------------------------------------------------*/
class BS_API_PLUGIN compl_n_frac_builder {

public:
	// ctors
	compl_n_frac_builder();

	compl_n_frac_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn);

	compl_n_frac_builder(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn,
		smart_ptr< well_pool_iface, true > src_well);

	void init(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn);
	void init(t_ulong nx, t_ulong ny, sp_obj trim_backend);
	void init(smart_ptr< well_pool_iface, true > src_well);

	// compdat build
	const cd_storage& compl_build(double date);

	// fractures build
	const frac_storage& frac_build(double date);

	// precalc all possible mesh-well intersections at once
	bool build_cache();

	// clear storage
	void clear();

	// storage getter
	const cd_storage& storage_compdat() const;

	// storage getter
	const frac_storage& storage_fracture() const;

	// use single cache among multiple instances
	void share_cache_with(const compl_n_frac_builder& rhs);

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

