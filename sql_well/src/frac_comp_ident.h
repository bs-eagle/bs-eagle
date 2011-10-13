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
#include "sql_well.h"

namespace blue_sky { namespace fci {

/*-----------------------------------------------------------------
 *  POD that holds info needed by COMPDAT
 *----------------------------------------------------------------*/
struct BS_API_PLUGIN compdat {
	typedef unsigned long ulong;
	typedef ulong cell_info[4];
	typedef ulong pos_i[3];

	std::string well_name;
	std::string branch_name;
	cell_info cell_pos;
	char dir;
	t_double kh_mult;

	// ctors
	compdat(const std::string& well_name_, const std::string& branch_name_,
		const pos_i& cell_pos_, const pos_i& mesh_size);
	// cell_pos decoded from cell_id
	compdat(const std::string& well_name_, const std::string& branch_name_,
		ulong cell_id, const pos_i& mesh_size);
	// for searching only
	compdat(ulong cell_id);

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
		smart_ptr< sql_well, true > src_well);

	void init(t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn);
	void init(smart_ptr< sql_well, true > src_well);

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
 * global functions
 *----------------------------------------------------------------*/
spv_float completions_ident(smart_ptr< sql_well > src_well, double date,
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn);

spv_float fractures_ident(smart_ptr< sql_well > src_well, double date,
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn);

}}  // eof blue_sky::fci

#endif /* end of include guard: FRAC_COMP_IDENT_3WM7KTWH */

