/// @file frac_comp_ient.cpp
/// @brief Implementation of finding fractions/completions intersections with mesh
/// @author uentity
/// @version 
/// @date 06.10.2011
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "frac_comp_ident.h"

namespace blue_sky { namespace fci {

spv_float completions_ident(smart_ptr< sql_well > sw, double date,
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn)
{
	// 1 find all unique well+branch names on given date
	// 2 build trimesh for given coord+zcorn
	// 3 for each well+branch combo do
	// 3.1 from 'branches' select well+branch_i trajectory (sql_well::get_branch_traj)
	// 3.2 find intersections of given branch with mesh (well_path_ident)
	// 3.3 select all completions that belong to well+branch_i
	// 3.4 for all completions do
	// 3.4.1 search for completion_j begin_j and end using md_j and lentgh_j
	// 3.4.2 consider all intersections between begin_j and end_j
	// 3.4.3.1 calc delta between consequent xpoint_k and xpoint_(k + 1)
	// 3.4.3.2 if delta == 1 mark direction as 'X'
	//         else if delta == dx direction = 'Y'
	//         else if delta = dx*dy direction = 'Z'
	// 3.4.3.3 put new record on each delta change
	// 3.4.4 end of intersections loop
	// 3.5 end of completions loop
	// 4 end of well+branch loop
}

}} /* blue_sky::fci */
