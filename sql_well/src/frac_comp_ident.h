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

spv_float completions_ident(smart_ptr< sql_well > sw, double date,
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn);

spv_float fractures_ident(smart_ptr< sql_well > sw, double date,
		t_ulong nx, t_ulong ny, spv_float coord, spv_float zcorn);

}}  // eof blue_sky::fci

#endif /* end of include guard: FRAC_COMP_IDENT_3WM7KTWH */

