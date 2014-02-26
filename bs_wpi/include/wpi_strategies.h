/// @file wpi_strategies.h
/// @brief Include all WPI strategies and make useful typedefs for each strategy of each trait
/// @author uentity
/// @version 1.0
/// @date 29.11.2013
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WPI_STRATEGIES_3MKRZELP
#define WPI_STRATEGIES_3MKRZELP

#include "wpi_strategy_traits.h"
#include "wpi_strategy_3d.h"
#include "wpi_strategy_2d.h"

namespace blue_sky { namespace wpi {

typedef strategy_3d_ex< carray_traits >              carray_3d;
typedef strategy_3d_ex< online_tops_traits >         onlinett_3d;
typedef strategy_3d_ex< online_tops_traits_bufpool > onlinett_bp_3d;
typedef strategy_3d_ex< sgrid_traits >               sgrid_3d;
typedef strategy_3d_ex< rgrid_traits >               rgrid_3d;

typedef strategy_2d_ex< carray_traits >              carray_2d;
typedef strategy_2d_ex< online_tops_traits >         onlinett_2d;
typedef strategy_2d_ex< online_tops_traits_bufpool > onlinett_bp_2d;
typedef strategy_2d_ex< sgrid_traits >               sgrid_2d;
typedef strategy_2d_ex< rgrid_traits >               rgrid_2d;

// shortcoming typedef
typedef strategy_2d_ex< online_tops_traits > strategy_2d;
typedef strategy_3d_ex< online_tops_traits > strategy_3d;

}} /* namespace blue_sky::wpi */

#endif /* end of include guard: WPI_STRATEGIES_3MKRZELP */

