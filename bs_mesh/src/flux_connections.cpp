ope/*! \file flux_connections.cpp
	\brief This file implement interface class which transfers mesh data to the reservoir simulation process
	\author Iskhakov Ruslan
	\date 2008-05-20 */
#include "bs_mesh_stdafx.h"

#include "flux_connections.h"
#include "strategy_name.h"


//! default constructor
template<typename strategy_t>
flux_connections<strategy_t>::flux_connections()
{
  std::string mtx_name = "bcsr_matrix_" + tools::strategy_name <strategy_t>::name ();
  conn_trans = BS_KERNEL.create_object(mtx_name);
  matrix_block_idx_plus = give_kernel::Instance().create_object(bs_array<i_type_t>::bs_type());
  matrix_block_idx_minus = give_kernel::Instance().create_object(bs_array<i_type_t>::bs_type());
}


BS_INST_STRAT(flux_connections);