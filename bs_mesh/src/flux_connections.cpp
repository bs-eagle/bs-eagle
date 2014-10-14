/*! \file flux_connections.cpp
	\brief This file implement interface class which transfers mesh data to the reservoir simulation process
	\author Iskhakov Ruslan
	\date 2008-05-20 */

#include "flux_connections.h"


//! default constructor

flux_connections::flux_connections()
{
  conn_trans = BS_KERNEL.create_object("bcsr_matrix");
  matrix_block_idx_plus = give_kernel::Instance().create_object(v_long::bs_type());
  matrix_block_idx_minus = give_kernel::Instance().create_object(v_long::bs_type());
}


//BS_INST_STRAT(flux_connections);
