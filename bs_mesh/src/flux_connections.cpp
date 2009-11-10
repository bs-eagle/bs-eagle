/*! \file flux_connections.cpp
	\brief This file implement interface class which transfers mesh data to the reservoir simulation process
	\author Iskhakov Ruslan
	\date 2008-05-20 */
#include "bs_mesh_stdafx.h"

#include "flux_connections.h"


//! default constructor
template<typename strategy_t>
flux_connections<strategy_t>::flux_connections()
{
  conn_trans = BS_KERNEL.create_object (csr_matrix_t::bs_type ());
}


BS_INST_STRAT(flux_connections);