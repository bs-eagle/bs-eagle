/*!
	\file mesh_base.cpp
	\brief This file implement base abstract class for bs meshes
  \author Mark Khait
  \date 2009-07-17
 */

#include "bs_mesh_stdafx.h"
#include "mesh_base.h"

template<class strategy_t>
mesh_base<strategy_t>::mesh_base()
{
  n_elements = 0;
  n_active_elements = 0;
  n_connections = 0;
}

template<class strategy_t>
void mesh_base<strategy_t>::check_data() const
{
  if (n_elements <= 0)
    bs_throw_exception (boost::format ("n_elements = %d is out of range")% n_elements);
  if (n_active_elements <= 0)
    bs_throw_exception (boost::format ("n_active_elements = %d is out of range")% n_active_elements);
  
//   if (n_connections <= 0)
//     bs_throw_exception (boost::format ("n_connections = %d is out of range")% n_connections);
    
  if (!volumes.size ())
    bs_throw_exception ("volumes array is not initialized");
  if (!ext_to_int.size ())
    bs_throw_exception ("ext_to_int array is not initialized");
}


BS_INST_STRAT(mesh_base);


