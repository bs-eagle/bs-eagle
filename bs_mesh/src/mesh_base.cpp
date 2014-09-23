/*!
	\file mesh_base.cpp
	\brief This file implement base abstract class for bs meshes
  \author Mark Khait
  \date 2009-07-17
 */
#include "mesh_base.h"


mesh_base::mesh_base()
{
#ifndef PURE_MESH
  volumes = give_kernel::Instance().create_object(v_float::bs_type());
  ext_to_int = give_kernel::Instance().create_object(v_long::bs_type());
  int_to_ext = give_kernel::Instance().create_object(v_long::bs_type());
#else
  volumes = 0;
  ext_to_int = 0;
  int_to_ext = 0;
#endif
  
  n_elements = 0;
  n_active_elements = 0;
  n_connections = 0;
}


void mesh_base::check_data() const
{

#ifndef PURE_MESH
  if (n_elements <= 0)
    bs_throw_exception (boost::format ("n_elements = %d is out of range")% n_elements);
  if (n_active_elements <= 0)
    bs_throw_exception (boost::format ("n_active_elements = %d is out of range")% n_active_elements);

  
//   if (n_connections <= 0)
//     bs_throw_exception (boost::format ("n_connections = %d is out of range")% n_connections);
    
  if (!volumes->size ())
    bs_throw_exception ("volumes array is not initialized");
  if (!ext_to_int->size ())
    bs_throw_exception ("ext_to_int array is not initialized");
 #endif
}


//BS_INST_STRAT(mesh_base);


