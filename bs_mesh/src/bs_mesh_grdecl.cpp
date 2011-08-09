/*!
	\file bs_mesh_grdecl.cpp
	\brief This file implement bs wrapper over mesh_grdecl
  \author Mark Khait
  \date 2009-07-21
 */
 
#include "bs_mesh_stdafx.h"
#include "bs_mesh_grdecl.h"

namespace blue_sky
  {
  
  bs_mesh_grdecl::bs_mesh_grdecl(bs_type_ctor_param)
  {

  }

  
  bs_mesh_grdecl::bs_mesh_grdecl(const bs_mesh_grdecl& src)
  : bs_refcounter (src)
  {
    // TODO: BUG:
    bs_throw_exception ("NOT IMPL YET");
    //*this = src;
  }
  

  BLUE_SKY_TYPE_IMPL(bs_mesh_grdecl, rs_smesh_iface, "bs_mesh_grdecl", "Mesh base (virtual)  class", "Mesh base (virtual) class");
//  BLUE_SKY_TYPE_IMPL_T_EXT(1 , (bs_mesh_grdecl<base_strategy_did>) , 1, (rs_smesh_iface), "bs_mesh_grdecl_did", "Mesh base (virtual)  class", "Mesh base (virtual)  class", false);
//  BLUE_SKY_TYPE_IMPL_T_EXT(1 , (bs_mesh_grdecl<base_strategy_dif>) , 1, (rs_smesh_iface), "bs_mesh_grdecl_dif", "Mesh base (virtual)  class", "Mesh base (virtual)  class", false);
  
//  BLUE_SKY_TYPE_IMPL_T_EXT(1 , (bs_mesh_grdecl<base_strategy_flf>) , 1, (rs_smesh_iface), "base_strategy_flf", "Mesh base (virtual)  class", "Mesh base (virtual) class", false);
//  BLUE_SKY_TYPE_IMPL_T_EXT(1 , (bs_mesh_grdecl<base_strategy_dld>) , 1, (rs_smesh_iface), "base_strategy_dld", "Mesh base (virtual)  class", "Mesh base (virtual)  class", false);
//  BLUE_SKY_TYPE_IMPL_T_EXT(1 , (bs_mesh_grdecl<base_strategy_dlf>) , 1, (rs_smesh_iface), "base_strategy_dlf", "Mesh base (virtual)  class", "Mesh base (virtual)  class", false);

  BLUE_SKY_TYPE_STD_CREATE(bs_mesh_grdecl);
  BLUE_SKY_TYPE_STD_COPY(bs_mesh_grdecl);
  
/*
  
  blue_sky::objbase* 
  bs_mesh_grdecl ::bs_create_copy (bs_type_cpy_ctor_param src) 
  {
    const bs_mesh_grdecl  *src_ptr = dynamic_cast <const bs_mesh_grdecl  *> (src.get ());
    if (!src_ptr)
      bs_throw_exception ("Can't cast to bs_mesh_grdecl");
      
    return new bs_mesh_grdecl  (src_ptr);
  }
*/
}; //namespace blue_sky
