/*!
	\file bs_mesh_ijk.cpp
	\brief This file implement bs wrapper over mesh_ijk
  \author Mark Khait
  \date 2009-07-21
 */
 
#include "bs_mesh_stdafx.h"
#include "bs_mesh_ijk.h"
#include "bs_flux_connections.h"


namespace blue_sky
  {
  template<class strategy_t>
  bs_mesh_ijk<strategy_t>::bs_mesh_ijk(bs_type_ctor_param)
  {

  }

  template<class strategy_t>
  bs_mesh_ijk<strategy_t>::bs_mesh_ijk(const bs_mesh_ijk<strategy_t>& src)
  : bs_refcounter (src)//, objbase (src)
  {
    // TODO: BUG:
    bs_throw_exception ("NOT IMPL YET");
    //*this = src;
  }
  
  BLUE_SKY_TYPE_IMPL_T_EXT(1 , (bs_mesh_ijk<base_strategy_fif>) , 1, (rs_smesh_iface), "bs_mesh_ijk_fi", "Mesh base (virtual)  class", "Mesh base (virtual) class", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1 , (bs_mesh_ijk<base_strategy_did>) , 1, (rs_smesh_iface), "bs_mesh_ijk_di", "Mesh base (virtual)  class", "Mesh base (virtual)  class", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1 , (bs_mesh_ijk<base_strategy_dif>) , 1, (rs_smesh_iface), "bs_mesh_ijk_mixi", "Mesh base (virtual)  class", "Mesh base (virtual)  class", false);

  BLUE_SKY_TYPE_STD_CREATE_T_DEF(bs_mesh_ijk, (class));
  //BLUE_SKY_TYPE_STD_COPY_T_DEF(bs_mesh_ijk, (class));
  
  template <typename strategy_t>
  blue_sky::objbase* 
  bs_mesh_ijk <strategy_t>::bs_create_copy (bs_type_cpy_ctor_param src) 
  {
    const bs_mesh_ijk <strategy_t> *src_ptr = dynamic_cast <const bs_mesh_ijk <strategy_t> *> (src.get ());
    if (!src_ptr)
      bs_throw_exception ("Can't cast to bs_mesh_ijk");
      
    return new bs_mesh_ijk <strategy_t> (src_ptr);
  }
}; //namespace blue_sky
