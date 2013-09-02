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
  
  bs_mesh_ijk::bs_mesh_ijk(bs_type_ctor_param)
  {

  }

  
  bs_mesh_ijk::bs_mesh_ijk(const bs_mesh_ijk& src)
  : bs_refcounter (src)//, objbase (src)
  {
    // TODO: BUG:
    bs_throw_exception ("NOT IMPL YET");
    //*this = src;
  }
  BLUE_SKY_TYPE_STD_CREATE(bs_mesh_ijk);
  BLUE_SKY_TYPE_STD_COPY(bs_mesh_ijk);

  BLUE_SKY_TYPE_IMPL(bs_mesh_ijk, rs_smesh_iface, "bs_mesh_ijk_fi", "Mesh base (virtual)  class", "Mesh base (virtual) class");

  
  /*
  blue_sky::objbase* 
  bs_mesh_ijk ::bs_create_copy (bs_type_cpy_ctor_param src) 
  {
    const bs_mesh_ijk  *src_ptr = dynamic_cast <const bs_mesh_ijk  *> (src.get ());
    if (!src_ptr)
      bs_throw_exception ("Can't cast to bs_mesh_ijk");
      
    return new bs_mesh_ijk  (src_ptr);
  }
*/
}; //namespace blue_sky
