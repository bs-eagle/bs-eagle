/*!
	\file bs_flux_connections.cpp
	\brief This file implement bs wrapper over flux_connections
  \author Mark Khait
  \date 2009-07-22
 */
 
#include "bs_mesh_stdafx.h"
#include "bs_flux_connections.h"

namespace blue_sky
  {
  bs_flux_connections::bs_flux_connections(bs_type_ctor_param)
  {

  }

  bs_flux_connections::bs_flux_connections(const bs_flux_connections& src)
  : bs_refcounter (src), objbase (src)
  {
    // TODO: BUG:
    bs_throw_exception ("NOT IMPL YET");
    //*this = src;
  }
  
  BLUE_SKY_TYPE_IMPL (bs_flux_connections, objbase, "bs_flux_connections", "bs_flux_connections", "bs_flux_connections");

  BLUE_SKY_TYPE_STD_CREATE (bs_flux_connections);
  //BLUE_SKY_TYPE_STD_COPY (bs_flux_connections);
  
  blue_sky::objbase* 
  bs_flux_connections::bs_create_copy (bs_type_cpy_ctor_param src) 
  {
    const bs_flux_connections *src_ptr = dynamic_cast <const bs_flux_connections *> (src.get ());
    if (!src_ptr)
      bs_throw_exception ("Can't cast to bs_flux_connections");
      
    return new bs_flux_connections (src_ptr);
  }

}; //namespace blue_sky
