/*!
	\file bs_flux_connections.cpp
	\brief This file implement bs wrapper over flux_connections
  \author Mark Khait
  \date 2009-07-22
 */
 
#include "bs_flux_connections.h"

namespace blue_sky
  {
  
  bs_flux_connections::bs_flux_connections(bs_type_ctor_param)
  {

  }

  bs_flux_connections::bs_flux_connections(const bs_flux_connections& src)
  : bs_refcounter (src)
  {
    wrapped = src.wrapped;
  }

  BLUE_SKY_TYPE_IMPL(bs_flux_connections, objbase, "flux_connections", "Mesh flux connections class", "Mesh flux connections class");

  BLUE_SKY_TYPE_STD_CREATE(bs_flux_connections);
  BLUE_SKY_TYPE_STD_COPY(bs_flux_connections);
  
}; //namespace blue_sky
