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
  template<class strategy_t>
  bs_flux_connections<strategy_t>::bs_flux_connections(bs_type_ctor_param)
  {

  }

  template<class strategy_t>
  bs_flux_connections<strategy_t>::bs_flux_connections(const bs_flux_connections<strategy_t>& src)
  : bs_refcounter (src), objbase (src)
  {
    // TODO: BUG:
    bs_throw_exception ("NOT IMPL YET");
    //*this = src;
  }
  
  BLUE_SKY_TYPE_IMPL_T_EXT(1 , (bs_flux_connections<base_strategy_fi>) , 1, (objbase), "bs_flux_connections_fi", "Mesh base (virtual)  class", "Mesh base (virtual) class", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1 , (bs_flux_connections<base_strategy_di>) , 1, (objbase), "bs_flux_connections_di", "Mesh base (virtual)  class", "Mesh base (virtual)  class", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1 , (bs_flux_connections<base_strategy_mixi>) , 1, (objbase), "bs_flux_connections_mixi", "Mesh base (virtual)  class", "Mesh base (virtual)  class", false);

  BLUE_SKY_TYPE_STD_CREATE_T_DEF(bs_flux_connections, (class));
  //BLUE_SKY_TYPE_STD_COPY_T_DEF(bs_flux_connections, (class));
  
  template <typename strategy_t>
  blue_sky::objbase* 
  bs_flux_connections <strategy_t>::bs_create_copy (bs_type_cpy_ctor_param src) 
  {
    const bs_flux_connections <strategy_t> *src_ptr = dynamic_cast <const bs_flux_connections <strategy_t> *> (src.get ());
    if (!src_ptr)
      bs_throw_exception ("Can't cast to bs_flux_connections");
      
    return new bs_flux_connections <strategy_t> (src_ptr);
  }

}; //namespace blue_sky
