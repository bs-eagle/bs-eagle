/*!
  \file rs_smesh_keys.cpp
  \brief This file defines keyword handlers for bs reservoir simulation mesh GRDECL
  \author Mark Khait
  \date 2009-08-07
*/

#include "bs_mesh_stdafx.h" 
#include "mesh_grdecl_keywords.h"
#include "bs_mesh_grdecl.h"

namespace blue_sky
  {
 
   
  template<class strategy_t>
  mesh_grdecl_keywords<strategy_t>::mesh_grdecl_keywords(bs_type_ctor_param)
  {

  }

  template<class strategy_t>
  mesh_grdecl_keywords<strategy_t>::mesh_grdecl_keywords(const mesh_grdecl_keywords<strategy_t>& src)
  : bs_refcounter (src), base_t (src)
  {
    // TODO: BUG:
    bs_throw_exception ("NOT IMPL YET");
    //*this = src;
  }
  
 
  template<class strategy_t>
  void mesh_grdecl_keywords<strategy_t>::register_keywords (sp_objbase &km, std::string provider) const
    {
      sp_km_iface_t keyword_manager (km, bs_dynamic_cast ());
      if (provider == "")
        {
          provider = "MESH_GRDECL";
          keyword_manager->register_keyword ("MESH_GRDECL", keyword_handler (&this_t::mesh_grdecl_handler));
        }
      keyword_manager->register_supported_keyword ("ZCORN", provider);
      keyword_manager->register_supported_keyword ("COORD", provider);
      base_t::register_keywords(km, provider);
    }
  
  template<class strategy_t>
  void mesh_grdecl_keywords<strategy_t>::activate_keywords(sp_objbase &km)
    {  
      sp_km_iface_t keyword_manager (km, bs_dynamic_cast ());
      index_t zcorn_dimens[6] = {2,0,2,0,2,0};
      index_t coord_dimens[6] = {1,1,1,1,0,6};
      item_t def_value = 0.0;
      keyword_manager->register_keyword ("ZCORN", keyword_handler (0, def_value, &zcorn_dimens[0]));  
      keyword_manager->register_keyword ("COORD", keyword_handler (0, def_value, &coord_dimens[0]));  
    }
    
  template <class strategy_t>
  void mesh_grdecl_keywords<strategy_t>::mesh_grdecl_handler(const std::string &keyword, keyword_params_t &params)
    {
      sp_idata_t idata (params.data, bs_dynamic_cast ());
      sp_bs_mesh_grdecl_t mesh_grdecl (BS_KERNEL.create_object (bs_mesh_grdecl <strategy_t>::bs_type ()), bs_dynamic_cast ());
      params.mesh = sp_objbase (mesh_grdecl);
      activate_keywords (params.km);
      base_t::activate_keywords (params.km);
    }
  
  BLUE_SKY_TYPE_STD_CREATE_T_DEF (mesh_grdecl_keywords, (class))
  BLUE_SKY_TYPE_STD_COPY_T_DEF (mesh_grdecl_keywords, (class))
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (mesh_grdecl_keywords<base_strategy_fi>), 1, (keyword_info_base<base_strategy_fi>), 
    "BOS Core mesh_grdecl keyword_info_fi", "MESH_GRDECL", "Reservoir sumulator structured mesh ijk keywords keywords", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (mesh_grdecl_keywords<base_strategy_di>), 1, (keyword_info_base<base_strategy_di>), 
    "BOS_Core mesh_grdecl keyword_info_di", "MESH_GRDECL", "Reservoir sumulator structured mesh ijk keywords keywords", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (mesh_grdecl_keywords<base_strategy_mixi>), 1, (keyword_info_base<base_strategy_mixi>), 
    "BOS_Core mesh_grdecl keyword_info_mixi", "MESH_GRDECL", "Reservoir sumulator structured mesh ijk keywords keywords", false)
    
}; //namespace blue_sky
