/*!
  \file rs_smesh_keys.cpp
  \brief This file defines keyword handlers for bs reservoir simulation mesh IJK
  \author Mark Khait
  \date 2009-08-07
*/

#include "bs_mesh_stdafx.h" 
#include "mesh_ijk_keywords.h"

namespace blue_sky
  {
 
   
  template<class strategy_t>
  mesh_ijk_keywords<strategy_t>::mesh_ijk_keywords(bs_type_ctor_param)
  {

  }

  template<class strategy_t>
  mesh_ijk_keywords<strategy_t>::mesh_ijk_keywords(const mesh_ijk_keywords<strategy_t>& src)
  : bs_refcounter (src), base_t (src)
  {
    // TODO: BUG:
    bs_throw_exception ("NOT IMPL YET");
    //*this = src;
  }
  
  template <class strategy_t>
  void mesh_ijk_keywords<strategy_t>::mesh_ijk_handler(const std::string &keyword, keyword_params_t &params)
    {
      sp_idata_t idata (params.data, bs_dynamic_cast ());
      sp_bs_mesh_ijk_t ijk_mesh (BS_KERNEL.create_object (bs_mesh_ijk <strategy_t>::bs_type ()), bs_dynamic_cast());
      params.mesh = sp_mesh_iface_t (ijk_mesh);
      activate_keywords (params.km);
      base_t::activate_keywords (params.km);
    }
  
  template<class strategy_t>
  void mesh_ijk_keywords<strategy_t>::register_keywords (sp_objbase &km, std::string provider) const
    {
      sp_km_iface_t keyword_manager (km, bs_dynamic_cast ());
      if (provider == "")
        {
          provider = "MESH_IJK";
          keyword_manager->register_keyword ("MESH_IJK", keyword_handler (&this_t::mesh_ijk_handler));
        }
      keyword_manager->register_supported_keyword ("DX", provider);
      keyword_manager->register_supported_keyword ("DY", provider);
      keyword_manager->register_supported_keyword ("DZ", provider);
      keyword_manager->register_supported_keyword ("TOPS", provider);
      base_t::register_keywords(km, provider);
    }
  
  template<class strategy_t>
  void mesh_ijk_keywords<strategy_t>::activate_keywords (sp_objbase &km)
    {
      sp_km_iface_t keyword_manager (km, bs_dynamic_cast ());
      index_t dx_dimens[] = {1,0,1,0,1,0};
      index_t tops_dimens[] = {1,0,1,0,0,1};
      item_t def_value = 0.0;
      
      keyword_manager->register_keyword ("DX", keyword_handler (0, def_value, &dx_dimens[0]));
      keyword_manager->register_keyword ("DY", keyword_handler (0, def_value, &dx_dimens[0]));
      keyword_manager->register_keyword ("DZ", keyword_handler (0, def_value, &dx_dimens[0]));
      keyword_manager->register_keyword ("TOPS", keyword_handler (0, def_value, &tops_dimens[0]));
    }  

    
  BLUE_SKY_TYPE_STD_CREATE_T_DEF (mesh_ijk_keywords, (class))
  BLUE_SKY_TYPE_STD_COPY_T_DEF (mesh_ijk_keywords, (class))
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (mesh_ijk_keywords<base_strategy_fi>), 1, (keyword_info_base<base_strategy_fi>), 
    "BOS Core mesh_ijk keyword_info_fi", "MESH_IJK", "Reservoir sumulator structured mesh ijk keywords keywords", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (mesh_ijk_keywords<base_strategy_di>), 1, (keyword_info_base<base_strategy_di>), 
    "BOS_Core mesh_ijk keyword_info_di", "MESH_IJK", "Reservoir sumulator structured mesh ijk keywords keywords", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (mesh_ijk_keywords<base_strategy_mixi>), 1, (keyword_info_base<base_strategy_mixi>), 
    "BOS_Core mesh_ijk keyword_info_mixi", "MESH_IJK", "Reservoir sumulator structured mesh ijk keywords keywords", false)
    
}; //namespace blue_sky
