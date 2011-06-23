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
 
   
  mesh_ijk_keywords::mesh_ijk_keywords(bs_type_ctor_param)
  {

  }

  mesh_ijk_keywords::mesh_ijk_keywords(const mesh_ijk_keywords& src)
  : bs_refcounter (src), base_t (src)
  {
    // TODO: BUG:
    bs_throw_exception ("NOT IMPL YET");
    //*this = src;
  }
  
  void mesh_ijk_keywords::mesh_ijk_handler(const std::string &keyword, keyword_params_t &params)
    {
      sp_idata_t idata (params.data, bs_dynamic_cast ());
      sp_bs_mesh_ijk_t ijk_mesh (BS_KERNEL.create_object (bs_mesh_ijk::bs_type ()), bs_dynamic_cast());
      params.mesh = sp_mesh_iface_t (ijk_mesh);
      activate_keywords (params.km);
      base_t::activate_keywords (params.km);
    }
  
  void mesh_ijk_keywords::register_keywords (sp_objbase &km, std::string provider) const
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
  
  void mesh_ijk_keywords::activate_keywords (sp_objbase &km)
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

    
  BLUE_SKY_TYPE_STD_CREATE (mesh_ijk_keywords)
  BLUE_SKY_TYPE_STD_COPY (mesh_ijk_keywords)
  BLUE_SKY_TYPE_IMPL (mesh_ijk_keywords, keyword_info_base, 
    "BOS Core mesh_ijk keyword_info", "MESH_IJK", "Reservoir sumulator structured mesh ijk keywords")
    
}; //namespace blue_sky
