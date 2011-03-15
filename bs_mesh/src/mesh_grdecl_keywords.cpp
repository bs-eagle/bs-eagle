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
 
   
  
  mesh_grdecl_keywords::mesh_grdecl_keywords(bs_type_ctor_param)
  {

  }

  
  mesh_grdecl_keywords::mesh_grdecl_keywords(const mesh_grdecl_keywords& src)
  : bs_refcounter (src), base_t (src)
  {
    // TODO: BUG:
    bs_throw_exception ("NOT IMPL YET");
    //*this = src;
  }
  
 
  
  void mesh_grdecl_keywords::register_keywords (sp_objbase &km, std::string provider) const
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
  
  
  void mesh_grdecl_keywords::activate_keywords(sp_objbase &km)
    {  
      sp_km_iface_t keyword_manager (km, bs_dynamic_cast ());
      t_int zcorn_dimens[6] = {2,0,2,0,2,0};
      t_int coord_dimens[6] = {1,1,1,1,0,6};
      t_float def_value = 0.0;
      
      keyword_manager->register_keyword ("ZCORN", keyword_handler (0, def_value, &zcorn_dimens[0]));  
      keyword_manager->register_keyword ("COORD", keyword_handler (0, def_value, &coord_dimens[0]));  
    }
    
  void mesh_grdecl_keywords::mesh_grdecl_handler(const std::string &keyword, keyword_params_t &params)
    {
      sp_idata_t idata (params.data, bs_dynamic_cast ());
      sp_bs_mesh_grdecl_t mesh_grdecl (BS_KERNEL.create_object (bs_mesh_grdecl ::bs_type ()), bs_dynamic_cast ());
      params.mesh = sp_objbase (mesh_grdecl);
      activate_keywords (params.km);
      base_t::activate_keywords (params.km);
    }
  
  BLUE_SKY_TYPE_STD_CREATE (mesh_grdecl_keywords)
  BLUE_SKY_TYPE_STD_COPY (mesh_grdecl_keywords)
  BLUE_SKY_TYPE_IMPL (mesh_grdecl_keywords, keyword_info_base, "BOS Core mesh_grdecl keyword_info", "MESH_GRDECL", "Reservoir sumulator structured mesh ijk keywords keywords")
  
    
}; //namespace blue_sky
