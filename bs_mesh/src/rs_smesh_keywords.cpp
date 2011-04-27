/*!
  \file rs_smesh_keys.cpp
  \brief This file difines keyword handlers for bs reservoir simulation structured meshes
  \author Mark Khait
  \date 2009-08-07
*/

#include "bs_mesh_stdafx.h" 
#include "rs_smesh_keywords.h"

namespace blue_sky
  {
 
  #define KH_COMMON_VARIABLES_DEF           \
    std::ostringstream out_s;               \
    size_t len;                             \
    len = 0;

  //! this macro throw exception and assert with message from out_s
  #define KH_ASSERT_EXCEPTION \
    BS_ASSERT(false) (out_s.str());\
    throw bs_exception("Structured mesh handlers class",out_s.str().c_str());
  
  
  smesh_keywords::smesh_keywords(bs_type_ctor_param)
  {

  }

  
  smesh_keywords::smesh_keywords(const smesh_keywords& src)
  : bs_refcounter (src), base_t (src)
  {
    // TODO: BUG:
    bs_throw_exception ("NOT IMPL YET");
    //*this = src;
  }
  
  
  void smesh_keywords::register_keywords (sp_objbase &km, std::string provider) const
    {
      smart_ptr <keyword_manager_iface , true> keyword_manager (km);
      keyword_manager->register_supported_keyword ("DIMENS", provider);
      keyword_manager->register_supported_keyword ("ACTNUM", provider);
      keyword_manager->register_supported_keyword ("PERMX", provider);
      keyword_manager->register_supported_keyword ("PERMY", provider);
      keyword_manager->register_supported_keyword ("PERMZ", provider);
      keyword_manager->register_supported_keyword ("PORO", provider);
      keyword_manager->register_supported_keyword ("NTG", provider);
      keyword_manager->register_supported_keyword ("MULTX", provider);
      keyword_manager->register_supported_keyword ("MULTY", provider);
      keyword_manager->register_supported_keyword ("MULTZ", provider);
      keyword_manager->register_supported_keyword ("MULTPV", provider);
    }
  
  
  void smesh_keywords::activate_keywords (sp_km_iface_t keyword_manager)
    {
      //keyword_manager->register_keyword ("DIMENS", keyword_handler (&this_t::DIMENS_reactor));
      std::vector<std::string> names(6);
      npy_intp array_dimens[6] = {1,0,1,0,1,0};
      t_float def_value_zero = 0.0;
      t_float def_value_one = 1.0;
      
      names[0] = "nx";
      names[1] = "ny";
      names[2] = "nz";

      keyword_manager->register_prop_keyword ("DIMENS", "iii", names, &this_t::DIMENS_reactor);
      names[0] = "minimal_pore_volume";
      keyword_manager->register_prop_keyword ("MINPV", "f", names, 0);
      names[0] = "minimal_splice_volume";
      keyword_manager->register_prop_keyword ("MINSV", "f", names, 0);
      names[0] = "maximum_splice_thickness";
      keyword_manager->register_prop_keyword ("MAXST", "f", names, 0);

      keyword_manager->register_i_pool_keyword ("ACTNUM", &array_dimens[0], 1, 0);
      keyword_manager->register_fp_pool_keyword ("PERMX", &array_dimens[0], def_value_zero, 0);
      keyword_manager->register_fp_pool_keyword ("PERMY", &array_dimens[0], def_value_zero, 0);
      keyword_manager->register_fp_pool_keyword ("PERMZ", &array_dimens[0], def_value_zero, 0);
      keyword_manager->register_fp_pool_keyword ("PORO", &array_dimens[0], def_value_zero, 0);
      keyword_manager->register_fp_pool_keyword ("NTG", &array_dimens[0], def_value_one, 0);
      keyword_manager->register_fp_pool_keyword ("MULTX", &array_dimens[0], def_value_one, 0);
      keyword_manager->register_fp_pool_keyword ("MULTY", &array_dimens[0], def_value_one, 0);
      keyword_manager->register_fp_pool_keyword ("MULTZ", &array_dimens[0], def_value_one, 0);
      keyword_manager->register_fp_pool_keyword ("MULTPV", &array_dimens[0], def_value_zero, 0);
    }  
  
  void smesh_keywords::DIMENS_reactor(const std::string &keyword, keyword_params_t &params)
    {
      KH_COMMON_VARIABLES_DEF
      t_long ndim = 0, nblock = 0;
      t_long itmp[3];
      

      itmp[0] = params.hdm->get_prop ()->get_i ("nx");
      itmp[1] = params.hdm->get_prop ()->get_i ("ny");
      itmp[2] = params.hdm->get_prop ()->get_i ("nz");
      
      params.hdm->get_pool()->set_pool_dims (itmp, 3);
      
      // Number of nodes
      ndim = itmp[0] * itmp[1] * itmp[2];
      // Number of blocks
      nblock = itmp[0] * itmp[1] * itmp[2];
      
      BOSOUT (section::read_data, level::medium) <<
        "Keyword " << keyword <<
        ": NX = " << itmp[0] <<
        ", NY = " << itmp[1] <<
        ", NZ = " << itmp[2] << bs_end;
      BOSOUT (section::read_data,level::medium) << keyword << bs_end;
    }
  
  
  BLUE_SKY_TYPE_STD_CREATE (smesh_keywords)
  BLUE_SKY_TYPE_STD_COPY (smesh_keywords)
  BLUE_SKY_TYPE_IMPL (smesh_keywords, keyword_info_base,"BOS Core struct_mesh keyword_info_base", "struct_mesh", "Reservoir sumulator structured struct mesh keywords keywords")
      
}; //namespace blue_sky
