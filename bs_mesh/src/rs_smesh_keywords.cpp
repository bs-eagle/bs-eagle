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
  
  
  void smesh_keywords::activate_keywords(sp_objbase &km)
    {
      smart_ptr <keyword_manager_iface , true> keyword_manager (km);
      keyword_manager->register_keyword ("DIMENS", keyword_handler (&this_t::DIMENS_handler));
      t_int array_dimens[6] = {1,0,1,0,1,0};
      t_float def_value_zero = 0.0;
      t_float def_value_one = 1.0;
      
      keyword_manager->register_keyword ("ACTNUM", keyword_handler (0, 1, &array_dimens[0]));
      keyword_manager->register_keyword ("PERMX", keyword_handler (0, def_value_zero, &array_dimens[0]));
      keyword_manager->register_keyword ("PERMY", keyword_handler (0, def_value_zero, &array_dimens[0]));
      keyword_manager->register_keyword ("PERMZ", keyword_handler (0, def_value_zero, &array_dimens[0]));
      keyword_manager->register_keyword ("PORO", keyword_handler (0, def_value_zero, &array_dimens[0]));
      keyword_manager->register_keyword ("NTG", keyword_handler (0, def_value_one, &array_dimens[0]));
      keyword_manager->register_keyword ("MULTX", keyword_handler (0, def_value_one, &array_dimens[0]));
      keyword_manager->register_keyword ("MULTY", keyword_handler (0, def_value_one, &array_dimens[0]));
      keyword_manager->register_keyword ("MULTZ", keyword_handler (0, def_value_one, &array_dimens[0]));
      keyword_manager->register_keyword ("MULTPV", keyword_handler (0, def_value_zero, &array_dimens[0]));
    }  
  
  void smesh_keywords::DIMENS_handler(const std::string &keyword, keyword_params_t &params)
    {
      KH_COMMON_VARIABLES_DEF
      t_long ndim = 0, nblock = 0;
      boost::array <t_long, 3> itmp;
      sp_idata_t idata (params.data, bs_dynamic_cast ());
      sp_reader_t reader (params.reader, bs_dynamic_cast ());
      
      if ((len = reader->read_array (keyword, itmp)) != 3)
        {
          out_s << "Error in " << reader->get_prefix() << ": not enough valid arguments for keyword " << keyword;
          //BOSERR << priority(section::read_data,level::err) << out_s << bs_end;
          KH_ASSERT_EXCEPTION
        }

      if (itmp[0] <= 0 || itmp[1] <= 0 || itmp[2] <= 0)
        {
          out_s << "Error in " << reader->get_prefix() << ": all dimensions in " << keyword << " must be positive";
          //BOSERR << priority(section::read_data,level::err) << out_s << bs_end;
          KH_ASSERT_EXCEPTION
        }

      idata->props->set_i ("nx", itmp[0]);
      idata->props->set_i ("ny", itmp[1]);
      idata->props->set_i ("nz", itmp[2]);
      
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
