/** 
 * @file cljp2_coarse.cpp
 * @brief implementation of Cleary-Luby-Jones-Plassmann grid coasering method
 * @author 
 * @version 
 * @date 2010-03-10
 */
#include "cljp2_coarse.h"

#include "coarse_tools.h"


namespace blue_sky
{
  cljp2_coarse::cljp2_coarse (bs_type_ctor_param) 
                : amg_coarse_iface(),
                  sp_graph (BS_KERNEL.create_object (v_long::bs_type ()))
    {
    }
  cljp2_coarse::cljp2_coarse (const this_t & /*src*/) : bs_refcounter (),
                sp_graph (BS_KERNEL.create_object (v_long::bs_type ()))
     {
     }
  
  t_long cljp2_coarse::build (sp_bcsr_t  /*matrix*/, 
                              spv_double /*meassure_array*/, 
                              spv_long /*cf_markers*/,
                              spv_long /*s_markers*/)
    {
      return 0;
    }
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (cljp2_coarse);
  BLUE_SKY_TYPE_STD_COPY (cljp2_coarse);

  BLUE_SKY_TYPE_IMPL (cljp2_coarse, amg_coarse_iface, "cljp2_coarse", "CLJP coarsering builder class", "Realization of CLJP coarsering");
}  // blue_sky namespace
