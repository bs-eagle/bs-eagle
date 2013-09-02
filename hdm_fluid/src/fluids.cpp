/** 
 * @file fluids.cpp
 * @brief fluid information storage
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-03-30
 */

#include "fluids.h"

using namespace std;
using namespace boost::python;


namespace blue_sky
{
  fluids::fluids (bs_type_ctor_param) 
    {
    }
  fluids::fluids (const fluids & /*src*/)
        : bs_refcounter ()
    {
    }


/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (fluids);
  BLUE_SKY_TYPE_STD_COPY (fluids);

  BLUE_SKY_TYPE_IMPL (fluids, fluids_iface, "fluids", "Fluids properties", "Realization of Fluids properties");
}  // blue_sky namespace
