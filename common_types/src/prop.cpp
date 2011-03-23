/** 
 * @file prop.cpp
 * @brief 
 * @author Oleg Borschuk
 * @date 2009-08-28
 */
#include "bs_kernel.h"
#include "prop.h"

using namespace std;
using namespace boost::python;


namespace blue_sky
{
  prop::prop (bs_type_ctor_param) 
        : prop_iface ()
    {
    }
  prop::prop (const prop& rhs) 
        : bs_refcounter (), prop_iface ()
    {
      *this = rhs;
    }

/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (prop);
  BLUE_SKY_TYPE_STD_COPY (prop);

  BLUE_SKY_TYPE_IMPL(prop, prop_iface, "prop", "Property storage", "realization of property storage");

}  // blue_sky namespace

