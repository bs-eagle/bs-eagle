/** 
 * @file prop.cpp
 * @brief 
 * @author Oleg Borschuk
 * @date 2009-08-28
 */
#ifdef BSPY_EXPORTING_PLUGIN
#include <boost/python.hpp>
#endif

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
        : bs_refcounter (rhs), prop_iface ()
    {
      *this = rhs;
    }

#ifdef BSPY_EXPORTING_PLUGIN
  std::string 
  prop::to_str () const
    {
      std::ostringstream oss;

      boost::archive::text_oarchive oar(oss);

      save (oar);
      return oss.str ();
    }
  void 
  prop::from_str (const std::string &s)
    {
      std::istringstream iss;

      iss.str (s);
      boost::archive::text_iarchive iar(iss);
      load (iar);
    }
#endif //BSPY_EXPORTING_PLUGIN
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (prop);
  BLUE_SKY_TYPE_STD_COPY (prop);

  BLUE_SKY_TYPE_IMPL(prop, prop_iface, "prop", "Property storage", "realization of property storage");

}  // blue_sky namespace

