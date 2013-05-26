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
#include "bs_serialize.h"

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

  std::string 
  prop::to_str () const
    {
      return serialize_to_str_indirect< prop, prop_iface >(this);
    }
  void 
  prop::from_str (const std::string &s)
    {
      smart_ptr< prop > pv = serialize_from_str_indirect< prop, prop_iface >(s);
      fp_impl = pv->fp_impl;
      i_impl = pv->i_impl;
      s_impl = pv->s_impl;
      b_impl = pv->b_impl;
    }
/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (prop);
  BLUE_SKY_TYPE_STD_COPY (prop);

  BLUE_SKY_TYPE_IMPL(prop, prop_iface, "prop", "Property storage", "realization of property storage");

}  // blue_sky namespace

