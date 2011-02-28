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
  template <class t_double, class i_type_t, class s_type_t, class b_type_t>
  prop<t_double, i_type_t, s_type_t, b_type_t>::prop (bs_type_ctor_param) 
        : prop_iface<t_double, i_type_t, s_type_t, b_type_t> ()
    {
    }
  template <class t_double, class i_type_t, class s_type_t, class b_type_t>
  prop<t_double, i_type_t, s_type_t, b_type_t>::prop (const prop& rhs) 
        : bs_refcounter (), prop_iface<t_double, i_type_t, s_type_t, b_type_t> ()
    {
      *this = rhs;
    }

/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE_T_DEF(prop, (class)(class)(class)(class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF(prop, (class)(class)(class)(class));

  BLUE_SKY_TYPE_IMPL_T_EXT(4, (prop<float, int, std::string, bool >), 4,  (prop_iface <float,  int, std::string, bool>), "prop", "Property storage", "realization of property storage", false);

}  // blue_sky namespace

