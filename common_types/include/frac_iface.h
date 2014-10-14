/** 
 * @file frac_iface.h
 * @brief interface for the fracture
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-29
 */
#ifndef FRAC_IFACE_LA5BNEMD

#define FRAC_IFACE_LA5BNEMD


#include <string>

#include "bs_object_base.h"
#include "conf.h"
#include "prop_iface.h"

namespace blue_sky
{
class BS_API_PLUGIN frac_iface : public objbase
  {
    public:
      typedef BS_SP (prop_iface)                      sp_prop_t;

    public:
      /** 
       * @brief destructor
       */
      virtual ~frac_iface ()
        {}

      /** 
       * @brief return SP to the property 
       */
      virtual sp_prop_t get_prop () = 0;

#ifdef BSPY_EXPORTING_PLUGIN
      /** 
       * @brief python print wrapper
       * 
       * @return return table description
       */
      virtual std::string py_str () const = 0;

#endif //BSPY_EXPORTING_PLUGIN
};

}  // end of bluesky name space

#endif /* end of include guard: FRAC_IFACE_LA5BNEMD */
