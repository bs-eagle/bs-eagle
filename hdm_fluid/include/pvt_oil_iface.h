/** 
 * @file pvt_oil_iface.h
 * @brief OIL PVT interface
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-03-25
 */
#ifndef PVT_OIL_IFACE_5STFVJQF

#define PVT_OIL_IFACE_5STFVJQF

#include "pvt_iface.h"
#include "table_iface.h"

namespace blue_sky
{
  /** 
   * @brief PVT OIL interface
   */
  class pvt_oil_iface: public pvt_iface
  {
    public:
      typedef smart_ptr<prop_iface, true>       sp_prop_t;
      typedef smart_ptr<table_iface, true>      sp_table_t;

      /** 
       * @brief destructor
       */
       virtual ~pvt_oil_iface () {}
  };
};

#endif /* end of include guard: PVT_OIL_IFACE_5STFVJQF */
