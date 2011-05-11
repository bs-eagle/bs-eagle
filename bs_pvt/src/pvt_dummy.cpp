/**
 * \file pvt_dummy.h
 * \brief dummy holder for PVT data to represent it ro GUI tree
 * \author Mark Khait
 * \date 10.05.2010
 * */

#include "bs_pvt_stdafx.h"
#include "pvt_dummy.h"


namespace blue_sky {

    
    
  pvt_dummy::pvt_dummy (bs_type_ctor_param param /* = NULL */)
  {
    pvt_data = BS_KERNEL.create_object ("table");
  }
  
  
  pvt_dummy::pvt_dummy (const pvt_dummy& s)
  : bs_refcounter (s), pvt_dummy_iface (s)
  {
    BS_ASSERT (false && "NOT IMPL YET");
  }
  
  void
  pvt_dummy::init (int pvt_type_)
  {
    pvt_type = pvt_type_;
  }
  
      
  BS_SP ( table_iface )
  pvt_dummy::get_table () const
  {
    return pvt_data;
  }
  
  int
  pvt_dummy::get_pvt_type () const
  {
    return pvt_type;
  }
  
  
  
   //bs stuff
  BLUE_SKY_TYPE_STD_CREATE(pvt_dummy)
  BLUE_SKY_TYPE_STD_COPY(pvt_dummy)

  BLUE_SKY_TYPE_IMPL(pvt_dummy, objbase, "pvt_dummy", "pvt_dummy", "pvt_dummy")
}