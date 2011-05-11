/**
 * \file scal_2p_dummy.h
 * \brief dummy holder for SCAL data to represent it ro GUI tree
 * \author Mark Khait
 * \date 10.05.2010
 * */

#include "bs_scal_stdafx.h"
#include "scal_2p_dummy.h"


namespace blue_sky {

    
    
  scal_2p_dummy::scal_2p_dummy (bs_type_ctor_param param /* = NULL */)
  {
    scal_data = BS_KERNEL.create_object ("table");
  }
  
  scal_3p_dummy::scal_3p_dummy (bs_type_ctor_param param /* = NULL */)
  {
    scal_data1 = BS_KERNEL.create_object ("table");
    scal_data2 = BS_KERNEL.create_object ("table");
  }
  
  scal_2p_dummy::scal_2p_dummy (const scal_2p_dummy& s)
  : bs_refcounter (s), scal_dummy_iface (s)
  {
    BS_ASSERT (false && "NOT IMPL YET");
  }
  
  scal_3p_dummy::scal_3p_dummy (const scal_3p_dummy& s)
  : bs_refcounter (s), scal_dummy_iface (s)
  {
    BS_ASSERT (false && "NOT IMPL YET");
  }
    
  std::pair <BS_SP( table_iface), BS_SP( table_iface)>
  scal_2p_dummy::get_table () const
  {
    return std::pair<BS_SP( table_iface), BS_SP( table_iface)> (scal_data, 0);
  }
  
  std::pair <BS_SP( table_iface), BS_SP( table_iface)>
  scal_3p_dummy::get_table () const
  {
    return std::pair<BS_SP( table_iface), BS_SP( table_iface)> (scal_data1, scal_data2);
  }
  
   //bs stuff
  BLUE_SKY_TYPE_STD_CREATE(scal_2p_dummy)
  BLUE_SKY_TYPE_STD_COPY(scal_2p_dummy)

  BLUE_SKY_TYPE_IMPL(scal_2p_dummy, objbase, "scal_2p_dummy", "scal_2p_dummy", "scal_2p_dummy")
  
  BLUE_SKY_TYPE_STD_CREATE(scal_3p_dummy)
  BLUE_SKY_TYPE_STD_COPY(scal_3p_dummy)

  BLUE_SKY_TYPE_IMPL(scal_3p_dummy, objbase, "scal_3p_dummy", "scal_3p_dummy", "scal_3p_dummy")

}