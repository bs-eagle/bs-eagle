#include "bs_bos_core_data_storage_stdafx.h"

#include "well_branch.h"
//#include "main_def.h"


 
namespace blue_sky
  {
  
  well_branch::~well_branch ()
  {

  }

  
  well_branch::well_branch(bs_type_ctor_param /*param*/)
  {
    sp_prop = BS_KERNEL.create_object ("prop");
    if (!sp_prop)
      {
        bs_throw_exception ("Type (prop) not registered");
      }
    sp_prop->add_property_s ("", "parent", L"name of parent branch");
    sp_prop->add_property_f (0, "md", L"measured depth of parent branch point");
  }

  well_branch::well_branch(const well_branch& src) 
    :bs_refcounter ()
  {
    *this = src;
  }
  

  
  //bs stuff
  BLUE_SKY_TYPE_STD_CREATE(well_branch)
  BLUE_SKY_TYPE_STD_COPY(well_branch)

  BLUE_SKY_TYPE_IMPL(well_branch, objbase, "well_branch", "BOS_Core well_branch class", "BOS_Core well_branch class")
}//ns bs
