#include "bs_bos_core_data_storage_stdafx.h"

#include "well.h"
#include "main_def.h"


 
namespace blue_sky
  {
  
  well::~well ()
  {

  }

  
  well::well(bs_type_ctor_param /*param*/)
  {
  }

  well::well(const well& src)
  {
    *this = src;
  }
  
  void well::add_branch (std::string branch_name, sp_branch_iface branch)
  {
    branches.insert(std::pair<std::string, sp_branch_iface>(branch_name, branch));
  }
     
  std::list <std::string> 
  well::get_branch_names ()
  {
    std::list <std::string> branch_names;
    std::map <std::string, sp_branch_iface>::iterator it, b, e;
    
    b = branches.begin();
    e = branches.end();
    
    for (it = b; it != e; it++)
      {
        branch_names.push_back(it->first);
      }
    
    return branch_names;
  }
  
  
  well::sp_branch_iface
  well::get_branch (std::string branch_name)
  {
    return branches.find(branch_name)->second;
  }

  
  //bs stuff
  BLUE_SKY_TYPE_STD_CREATE(well)
  BLUE_SKY_TYPE_STD_COPY(well)

  BLUE_SKY_TYPE_IMPL(well, objbase, "well", "BOS_Core well class", "BOS_Core well class")
}//ns bs
