#include "bs_bos_core_data_storage_stdafx.h"

#include "well.h"
#include "main_def.h"


 
namespace blue_sky
  {
  
  well::well (bs_type_ctor_param /*param*/)
  {
      sp_prop = BS_KERNEL.create_object ("prop");
      if (!sp_prop)
        {
          bs_throw_exception ("Type (prop) not refractered");
        }
  }

  well::well(const well& src): bs_refcounter ()
  {
    *this = src;
  }
  
  void 
  well::add_branch (const std::string &branch_name, sp_branch_t branch)
  {
     branches.insert (pair_t(branch_name, branch));
  }
     
  well::list_t 
  well::get_branch_names () const
  {
    list_t branch_names;
    map_t::const_iterator it, b, e;
    
    b = branches.begin();
    e = branches.end();
    
    for (it = b; it != e; it++)
      {
        branch_names.push_back (it->first);
      }
    
    return branch_names;
  }
  
  
  well::sp_branch_t 
  well::get_branch (const std::string &branch_name)
  {
    map_t::iterator i = branches.find(branch_name);
    if (i != branches.end ())
      return i->second;
    else
      return sp_branch_t (0);
  }

#ifdef BSPY_EXPORTING_PLUGIN
  std::string 
  well::py_str () const
    {
      std::stringstream s;
      s << sp_prop->py_str () << "\n";
      return s.str ();
    }
#endif //BSPY_EXPORTING_PLUGIN
  
  //bs stuff
  BLUE_SKY_TYPE_STD_CREATE(well)
  BLUE_SKY_TYPE_STD_COPY(well)

  BLUE_SKY_TYPE_IMPL(well, objbase, "well", "BOS_Core well class", "BOS_Core well class")
}//ns bs
