#include "well.h"
#include "bs_misc.h"
//#include "main_def.h"


 
namespace blue_sky
  {
  
  well_obj::well_obj (bs_type_ctor_param /*param*/)
  {
      sp_prop = BS_KERNEL.create_object ("prop");
      if (!sp_prop)
        {
          bs_throw_exception ("Type (prop) not refractered");
        }
  }

  well_obj::well_obj(const well_obj& src): bs_refcounter ()
  {
    *this = src;
  }
  
  void 
  well_obj::add_branch (const std::string &branch_name, sp_branch_t branch)
  {
     branches.insert (pair_t(branch_name, branch));
  }
     
  well_obj::list_t 
  well_obj::get_branch_names () const
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
  
  
  well_obj::sp_branch_t 
  well_obj::get_branch (const std::string &branch_name)
  {
    map_t::iterator i = branches.find(branch_name);
    if (i != branches.end ())
      return i->second;
    else
      return sp_branch_t (0);
  }

#ifdef BSPY_EXPORTING_PLUGIN
  std::string 
  well_obj::py_str () const
    {
      std::stringstream s;
      s << wstr2str(sp_prop->py_str ()) << "\n";
      return s.str ();
    }
#endif //BSPY_EXPORTING_PLUGIN
  
  //bs stuff
  BLUE_SKY_TYPE_STD_CREATE(well_obj)
  BLUE_SKY_TYPE_STD_COPY(well_obj)

  BLUE_SKY_TYPE_IMPL(well_obj, well_obj_iface, "well_obj", "BOS_Core well_obj class", "BOS_Core well_obj class")
}//ns bs
