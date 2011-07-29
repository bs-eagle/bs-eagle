#ifndef BS_WELL_H
#define BS_WELL_H
/**
 * \file well.h
 * \brief interface for well class
 * \author Mark Khait
 * \date 29.07.2011
 * */
 
#include "well_iface.h"
#include "well_branch_iface.h"
 
 
namespace blue_sky {

  class gis;
  class perforation;
  class fracture;
  
  class BS_API_PLUGIN well: public well_iface
    {
     
    public:

      //METHODS
     ~well ();
     
     void add_branch (std::string branch_name, sp_branch_iface branch);
     
     std::list <std::string> get_branch_names ();
     
     sp_branch_iface get_branch (std::string branch_name);
     
    public:
      BLUE_SKY_TYPE_DECL (well)

    public:
    
      std::map <std::string, sp_branch_iface> branches;
    
    };

} //ns bs
#endif //BS_WELL_H
