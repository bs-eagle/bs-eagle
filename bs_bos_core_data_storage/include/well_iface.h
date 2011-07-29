#ifndef BS_WELL_IFACE_H
#define BS_WELL_IFACE_H
/**
 * \file well_iface.h
 * \brief interface for well class
 * \author Mark Khait
 * \date 29.07.2011
 * */
 
#include "well_branch_iface.h" 
 
namespace blue_sky {

  
  class BS_API_PLUGIN well_iface: public objbase
    {
    public:
      typedef BS_SP (well_branch_iface)             sp_branch_iface;
      
    public:  

      //METHODS
     virtual ~well_iface () {};
     
     virtual void add_branch (std::string branch_name, sp_branch_iface branch) = 0;
     
     virtual std::list <std::string> get_branch_names () = 0;
     
     virtual sp_branch_iface get_branch (std::string branch_name) = 0;

    public:
    
    };

} //ns bs
#endif //BS_WELL_IFACE_H
