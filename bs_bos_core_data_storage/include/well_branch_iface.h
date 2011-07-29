#ifndef BS_WELL_BRANCH_IFACE_H
#define BS_WELL_BRANCH_IFACE_H
/**
 * \file well_branch_iface.h
 * \brief interface for well_branch class represents well branch element 
 * \author Mark Khait
 * \date 29.07.2011
 * */
 
#include "conf.h"
#include "prop_iface.h"
#include "gis_iface.h"
#include "traj_iface.h"
 
 
namespace blue_sky {

  class gis;
  class perforation;
  class fracture;
  
  class BS_API_PLUGIN well_branch_iface: public objbase
    {
     public:
     
      typedef BS_SP (prop_iface)             sp_prop_iface;
      typedef BS_SP (gis_iface)              sp_gis_iface;
      typedef BS_SP (traj_iface)              sp_traj_iface;
      
    public:

      //METHODS
     virtual ~well_branch_iface () {};
    
     virtual sp_prop_iface get_prop () = 0;
     
     virtual sp_gis_iface get_gis () = 0;
     
     virtual sp_traj_iface get_traj () = 0;
     
     virtual void set_gis (sp_gis_iface new_gis) = 0;
     
     virtual void set_traj (sp_traj_iface new_traj) = 0;
      

    public:
    
    };

} //ns bs
#endif //BS_WELL_BRANCH_H
