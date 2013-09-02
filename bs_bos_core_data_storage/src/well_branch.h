#ifndef BS_WELL_BRANCH_H
#define BS_WELL_BRANCH_H
/**
 * \file well_branch.h
 * \brief well_branch class represents well branch element 
 * \author Mark Khait
 * \date 28.07.2011
 * */

#include "well_branch_iface.h" 
#include "conf.h"

 
 
namespace blue_sky {

  class gis;
  class perforation;
  class fracture;
  
  class BS_API_PLUGIN well_branch: public well_branch_iface
    {
    public:
    
      typedef BS_SP (prop_iface)             sp_prop_t;
      typedef BS_SP (gis_iface)              sp_gis_t;
      typedef BS_SP (traj_iface)             sp_traj_t;
      typedef std::map <t_float, perforation *> perforation_map;
      typedef std::map <t_float, fracture *> fracture_map;
      typedef BS_SP (prop_iface)             sp_prop_iface;
      
    public:

      //METHODS
      ~well_branch();
    
      sp_prop_t get_prop ()
        {
          return sp_prop;
        }
     
      sp_gis_t get_gis ()
        {
          return sp_gis;
        }
     
      sp_traj_t get_traj ()
        {
          return sp_traj;
        }
     
      void set_gis (sp_gis_t new_gis)
        {
          sp_gis = new_gis;
        }
      
      void set_traj (sp_traj_t new_traj)
        {
          sp_traj = new_traj;
        }
      
      std::string py_str() const
        {
          return "well_branch";
        }
     
    public:
      BLUE_SKY_TYPE_DECL (well_branch)

    public:
    
      perforation_map perfs;
      fracture_map fracs;
      
      sp_prop_t sp_prop;
      
      sp_gis_t sp_gis;
      sp_traj_t sp_traj;
    };

} //ns bs
#endif //BS_WELL_BRANCH_H
