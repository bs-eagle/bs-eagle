#ifndef BS_WELL_H
#define BS_WELL_H
/**
 * \file well.h
 * \brief interface for well class
 * \author Mark Khait
 * \date 29.07.2011
 * */
#include <map>
#include <string>

#include "well_iface.h"
#include "well_branch_iface.h"
 
 
namespace blue_sky {

  class gis;
  class perforation;
  class fracture;
  
  class BS_API_PLUGIN well_obj: public well_obj_iface
    {
     
    public:
      typedef BS_SP (well_branch_iface)                 sp_branch_t;
      typedef std::map <std::string, sp_branch_t>       map_t;
      typedef std::pair<std::string, sp_branch_t>       pair_t;
      typedef std::list<std::string>                    list_t;
      typedef BS_SP (prop_iface)                        sp_prop_t;

    public:
      
      ~well_obj ()
        {}

     /** 
      * @brief add new branch to the well
      * 
      * @param branch_name  -- <INPUT> new branch name
      * @param branch       -- <INPUT> new branch
      */
     virtual void add_branch (const std::string &branch_name, sp_branch_t branch);
     
     /** 
      * @brief return list of branch names
      */
     virtual list_t get_branch_names () const;
     
     /** 
      * @brief return branch by name
      * 
      * @param branch_name  -- <INPUT> branch name
      * 
      * @return smatr pointer ti branch
      */
     virtual sp_branch_t get_branch (const std::string &branch_name);

      /** 
       * @brief return SP to the property 
       */
      virtual sp_prop_t get_prop ()
        {
          return sp_prop;
        }

#ifdef BSPY_EXPORTING_PLUGIN
      /** 
       * @brief python print wrapper
       * 
       * @return return table description
       */
      virtual std::string py_str () const;

#endif //BSPY_EXPORTING_PLUGIN
      //METHODS
    public:
      BLUE_SKY_TYPE_DECL (well_obj)

    public:
      sp_prop_t sp_prop; 
      map_t branches;
    };

} //ns bs
#endif //BS_WELL_H
