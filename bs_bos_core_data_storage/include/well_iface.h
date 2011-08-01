#ifndef BS_WELL_IFACE_H
#define BS_WELL_IFACE_H
/**
 * \file well_iface.h
 * \brief interface for well class
 * \author Mark Khait
 * \date 29.07.2011
 * */
#include <string> 
#include "bs_object_base.h"
#include "conf.h"
#include "well_branch_iface.h" 
#include "prop_iface.h"
 
namespace blue_sky {

  
  class well_iface: public objbase
    {
    public:
      typedef BS_SP (well_branch_iface)                 sp_branch_t;
      typedef std::map <std::string, sp_branch_t>       map_t;
      typedef std::list<std::string>                    list_t;
      typedef BS_SP (prop_iface)                        sp_prop_t;

      
    public:  

      //METHODS
      /** 
       * @brief destructor
       */
     virtual ~well_iface () {};
     
     /** 
      * @brief add new branch to the well
      * 
      * @param branch_name  -- <INPUT> new branch name
      * @param branch       -- <INPUT> new branch
      */
     virtual void add_branch (const std::string &branch_name, sp_branch_t branch) = 0;
     
     /** 
      * @brief return list of branch names
      */
     virtual list_t get_branch_names () const = 0;
     
     /** 
      * @brief return branch by name
      * 
      * @param branch_name  -- <INPUT> branch name
      * 
      * @return smatr pointer ti branch
      */
     virtual sp_branch_t get_branch (const std::string &branch_name) = 0;

      /** 
       * @brief return SP to the property 
       */
      virtual sp_prop_t get_prop () = 0;

#ifdef BSPY_EXPORTING_PLUGIN
      /** 
       * @brief python print wrapper
       * 
       * @return return table description
       */
      virtual std::string py_str () const = 0;

#endif //BSPY_EXPORTING_PLUGIN

    };

} //ns bs
#endif //BS_WELL_IFACE_H
