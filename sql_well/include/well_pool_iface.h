/** 
 * @file well_pool_iface.h
 * @brief interface for the well storage
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-29
 */
#ifndef WELL_POOL_IFACE_LA5BNEMD

#define WELL_POOL_IFACE_LA5BNEMD


#include <string>

#include "bs_object_base.h"
#include "conf.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include BS_STOP_PLUGIN_IMPORT ()


#include "prop_iface.h"
#include "table_iface.h"
#include "gis_iface.h"
#include "traj_iface.h"


namespace blue_sky
{
class well_pool_iface : public objbase
  {
    public:
      typedef std::list<std::string>                    list_t;
      typedef BS_SP (prop_iface)                        sp_prop_t;
      typedef BS_SP (table_iface)                       sp_table_t;
      typedef BS_SP (gis_iface)                         sp_gis_t;
      typedef BS_SP (traj_iface)                        sp_traj_t;

    public:
      /** 
       * @brief destructor
       */
      virtual ~well_pool_iface ()
        {}

      /** 
       * @brief return SP to the property 
       */
      //virtual sp_prop_t get_prop () = 0;

      /** 
       * @brief this method should be call first
       * 
       * @param file  -- <INPUT> data base file
       * 
       * @return 0 if success
       */
      virtual int open_db (const std::string &file) = 0; 

      /** 
       * @brief close database connection
       */
      virtual void close_db () = 0; 


      // wels
      /** 
       * @brief add well to the storage
       * 
       * @param well_name -- <INPUT> name for the new well
       * 
       * @return 0 if OK
       */
      virtual int add_well (const std::string &well_name) = 0;
      virtual list_t get_well_names () const = 0;
      //virtual int set_well_param (const std::string &wname, double date, const std::string param, double value) = 0;
      //virtual double get_well_param (const std::string &wname, double date, const std::string param) = 0;
      
      // branches
      //virtual list_t get_branches_names (const std::string &well_name) const = 0;
      //virtual int add_branch (const std::string &wname, const std::string &branch, 
      //                        t_double md, const std::string &parent) = 0;
      //virtual int add_branch_prop (const std::string &wname, const std::string &branch,
      //                             sp_table_t tbl) = 0;
      virtual int add_branch_gis (const std::string &wname, const std::string &branch,
                                  sp_gis_t g) = 0;
      virtual int add_branch_traj (const std::string &wname, const std::string &branch,
                                   sp_traj_t t) = 0;
      //virtual sp_table_t get_branch_prop (const std::string &wname, const std::string &branch) const = 0;
      virtual sp_gis_t get_branch_gis (const std::string &wname, const std::string &branch) const = 0;
      virtual sp_traj_t get_branch_traj (const std::string &wname, const std::string &branch) const = 0;
      //virtual void remove_branch (const std::string &wname, const std::string &branch) = 0;


      /** 
       * @brief create internal structure of the database
       * 
       * @return 0 if success
       */
      virtual int create_db_struct () = 0;

      virtual void fill_db () = 0;
#ifdef BSPY_EXPORTING_PLUGIN
      /** 
       * @brief python print wrapper
       * 
       * @return return table description
       */
      virtual std::string py_str () const = 0;

#endif //BSPY_EXPORTING_PLUGIN
};

}  // end of bluesky name space

#endif /* end of include guard: well_pool_IFACE_LA5BNEMD */
