/** 
 * @file sql_well.h
 * @brief implementation of #well_pool_iface on sqlite database
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-29
 */
#ifndef SQL_WELL_YSG17OI0

#define SQL_WELL_YSG17OI0


#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <sqlite3.h>

#include "well_pool_iface.h"

namespace blue_sky
{
  
  class BS_API_PLUGIN sql_well : public well_pool_iface
    {
    
    public: 
      typedef std::list<std::string>                    list_t;
      typedef BS_SP (prop_iface)                        sp_prop_t;
      typedef BS_SP (table_iface)                       sp_table_t;
      typedef BS_SP (gis_iface)                         sp_gis_t;
      typedef BS_SP (traj_iface)                        sp_traj_t;

      // ------------------------------------
      // METHODS
      // ------------------------------------
    public:
      // destructor
      virtual ~sql_well ()
        {}
      // ------------------------------------
      // INTERFACE METHODS
      // ------------------------------------

      /** 
       * @brief return SP to the property 
       */
      //virtual sp_prop_t get_prop ()
      //  {
      //    return sp_prop;
      //  }

      /** 
       * @brief this method should be call first
       * 
       * @param file  -- <INPUT> data base file
       * 
       * @return 0 if success
       */
      virtual int open_db (const std::string &file); 

      /** 
       * @brief close database connection
       */
      virtual void close_db (); 


      /** 
       * @brief create internal structure of the database
       * 
       * @return 0 if success
       */
      virtual int create_db_struct ();

      virtual void fill_db ();

      // wels
      /** 
       * @brief add well to the storage
       * 
       * @param well_name -- <INPUT> name for the new well
       * 
       * @return 0 if OK
       */
      virtual int add_well (const std::string &well_name);
      virtual list_t get_well_names () const;
      //virtual int set_well_param (const std::string &wname, double date, const std::string param, double value);
      //virtual double get_well_param (const std::string &wname, double date, const std::string param);
      
      // branches
      //virtual list_t get_branches_names (const std::string &well_name) const;
      //virtual int add_branch (const std::string &wname, const std::string &branch, 
      //                        t_double md, const std::string &parent);
      //virtual int add_branch_prop (const std::string &wname, const std::string &branch,
      //                             sp_table_t tbl);
      virtual int add_branch_gis (const std::string &wname, const std::string &branch,
                                  sp_gis_t g);
      virtual int add_branch_traj (const std::string &wname, const std::string &branch,
                                   sp_traj_t t);
      //virtual sp_table_t get_branch_prop (const std::string &wname, const std::string &branch) const;
      virtual sp_gis_t get_branch_gis (const std::string &wname, const std::string &branch) const;
      virtual sp_traj_t get_branch_traj (const std::string &wname, const std::string &branch) const;
      //virtual void remove_branch (const std::string &wname, const std::string &branch);
    public:
#ifdef BSPY_EXPORTING_PLUGIN
      /** 
       * @brief python print wrapper
       * 
       * @return return table description
       */
      virtual std::string py_str () const;

#endif //BSPY_EXPORTING_PLUGIN
      
    protected:
      int create_db (sqlite3 *db_in);
      int insert_or_update (const std::string &select_sql,
                            const std::string &insert_sql,
                            const std::string &update_sql);

      // ------------------------------
      // VARIABLES
      // ------------------------------
    protected:
      //sp_prop_t sp_prop;        //!< ptoperties pointer
      std::string   file_name;          //!< database filename
      sqlite3       *db;                //!< database pointer

      BLUE_SKY_TYPE_DECL (sql_well);
    };

}; //end of blue_sky namespace

#endif /* end of include guard: sql_well_YSG17OI0 */
