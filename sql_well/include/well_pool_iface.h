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

#ifdef BSPY_EXPORTING_PLUGIN
#include <boost/python/list.hpp>
#endif //BSPY_EXPORTING_PLUGIN

#include "bs_object_base.h"
#include "conf.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include BS_STOP_PLUGIN_IMPORT ()

#include "prop_iface.h"
#include "table_iface.h"
#include "gis_iface.h"
#include "traj_iface.h"
#include "pool_iface.h"


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
      typedef BS_SP (h5_pool_iface)           					sp_pool_t;

      static const int CTRL_P_BHP = 1;
      static const int CTRL_P_LRATE = 2;
      static const int CTRL_P_ORATE = 3;
      static const int CTRL_P_WRATE = 4;
      static const int CTRL_P_GRATE = 5;

      static const int CTRL_I_BHP = -1;
      static const int CTRL_I_WRATE = -2;
      static const int CTRL_I_GRATE = -3;
      static const int CTRL_I_ORATE = -4;
      static const int STATUS_SHUT = 0;
      static const int STATUS_CLOSE = 1;
      static const int STATUS_OPEN = 2;

      static const int STATUS_CON_SHUT = 0;
      static const int STATUS_CON_OPEN = 1;

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


      virtual int prepare_sql (const std::string &sql) = 0;
      virtual int step_sql () = 0;
      virtual int finalize_sql () = 0;
      virtual t_int get_sql_int (t_int col) = 0;
      virtual t_double get_sql_real (t_int col) = 0;
      virtual bool get_sql_bool (t_int col) = 0;
      virtual std::string get_sql_str (t_int col) = 0;
	  virtual bool get_sql_exist (t_int col) = 0;
      virtual int exec_sql (const std::string &sql) = 0;
      virtual int insert_or_update (const std::string &select_sql,
                                    const std::string &insert_sql,
                                    const std::string &update_sql) = 0;
      virtual spv_double get_table (const std::string &table_name, boost::python::list &table_columns, const std::string &filter) = 0;

      /** 
       * @brief create internal structure of the database
       * 
       * @return 0 if success
       */
      virtual int create_db_struct () = 0;

      virtual void fill_db () = 0;

      /** 
       * @brief read from ascii file in new format
       * 
       * @param fname -- <INPUT> input file name
       * @param starting_date -- <INPUT> starting date
       *
       * @return 0 if success
       */
      virtual int read_from_ascii_file (const std::string &fname, 
                                        double starting_date) = 0;

      /** 
       * @brief Save all data from db to BOS ascii format
       * 
       * @param fname   -- <INPUT> file name
       * 
       * @return 0 if success
       */
      virtual int save_to_bos_ascii_file (const std::string &fname, sp_pool_t pool) = 0;
#ifdef BSPY_EXPORTING_PLUGIN
      virtual boost::python::list d2date (double d) const = 0;
      virtual double date2d (int year, int month, int day, int hour, int minute, int second) const = 0;
      virtual std::string d2str (double d) const = 0;
      virtual std::string t2str (double d) const = 0;

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
