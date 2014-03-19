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

#include "bos_reader_iface.h"
#include "well_pool_iface.h"

#include <boost/serialization/access.hpp>
#include "bs_serialize_decl.h"

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
      typedef BS_SP (h5_pool_iface)           			sp_pool_t;
      typedef BS_SP (bos_reader_iface)                  sp_reader_t;
      typedef BS_SP (dt_tools_iface)                    sp_dt_t;

      // ------------------------------------
      // METHODS
      // ------------------------------------
    public:
      // destructor
      virtual ~sql_well ();
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
                                  sp_gis_t g, std::string wlog_name = "", uint wlog_type = 0);
      virtual int add_branch_traj (const std::string &wname, const std::string &branch,
                                   sp_traj_t t);
      //virtual sp_table_t get_branch_prop (const std::string &wname, const std::string &branch) const;
      virtual sp_gis_t get_branch_gis (const std::string &wname, const std::string &branch,
                                   std::string wlog_name = "", uint wlog_type = 0);
      virtual sp_traj_t get_branch_traj (const std::string &wname, const std::string &branch);
      //virtual void remove_branch (const std::string &wname, const std::string &branch);

      virtual int update_branch_traj (const std::string &wname, const std::string &branch, sp_traj_t t);

      virtual void backup_to_file (const std::string &filename);

      virtual int prepare_sql (const std::string &sql);
      virtual int step_sql ();
      virtual int finalize_sql ();
      virtual t_int get_sql_int (t_int col);
      virtual t_double get_sql_real (t_int col);
      virtual bool get_sql_bool (t_int col);
      virtual std::string get_sql_str (t_int col);
	  virtual bool get_sql_exist (t_int col);
      virtual int exec_sql (const std::string &sql);
      virtual int exec_sql_and_return_rowid (const std::string &sql);
      virtual int merge_with_db (const std::string &dbname);
      virtual int insert_or_update (const std::string &select_sql,
                                    const std::string &insert_sql,
                                    const std::string &update_sql);
      /**
       * @brief get data from SQL table
       * @param table_name     -- <INPUT> name of table
       * @param table_columns  -- <INPUT> list of column names
       * @param filter         -- <INPUT> specification for SQL.SELECT
       * @return smart pointer to 2D array
       */
      virtual spv_double get_table (const std::string &table_name, boost::python::list &table_columns, const std::string &filter);

      /**
       * @brief read from ascii file in new format
       *
       * @param fname -- <INPUT> input file name
       *
       * @return 0 if success
       */
      virtual int read_from_ascii_file (const std::string &fname, double starting_date);

      /**
       * @brief Save all data from db to BOS ascii format
       *
       * @param fname   -- <INPUT> file name
       *
       * @return 0 if success
       */
      virtual int save_to_bos_ascii_file (const std::string &fname, sp_pool_t pool, sp_prop_t prop);

      // return list of cutom well logs added via add_branch_gis() with nonempty
      // wlog_name parameter
      std::vector< std::string > get_wlog_names(const std::string &wname, const std::string &branch);

    public:
#ifdef BSPY_EXPORTING_PLUGIN
#if 0
      virtual boost::python::list d2date (double d) const
        {
          boost::python::list l;
          int day, month, year, hour, minute, second;
          double dd;
          dd = get_date_day_month_year (d, day, month, year);
          d2hms (dd, hour, minute, second);
          l.append (year);
          l.append (month);
          l.append (day);
          l.append (hour);
          l.append (minute);
          l.append (second);
          return l;
        }
      virtual double date2d (int year, int month, int day, int hour, int minute, int second) const
        {
          double tm;
          tm = ((double)hour + ((double)minute + (date_sim)second / 60.0) / 60.0) / 24.0;
          return ymd2d (year, month, day) + tm;
        }
      virtual std::string d2str (double d) const
        {
          char b[1024];
          int day, month, year;
          get_date_day_month_year (d, day, month, year);
          sprintf (b, "%02d.%02d.%04d", day, month, year);
          return std::string (b);
        }
      virtual std::string t2str (double d) const
        {
          char b[1024];
          int hour, minute, second;
          d2hms (d, hour, minute, second);
          sprintf (b, "Time: %02d:%02d:%02d", hour, minute, second);
          return std::string (b);
        }
#endif //0
      /**
       * @brief python print wrapper
       *
       * @return return table description
       */
      virtual std::string py_str () const;

#endif //BSPY_EXPORTING_PLUGIN

      std::string   file_name;          //!< database filename

    protected:
      int create_db (sqlite3 *db_in);

      int read_date_and_time (char *buf, char **next_start, double *dd);
      int read_w_spec (char *buf);
      int read_w_branch_f (char *buf);
      int read_w_comp (char *buf, double d);
      int read_w_frac (char *buf, double d);
      int read_w_prod (char *buf, double d);
      int read_w_inj (char *buf, double d);


      // ------------------------------
      // VARIABLES
      // ------------------------------
    protected:
      //sp_prop_t sp_prop;        //!< ptoperties pointer
      sqlite3       *db;                //!< database pointer
      sqlite3_stmt  *stmp_sql;
      sp_reader_t   fr_file;            //!< read from ascii helper

      BLUE_SKY_TYPE_DECL (sql_well);

      friend class boost::serialization::access;
      friend class blue_sky::bs_serialize;
    };

}; //end of blue_sky namespace

#endif /* end of include guard: sql_well_YSG17OI0 */
