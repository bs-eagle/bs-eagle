/**
 * @file sql_well.cpp
 * @brief implementation of frac storage
 * @author Oleg Borschuk
 * @version
 * @date 2011-07-29
 */

// d REAL NOT NULL REFERENCES dates(d) ON UPDATE CASCADE ON DELETE CASCADE,
//
// d REAL NOT NULL,

#include "bs_kernel.h"
#include "bs_misc.h"
#include "sql_well.h"

#include <iostream>
#include <stdio.h>
#include <ctype.h>
#include <boost/lexical_cast.hpp>
//#include <boost/format.hpp>

#ifdef BSPY_EXPORTING_PLUGIN

#include <boost/python.hpp>
using namespace boost::python;

#endif //BSPY_EXPORTING_PLUGIN

#define DB_FORMAT_PROP L"DB_format"

using namespace boost;

namespace blue_sky {

// hidden details
namespace {
static const char* spaces = " \n\r\t";

std::string trim(const std::string& ss) {
  std::string s = ss;
  while(s.size() > 0 && strchr(spaces, s[0]) != NULL)
   s.erase(s.begin());
  while(s.size() > 0 && strchr(spaces, s[s.size() - 1]) != NULL)
   s.erase(s.size() - 1);
  return s;
}

std::string to_upper(const std::string& s) {
  std::string us(s.size(), ' ');
  std::transform(s.begin(), s.end(), us.begin(), ::toupper);
  return us;
}

std::string to_lower(const std::string& s) {
  std::string us(s.size(), ' ');
  std::transform(s.begin(), s.end(), us.begin(), ::tolower);
  return us;
}

} // eof hidden namespace

  sql_well::sql_well (bs_type_ctor_param)
    {
      db = 0;
      stmp_sql = 0;
      fr_file = 0;

    }
  sql_well::sql_well (const sql_well& rhs)
        : bs_refcounter (), file_name(rhs.file_name), db(rhs.db)
    {
      //*this = rhs;
      // don't copy pending statement
      stmp_sql = 0;
      fr_file = 0;
    }
  sql_well::~sql_well ()
    {
      if (db)
        close_db ();
    }


  int
  sql_well::add_well (const std::string &well_name)
    {
      if (!db)
        return -2;
      if (stmp_sql)
        finalize_sql ();

      int rc = 0;
      char *zErrMsg = 0;
      char buf[2048];

      sprintf (buf, "INSERT INTO wells (name) VALUES('%s');",
               well_name.c_str ());
      rc = sqlite3_exec (db, buf, NULL, 0, &zErrMsg);
      if( rc != SQLITE_OK )
        {
          fprintf (stderr, "SQL error (add_well): %s\n", zErrMsg);
          sqlite3_free (zErrMsg);
        }
      return 0;
    }

  sql_well::list_t
  sql_well::get_well_names () const
    {
      list_t lst;
      if (!db)
        return lst;

      int rc = 0;
      //char *zErrMsg = 0;
      const char *ttt;
      sqlite3_stmt *stmp;
      std::string sql = "SELECT name FROM wells ORDER BY name ASC";
      rc = sqlite3_prepare_v2 (db, sql.c_str (), sql.length () + 1, &stmp, &ttt);
      if (rc)
        {
          fprintf (stderr, "Can't make select: %s\n", sqlite3_errmsg (db));
          return lst;
        }
      while (sqlite3_step (stmp) == SQLITE_ROW) // UPDATE
        {

          lst.push_back (std::string ((const char *)sqlite3_column_text (stmp, 0)));
        }
      sqlite3_finalize (stmp);
      return lst;
    }

  int
  sql_well::add_branch_gis (const std::string &wname, const std::string &branch,
                            sp_gis_t g, std::string wlog_name, uint wlog_type,
                            bool replace_existing
  ) {
      // helper
      // bind well log data blob and execute PREPARED sql statement
      struct dump_wlog_data {
        static int go(sql_well& sqw, const sp_gis_t& g) {
          // bind well log data
          const std::string log_data = g->to_str();
          if (sqlite3_bind_blob(sqw.stmp_sql, 1, &log_data.c_str()[0], log_data.size (), SQLITE_STATIC))
            {
              fprintf (stderr, "Can't make select: %s\n", sqlite3_errmsg (sqw.db));
              sqw.finalize_sql();
              return -3;
            }

          // exec query
          sqw.step_sql();
          sqw.finalize_sql();
          return 0;
        }
      };

      if (!db)
        return -1;
      if (stmp_sql)
        finalize_sql ();

      std::string q;
      int res;

      if(wlog_name.size() == 0) {
        // 0. Check wlog property indicating if it has been converted to new format before
        sp_prop_t log_prop = g->get_prop();
        std::vector< std::wstring > names = log_prop->get_names_i();
        // ensure DB format prop exists
        if(std::find(names.begin(), names.end(), DB_FORMAT_PROP) == names.end()) {
          log_prop->add_property_i(0, DB_FORMAT_PROP, L"");
        }
        // mark that we store logs in new format
        log_prop->set_i(DB_FORMAT_PROP, 1);

        // 1. convert old representation to new
        // split logs table into sequence of DEPT-LOG data tables
        sp_table_t log_data = g->get_table();
        names = log_data->get_col_names();
        // process only tables with >= 2 columns
        if(names.size() < 2)
          return res;

        // search DEPTH values
        // take first column by default
        ulong dept_idx = 0;
        for(ulong i = 0; i < names.size(); ++i) {
          const std::string cur_col = trim(to_lower(wstr2str(names[i], "utf-8")));
          if(cur_col.find("dept", 0, 4) != std::string::npos) {
            dept_idx = i;
            break;
          }
        }
        // extract depth vector
        spv_double dept_data = BS_KERNEL.create_object(v_double::bs_type());
        dept_data->init(log_data->get_col_vector(dept_idx));

        // prepare string representation of properties from parent gis
        // TODO: find better way to copy props from one object to another
        const std::string log_prop_dump = g->get_prop()->to_str();

        // loop over all well log columns
        spv_double cur_log = BS_KERNEL.create_object(v_double::bs_type());
        sp_gis_t cur_gis = BS_KERNEL.create_object("gis");
        for(ulong i = 0; i < names.size(); ++i) {
          // skip depth
          if(i == dept_idx) continue;
          // create and fill new gis object
          cur_gis->get_prop()->from_str(log_prop_dump);
          sp_table_t cur_data = cur_gis->get_table();
          cur_data->init(0, 2);
          cur_data->add_col_vector(0, names[dept_idx], dept_data);
          cur_log->init(log_data->get_col_vector(i));
          cur_data->add_col_vector(1, names[i], cur_log);

          // write it to DB
          add_branch_gis(
            wname, branch, cur_gis, wstr2str(names[i], "utf-8"), wlog_type, replace_existing
          );
        }

        // 2. Write gis in old format if no wlog_name specified
        // TODO: deprecate and disable writing to branches table at all
        q = "UPDATE ";
        //if(!replace_existing)
        //  q += "OR IGNORE ";
        q += "branches SET well_log = ?1 WHERE well_name = '";
        q += wname + "' AND branch_name = '" + branch + "'";

        if(prepare_sql(q) < 0)
          return -1;
        // exec query
        return dump_wlog_data::go(*this, g);
      }

      // new implementation writes to separate well logs table
      // format query
      q = "INSERT OR ";
      if(replace_existing)
        q += "REPLACE";
      else
        q += "IGNORE";
      q += " INTO well_logs (well_name, branch_name, wlog_name, wlog_type, wlog_data) VALUES ('";
      q += wname + "', '" + branch + "', '" + wlog_name + "', " +
        boost::lexical_cast< std::string >(wlog_type) + ", ?1)";

      // prepare query
      if(prepare_sql(q) < 0)
        return -1;
      // exec query
      return dump_wlog_data::go(*this, g);
    }

  std::vector< std::string > sql_well::get_wlog_names(
      const std::string &wname, const std::string &branch, uint wlog_type
  ) {
      std::vector< std::string > res;
      if (!db)
        return res;

      std::string q = "SELECT wlog_name FROM well_logs WHERE well_name = '" + wname +
        "' AND branch_name = '" + branch + "' AND wlog_type = " +
        boost::lexical_cast< std::string >(wlog_type);

      // exec sql
      if(prepare_sql(q) < 0) {
        finalize_sql();
        return res;
      }

      while(step_sql() == 0) {
        res.push_back(get_sql_str(0));
      }

      finalize_sql();
      return res;
  }

  sql_well::sp_gis_t
  sql_well::get_branch_gis (const std::string &wname, const std::string &branch,
      std::string wlog_name, uint wlog_type)
    {
      // helper function
      struct extract_str_blob {
        static std::string go(sql_well& sqw, const int blob_idx) {
          // extract blob
          std::string s;
          const int n = sqlite3_column_bytes (sqw.stmp_sql, blob_idx);
          const char *b = (const char *)sqlite3_column_blob (sqw.stmp_sql, blob_idx);
          s.assign(b, n);
          return s;
        }
      };

      // result
      sp_gis_t sp_gis;

      if (!db)
        return sp_gis;
      finalize_sql();

      std::string q;
      const std::string select_filter = " WHERE well_name = '" + wname +
        "' AND branch_name = '" + branch + "'";
      // what blobs can we extract from result?
      bool has_data = false;

      if(wlog_name.size()) {
        // check well_logs table
        // format query
        q = "SELECT wlog_data FROM well_logs" + select_filter +
          " AND wlog_name = '" + wlog_name + "' AND wlog_type = " +
          boost::lexical_cast< std::string >(wlog_type);
        // exec sql
        if(prepare_sql(q) == 0 && step_sql() == 0) {
          has_data = true;
        }
      }
      else {
        // make old-fashioned query to branches table
        // format query
        q = "SELECT well_log FROM branches" + select_filter;
        // exec sql
        if(prepare_sql(q) == 0 && step_sql() == 0) {
          has_data = true;
        }
      }

      // 3. read traj BLOB
      if(has_data) {
        // leave this for debugging purposes
        std::cout << "READ WELL LOG: " << wname;

        // extract well log data
        q = extract_str_blob::go(*this, 0);
        if(!q.empty()) {
          // we have some to read
          sp_gis = BS_KERNEL.create_object ("gis");
          sp_gis->from_str(q);
        }
        std::cout << ", DATA = " << q.size() << std::endl;
      }
      finalize_sql();

      return sp_gis;
    }

   bool sql_well::rename_well_log(
       const std::string &wname, const std::string &branch,
       const std::string& old_name, const std::string& new_name
   ) {
      if (!db)
        return -1;
      if (stmp_sql)
        finalize_sql ();

      const std::string q = "UPDATE well_logs SET wlog_name = '" + new_name +
      "' WHERE well_name = '" + wname + "' AND branch_name = '" + branch +
      "' AND wlog_name = '" + old_name + "'";

      bool res = (exec_sql(q) == 0);

      // try to rename also in old-fashioned logs table inside branches
      sp_gis_t g = get_branch_gis(wname, branch);
      if(!g) return res;

      sp_table_t wlog_data = g->get_table();
      const std::wstring wold_name = str2wstr(old_name, "utf-8");
      const ulong n_wlogs(wlog_data->get_n_cols());
      for(ulong i = 0; i < n_wlogs; ++i) {
        if(wlog_data->get_col_name(i) != wold_name) continue;
        // we found a matching column
        // rename it
        wlog_data->set_col_name(i, str2wstr(new_name, "utf-8"));
        // and write result to DB
        add_branch_gis(wname, branch, g);
        return true;
      }

      return res;
   }

   bool sql_well::delete_well_log(
       const std::string &wname, const std::string &branch, std::string wlog_name
   ) {
      if (!db)
        return -1;
      if (stmp_sql)
        finalize_sql ();

      std::string q;
      const std::string select_filter = " WHERE well_name = '" + wname +
        "' AND branch_name = '" + branch + "'";
      bool res = false;

      if(wlog_name.size() == 0) {
        // old-fashioned query
        q = "UPDATE branches SET well_log = NULL" + select_filter;
        res = (exec_sql(q) == 0);
      }
      else {
        // delete from well_logs table
        q = "DELETE FROM well_logs" + select_filter +
          " AND wlog_name = '" + wlog_name + "'";
        res = (exec_sql(q) == 0);
      }

      // try to delete also from old-fashioned logs table inside branches
      sp_gis_t g = get_branch_gis(wname, branch);
      if(!g) return res;

      sp_table_t wlog_data = g->get_table();
      // delete corresponding column
      const std::wstring w_name = str2wstr(wlog_name, "utf-8");
      const ulong n_wlogs(wlog_data->get_n_cols());
      for(ulong i = 0; i < n_wlogs; ++i) {
        if(wlog_data->get_col_name(i) != w_name) continue;
        // we found a matching column
        // delete it
        wlog_data->remove_col(i);
        // and write result to DB
        add_branch_gis(wname, branch, g);
        return true;
      }

      return res;
   }

   int
   sql_well::add_branch_traj (const std::string &wname, const std::string &branch,
                              sp_traj_t t)
     {
      if (!db)
        return -1;
      if (stmp_sql)
        finalize_sql ();

      const std::string q = "UPDATE branches SET traj = ?1 WHERE well_name = '" +
      wname + "' AND branch_name = '" + branch + "'";

      // prepare query
      if(prepare_sql(q) < 0)
        return -1;

      // bind well log data
      const std::string traj_data = t->to_str();
      if (sqlite3_bind_blob (stmp_sql, 1, &traj_data.c_str()[0], traj_data.size (), SQLITE_STATIC))
        {
          fprintf (stderr, "Can't make select: %s\n", sqlite3_errmsg (db));
          finalize_sql();
          return -3;
        }

      // exec query
      step_sql();
      finalize_sql();

      return 0;
     }

   sql_well::sp_traj_t
   sql_well::get_branch_traj (const std::string &wname, const std::string &branch)
     {
      // helper function
      struct extract_str_blob {
        static std::string go(sql_well& sqw, const int blob_idx) {
          // extract blob
          std::string s;
          const int n = sqlite3_column_bytes (sqw.stmp_sql, blob_idx);
          const char *b = (const char *)sqlite3_column_blob (sqw.stmp_sql, blob_idx);
          s.assign(b, n);
          return s;
        }
      };

      sp_traj_t sp_traj;
      if (!db)
        return sp_traj;
      finalize_sql();

      // format sql
      std::string q = "SELECT traj FROM branches WHERE well_name = '" + wname +
        "' AND branch_name = '" + branch + "'";

      // exec sql
      if(prepare_sql(q) == 0 && step_sql() == 0) {
        // leave this for debugging purposes
        std::cout << "READ WELL TRAJ: " << wname;
        // extract trajectory
        q = extract_str_blob::go(*this, 0);
        if(!q.empty()) {
          sp_traj = BS_KERNEL.create_object ("traj");
          sp_traj->from_str(q);
        }
        std::cout << ", DATA = " << q.size() << std::endl;
      }
      finalize_sql();

      return sp_traj;
   }

#ifdef BSPY_EXPORTING_PLUGIN
  std::string
  sql_well::py_str () const
    {
      std::stringstream s;
      s << file_name << "\n";
      return s.str ();
    }
#endif //BSPY_EXPORTING_PLUGIN

/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (sql_well);
  BLUE_SKY_TYPE_STD_COPY (sql_well);

  BLUE_SKY_TYPE_IMPL(sql_well,  well_pool_iface, "sql_well", "sql_well storage", "realization of well sql_well storage");

}  // blue_sky namespace

