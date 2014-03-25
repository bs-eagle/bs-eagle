/// @file sql_well_db.cpp
/// @brief sql_well DB-related functionality implementation
/// @author uentity
/// @version 1.0
/// @date 25.03.2014
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifdef BSPY_EXPORTING_PLUGIN

#include <boost/python.hpp>
using namespace boost::python;

#endif //BSPY_EXPORTING_PLUGIN

#include "sql_well.h"
#include <stdio.h>
#include <fstream>
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/fstream.hpp"

namespace blue_sky {

// hidden details
namespace {

bool wlogs_table_exists(sql_well& sqw) {
  std::string q = "SELECT count(*) FROM sqlite_master WHERE type='table' AND name='well_logs'";
  sqw.prepare_sql(q);

  bool res = false;
  if(sqw.step_sql() == 0) {
    res = bool(sqw.get_sql_int(0));
  }
  sqw.finalize_sql();
  return res;
}

// retrive columns names of given table
std::set< std::string > table_columns(sql_well& sqw, const std::string& tbl_name) {
    // retrive wlog table columns
    std::set< std::string > col_names;
    const std::string q = "PRAGMA table_info(" + tbl_name + ")";
    //std::cout << tbl_name << " columns:" << std::endl;
    if(sqw.prepare_sql(q) == 0) {
      while(sqw.step_sql() == 0) {
        col_names.insert(sqw.get_sql_str(1));
        //std::cout << sqw.get_sql_str(1) << std::endl;
      }
    }
    return col_names;
}

// create well_log table for storing multiple well logs with possibly different depths
// i.e. we can store mmultiple gis objects for single well
bool create_wlogs_table(sql_well& sqw) {
  std::string q;
  bool res = true;
  //std::cout << "*** create_wlogs_table: " << sqw.file_name << std::endl;
  if(!wlogs_table_exists(sqw)) {
    //std::cout << "create_wlogs_table: no well_logs table in " << sqw.file_name << std::endl;
    // create table from scratch
    q = "CREATE TABLE IF NOT EXISTS \
      well_logs(\
      well_name TEXT NOT NULL,\
      branch_name TEXT NOT NULL,\
      wlog_name TEXT NOT NULL,\
      wlog_data BLOB,\
      wlog_type INT DEFAULT 0,\
      PRIMARY KEY(well_name, branch_name, wlog_name));\
      CREATE INDEX IF NOT EXISTS iwlname ON well_logs (wlog_name ASC);\
      CREATE UNIQUE INDEX IF NOT EXISTS iwlpkey ON well_logs (well_name, branch_name, wlog_name ASC);\
      ";
    res &= (sqw.exec_sql(q) != 0);
    //std::cout << "create_wlogs_table: well_logs table created = " << res << std::endl;
  }
  else {
    const std::set< std::string >& col_names = table_columns(sqw, "well_logs");
    if(!col_names.size())
      return false;

    // ensure that we always have wlog_type column
    // for old DBs
    //if(col_names.find("wlog_prop") == col_names.end()) {
    //  q = "ALTER TABLE well_logs ADD COLUMN wlog_prop BLOB";
    //  sqw.exec_sql(q);
    //}
    if(col_names.find("wlog_type") == col_names.end()) {
      q = "ALTER TABLE well_logs ADD COLUMN wlog_type INT DEFAULT 0";
      sqw.exec_sql(q);
    }
  }

  // fix wells table
  const std::set< std::string >& wells_cols = table_columns(sqw, "wells");
  if(!wells_cols.size())
    return false;
  if(wells_cols.find("KB") == wells_cols.end()) {
    q = "ALTER TABLE wells ADD COLUMN KB REAL DEFAULT -1";
    sqw.exec_sql(q);
  }
  //if(wells_cols.find("uid") == wells_cols.end()) {
  //  q = "ALTER TABLE wells ADD COLUMN uid TEXT DEFAULT ''";
  //  sqw.exec_sql(q);
  //}

  // fix branches table
  //const std::set< std::string >& br_cols = table_columns(sqw, "branches");
  //if(!br_cols.size())
  //  return false;
  //if(br_cols.find("traj_prop") == br_cols.end()) {
  //  q = "ALTER TABLE branches ADD COLUMN traj_prop BLOB";
  //  sqw.exec_sql(q);
  //}

  // convert well logs from old format to new
  q = "SELECT well_name, branch_name FROM branches";
  if(sqw.prepare_sql(q) < 0)
    return true;

  while(sqw.step_sql() == 0) {
    const std::string wname = sqw.get_sql_str(0);
    const std::string branch_name = sqw.get_sql_str(1);
    sql_well::sp_gis_t old_g = sqw.get_branch_gis(wname, branch_name);
    if(old_g)
      // don't overwrite existing logs
      sqw.add_branch_gis(wname, branch_name, old_g, "", 0, false);
  }
  sqw.finalize_sql();

  return true;
}

} // eof hidden namespace

  int
  sql_well::open_db (const std::string &file)
    {
      int rc = 0;
      char *zErrMsg = 0;

      if (db)
        close_db ();
      //const std::string file = wstr2str (file_);
      printf ("SQL open_db %s\n", file.c_str ());
      if (file == ":memory:" || !boost::filesystem::exists (file))
        {
          rc = sqlite3_open (file.c_str (), &db);
          if (rc)
            {
              fprintf (stderr, "Can't open database: %s (%s)\n", sqlite3_errmsg (db), file.c_str());
              sqlite3_close (db);
              db = 0;
              return -1;
            }
          create_db (db);

          rc = sqlite3_exec (db, "INSERT INTO groups(name) VALUES ('field');", NULL, NULL, &zErrMsg);
          if (rc != SQLITE_OK)
            {
              fprintf (stderr, "SQL error: %s\n", zErrMsg);
              sqlite3_free (zErrMsg);
            }
          //sqlite3_close (db);
        }
      else
        {
          rc = sqlite3_open (file.c_str (), &db);
          if (rc != SQLITE_OK)
            {
              fprintf (stderr, "Can't open database: %s (%s) - 2\n", sqlite3_errmsg (db), file.c_str());
              sqlite3_close (db);
              db = 0;
              return -1;
            }
        }

      file_name = file;
      return create_wlogs_table(*this) ? 0 : -1;
    }

  void
  sql_well::backup_to_file (const std::string &filename)
    {
      if (!db)
        return;
      sqlite3 *ddb = 0;
      int rc = sqlite3_open (filename.c_str(), &ddb);
      if (rc)
        {
          fprintf (stderr, "Can't open database: %s\n", sqlite3_errmsg (ddb));
          sqlite3_close (ddb);
          ddb = 0;
          return;
        }
      sqlite3_backup *b = sqlite3_backup_init(ddb, "main", db, "main");
      if (b)
        {
          while(sqlite3_backup_step(b, -1) == SQLITE_OK) {}
	        //while (sqlite3_backup_remaining(b) > 0)
          sqlite3_backup_finish(b);
        }

      rc = sqlite3_errcode(ddb);
      if( rc != SQLITE_OK )
      {
        fprintf (stderr, "SQL error with backup_to_file: %d\n", rc);
      }
      sqlite3_close(ddb);
      //db->sqlite_backup_to_file(filename);
    }


  int
  sql_well::merge_with_db(const std::string& dbname)
    {
      if (!db)
        return 0;
      if (stmp_sql)
        finalize_sql ();
      //std::string dbname = wstr2str (dbname_);

      int rc = 0;
      char *zErrMsg = 0;

      std::string sql = "attach '" + dbname + "' as tomerge;    insert or ignore into groups select * from tomerge.groups; \
                                                                insert or ignore into dates select * from tomerge.dates; \
                                                                insert or ignore into wells select * from tomerge.wells; \
                                                                insert or replace into branches select * from tomerge.branches; \
                                                                insert or ignore into completions select * from tomerge.completions; \
                                                                insert or ignore into fractures select * from tomerge.fractures; \
                                                                insert or ignore into permfrac select * from tomerge.permfrac; \
                                                                insert or ignore into well_hist select * from tomerge.well_hist; \
                                                                insert or ignore into wells_in_group select * from tomerge.wells_in_group; \
                                                                detach database tomerge";

      rc = sqlite3_exec (db, sql.c_str (), NULL, 0, &zErrMsg);
      if( rc != SQLITE_OK )
      {
        fprintf (stderr, "SQL error with tomerge: %s\n", zErrMsg);
        sqlite3_free (zErrMsg);
        return -4;
      }
      return 0;
    }

  int
  sql_well::create_db (sqlite3 *db_in)
    {
      int rc = 0;
      char *zErrMsg = 0;
      const char *sql = \
"\
PRAGMA foreign_keys=ON;\
BEGIN;\
CREATE TABLE wells(name TEXT UNIQUE PRIMARY KEY, \
				    x REAL DEFAULT -1, \
				    y REAL DEFAULT -1, \
				    horiz INTEGER DEFAULT 0, \
				    KB REAL DEFAULT -1, \
				    src INT DEFAULT 0);\
CREATE TABLE groups(name TEXT UNIQUE PRIMARY KEY);\
CREATE TABLE wells_in_group(gr_name TEXT NOT NULL REFERENCES groups(name) ON UPDATE CASCADE ON DELETE CASCADE,\
						    well_name TEXT NOT NULL REFERENCES wells(name) ON UPDATE CASCADE ON DELETE CASCADE);\
CREATE INDEX i4 ON wells_in_group (gr_name ASC);\
CREATE INDEX i5 ON wells_in_group (well_name ASC);\
CREATE TABLE dates(d REAL UNIQUE PRIMARY KEY);\
CREATE TABLE well_hist(well_name TEXT NOT NULL REFERENCES wells(name) ON UPDATE CASCADE ON DELETE CASCADE,\
					 d REAL NOT NULL, \
					 p_or REAL DEFAULT -1,\
					 p_wr REAL DEFAULT -1,\
					 p_gr REAL DEFAULT -1,\
					 p_lr REAL DEFAULT -1,\
					 p_bhp REAL DEFAULT -1,\
                     p_fgr REAL DEFAULT -1,\
					 i_or REAL DEFAULT -1,\
					 i_wr REAL DEFAULT -1,\
					 i_gr REAL DEFAULT -1,\
					 i_bhp REAL DEFAULT -1,\
                     wefac REAL DEFAULT 1.0,\
					 ctrl INTEGER DEFAULT 0,\
					 status INTEGER DEFAULT 0,\
					 lim_p_or REAL DEFAULT -1,\
					 lim_p_wr REAL DEFAULT -1,\
					 lim_p_gr REAL DEFAULT -1,\
                     lim_p_lr REAL DEFAULT -1,\
					 lim_p_bhp REAL DEFAULT -1,\
					 lim_i_or REAL DEFAULT -1,\
					 lim_i_wr REAL DEFAULT -1,\
					 lim_i_gr REAL DEFAULT -1,\
					 lim_i_bhp REAL DEFAULT -1);\
CREATE INDEX i1 ON well_hist (well_name ASC);\
CREATE INDEX i2 ON well_hist (d ASC);\
CREATE UNIQUE INDEX i3 ON well_hist (well_name, d ASC);\
CREATE TABLE well_res(well_name TEXT NOT NULL REFERENCES wells(name) ON UPDATE CASCADE ON DELETE CASCADE,\
					 d REAL NOT NULL, \
					 p_or REAL DEFAULT -1,\
					 p_wr REAL DEFAULT -1,\
					 p_gr REAL DEFAULT -1,\
					 p_bhp REAL DEFAULT -1,\
					 i_or REAL DEFAULT -1,\
					 i_wr REAL DEFAULT -1,\
					 i_gr REAL DEFAULT -1,\
					 i_bhp REAL DEFAULT -1,\
					 wefac REAL DEFAULT 1,\
					 p_gor REAL DEFAULT -1,\
					 p_fgr REAL DEFAULT -1,\
					 ctrl INTEGER DEFAULT 0,\
					 status INTEGER DEFAULT 0,\
					 tot_p_or REAL DEFAULT -1,\
					 tot_p_wr REAL DEFAULT -1,\
					 tot_p_gr REAL DEFAULT -1,\
					 tot_p_fgr REAL DEFAULT -1,\
					 tot_i_or REAL DEFAULT -1,\
					 tot_i_wr REAL DEFAULT -1,\
					 tot_i_gr REAL DEFAULT -1);\
CREATE INDEX i6 ON well_res (well_name ASC);\
CREATE INDEX i7 ON well_res (d ASC);\
CREATE UNIQUE INDEX i8 ON well_res (well_name, d ASC);\
CREATE TABLE branches(well_name TEXT NOT NULL REFERENCES wells(name) ON UPDATE CASCADE ON DELETE CASCADE,\
					   branch_name TEXT NOT NULL DEFAULT 'main', \
                       md REAL DEFAULT -1,\
                       parent TEXT DEFAULT '',\
					   traj BLOB, \
					   well_log BLOB, \
					   PRIMARY KEY (well_name, branch_name));\
CREATE INDEX i9 ON branches (well_name ASC);\
CREATE UNIQUE INDEX i10 ON branches (well_name, branch_name ASC);\
CREATE TRIGGER tr1 AFTER INSERT ON wells\
	BEGIN\
		INSERT INTO branches(well_name, branch_name) VALUES(new.name, 'main');\
	END;\
CREATE TRIGGER tr2 AFTER INSERT ON wells\
	BEGIN\
		INSERT INTO wells_in_group(gr_name, well_name) VALUES ('field', new.name);\
	END;\
CREATE TABLE fractures(well_name TEXT NOT NULL,\
					     branch_name TEXT DEFAULT 'main', \
					     md REAL NOT NULL, \
					     d REAL NOT NULL, \
					     status  INTEGER DEFAULT 0,\
					     half_up REAL DEFAULT 5,\
					     half_down REAL DEFAULT 5,\
					     angle REAL DEFAULT 0,\
					     half_length_1 REAL DEFAULT 50,\
					     half_length_2 REAL DEFAULT 50,\
					     perm REAL DEFAULT -1,\
					     half_thin REAL DEFAULT 0.005,\
               skin REAL DEFAULT 0,\
					     FOREIGN KEY (well_name, branch_name) REFERENCES branches(well_name, branch_name) ON UPDATE CASCADE ON DELETE CASCADE\
					     );\
CREATE INDEX i11 ON fractures (well_name ASC);\
CREATE INDEX i12 ON fractures (well_name, branch_name ASC);\
CREATE TABLE completions(well_name TEXT NOT NULL, \
					     branch_name TEXT NOT NULL DEFAULT 'main', \
					     md REAL NOT NULL, \
					     d REAL NOT NULL, \
					     status  INTEGER DEFAULT 0,\
					     length REAL DEFAULT 1,\
					     rw REAL DEFAULT 0.08,\
					     r0 REAL DEFAULT -1,\
					     kh REAL DEFAULT -1,\
					     kh_mult REAL DEFAULT 1,\
               skin REAL DEFAULT 0,\
					     FOREIGN KEY (well_name, branch_name) REFERENCES branches(well_name, branch_name) ON UPDATE CASCADE ON DELETE CASCADE\
					     );					     \
CREATE INDEX i13 ON completions (well_name ASC);\
CREATE INDEX i14 ON completions (well_name, branch_name ASC);\
CREATE TRIGGER tr3 BEFORE INSERT ON fractures\
	BEGIN\
		INSERT OR REPLACE INTO dates(d) VALUES(new.d);\
	END;\
CREATE TRIGGER tr4 BEFORE INSERT ON completions\
	BEGIN\
		INSERT OR REPLACE INTO dates(d) VALUES(new.d);\
	END;\
CREATE TRIGGER tr5 BEFORE INSERT ON well_hist\
	BEGIN\
		INSERT OR REPLACE INTO dates(d) VALUES(new.d);\
	END;\
CREATE TABLE permfrac(d REAL NOT NULL, \
					     status  INTEGER DEFAULT 0,\
					     x1 REAL NOT NULL,\
					     y1 REAL NOT NULL,\
					     x2 REAL NOT NULL,\
					     y2 REAL NOT NULL,\
					     z_top REAL NOT NULL,\
					     z_bottom REAL NOT NULL,\
					     dx REAL DEFAULT 0.01,\
					     dy REAL DEFAULT 0.01,\
					     perm REAL NOT NULL,\
					     i1  INTEGER DEFAULT -1,\
					     j1  INTEGER DEFAULT -1,\
					     i2  INTEGER DEFAULT -1,\
					     j2  INTEGER DEFAULT -1,\
					     k1  INTEGER DEFAULT -1,\
					     k2  INTEGER DEFAULT -1);\
COMMIT;\
";
      if (!db_in)
        return -1;
      rc = sqlite3_exec (db_in, sql, NULL, 0, &zErrMsg);
      if (rc != SQLITE_OK)
        {
          fprintf(stderr, "SQL error: %s\n", zErrMsg);
          sqlite3_free (zErrMsg);
        }
      return 0;
    }


  void
  sql_well::close_db ()
    {
      if (db)
        {


          if (stmp_sql)
            finalize_sql ();

#if 0
          int rc = 0;
          char *zErrMsg = 0;
          char buf[2048];
          sprintf (buf, "ATTACH DATABASE '%s' as backup; BEGIN", file_name.c_str ());
          rc = sqlite3_exec (db, buf, NULL, NULL, &zErrMsg);
          if (rc != SQLITE_OK)
            {
              fprintf (stderr, "SQL error: %s\n", zErrMsg);
              sqlite3_free (zErrMsg);
            }
          rc = sqlite3_exec (db, "SELECT name FROM backup.sqlite_master WHERE type='table'",
                            &clear_table, db, &zErrMsg);
          if (rc != SQLITE_OK)
            {
              fprintf (stderr, "SQL error: %s\n", zErrMsg);
              sqlite3_free (zErrMsg);
            }

          rc = sqlite3_exec (db, "SELECT name FROM main.sqlite_master WHERE type='table'",
                             &main_to_backup, db, &zErrMsg);
          if (rc != SQLITE_OK)
            {
              fprintf (stderr, "SQL error: %s\n", zErrMsg);
              sqlite3_free (zErrMsg);
            }
          sqlite3_exec(db, "COMMIT; DETACH DATABASE backup", NULL, NULL, NULL);
#endif //0
          sqlite3_close (db);
        }

      db = 0;
      file_name = "";
    }


  int
  sql_well::create_db_struct ()
    {
      if(create_db (db)) {
        return create_wlogs_table(*this) ? 0 : -1;
      }
      return -1;
    }

  void
  sql_well::fill_db ()
    {
      //int rc = 0;
      //char *zErrMsg = 0;
      char sw_name[4096];
      char sw_sel[4096];
      char sw_ins[4096];
      char sw_up[4096];
      if (!db)
        return;

      if (stmp_sql)
        finalize_sql ();
      //rc = sqlite3_exec (db, "BEGIN TRANSACTION", NULL, 0, &zErrMsg);
      for (int i = 0; i < 500; ++i)
        {
          sprintf (sw_name, "well_%d", i);
          add_well (std::string (sw_name));
          for (double j = 0; j < 240; j = j + 1.0)
            {
              sprintf (sw_ins, "INSERT INTO well_dynamic (well_name, date, h_wrate) VALUES ('well_%d', %lf, %lf)",
                       i, j, 50.0);
              sprintf (sw_up, "UPDATE well_dynamic SET h_orate = %lf WHERE well_name='well_%d' and date=%lf",
                       j, i, j);
              sprintf (sw_sel, "SELECT * FROM well_dynamic WHERE well_name='well_%d' and date=%lf",
                       i, j);
              //printf ("%s\n", sw);
              if (insert_or_update (std::string (sw_sel),
                                    std::string (sw_ins),
                                    std::string (sw_up)))
                {
                  return;
                }
            }
        }
      //rc = sqlite3_exec (db, "COMMIT TRANSACTION", NULL, 0, &zErrMsg);
      return;
    }
  int
  sql_well::insert_or_update (const std::string &select_sql,
                              const std::string &insert_sql,
                              const std::string &update_sql)
    {
      if (!db)
        return -2;

      if (stmp_sql)
        finalize_sql ();
      int rc = 0;
      char *zErrMsg = 0;
      const char *ttt;
      sqlite3_stmt *stmp;

      rc = sqlite3_prepare_v2 (db, select_sql.c_str (), select_sql.length () + 1, &stmp, &ttt);
      if (rc)
        {
          fprintf (stderr, "Can't make select: %s\n", sqlite3_errmsg (db));
          return -1;
        }
      if (sqlite3_step (stmp) == SQLITE_ROW) // UPDATE
        {
          sqlite3_finalize (stmp);

          rc = sqlite3_exec (db, update_sql.c_str (), NULL, 0, &zErrMsg);
          if( rc != SQLITE_OK )
            {
              fprintf (stderr, "SQL error: %s\n", zErrMsg);
              sqlite3_free (zErrMsg);
              return -4;
            }
        }
      else
        {
          sqlite3_finalize (stmp);
          rc = sqlite3_exec (db, insert_sql.c_str (), NULL, 0, &zErrMsg);
          if( rc != SQLITE_OK )
            {
              fprintf (stderr, "SQL error: %s\n", zErrMsg);
              sqlite3_free (zErrMsg);
              return -4;
            }
        }

      return 0;
    }

  int
  sql_well::prepare_sql (const std::string &sql)
    {
      if (!db)
        return -1;
      if (stmp_sql)
        finalize_sql ();
      const char *ttt;
      int rc = 0;
      
      rc = sqlite3_prepare_v2 (db, sql.c_str (), sql.length () + 1, &stmp_sql, &ttt);
      if (rc)
        {
          fprintf (stderr, "Can't make select: %s (%d)\n", sqlite3_errmsg (db), rc);
          fprintf (stderr, "Query: %s\n", sql.c_str());
          finalize_sql ();
          return -1;
        }
      return 0;
    }

//#ifdef BSPY_EXPORTING_PLUGIN
  spv_double
  sql_well::get_table (const std::string &table_name, boost::python::list &table_columns, const std::string &filter)
    {
      std::string request_str = "SELECT COUNT(*) FROM " + table_name + " " + filter;
      prepare_sql (request_str);
      step_sql ();
      int n_rows = get_sql_int(0);
      finalize_sql ();

      if (n_rows < 1)
        return BS_KERNEL.create_object (v_double::bs_type());

      int n_cols = boost::python::len (table_columns);
      spv_double table = BS_KERNEL.create_object (v_double::bs_type());
      table->resize (n_rows * n_cols);
      double *table_values = &(*table)[0];

      request_str = "SELECT ";
      for (int j = 0; j < n_cols; j++)
        {
          std::string col_name = extract<std::string>(table_columns[j]);
          //TODO: check col_name
          request_str = request_str + col_name;
          if (j < n_cols - 1)
            request_str = request_str + ", ";
        }
      request_str = request_str + " FROM " + table_name + " " + filter;
      prepare_sql (request_str);


      for (int i = 0; (i < n_rows) && (!step_sql ()); ++i)
        {
          for (int j = 0; j < n_cols; j++)
            {
              table_values[i * n_cols + j] = get_sql_real(j);
            }
        }
      finalize_sql ();

      npy_intp dims[] = {n_rows, n_cols};
      table->reshape (2, dims);
      return table;
    }
//#endif //BSPY_EXPORTING_PLUGIN

  int
  sql_well::step_sql ()
    {
      if (!db)
        return -1;
      if (!stmp_sql)
        return -1;
      if (sqlite3_step (stmp_sql) == SQLITE_ROW) // UPDATE
        return 0;
      else
        return 2;
    }
  int
  sql_well::finalize_sql ()
    {
      if (stmp_sql)
        {
          sqlite3_finalize (stmp_sql);
          stmp_sql = 0;
        }
      return 0;
    }
  t_int
  sql_well::get_sql_int (t_int col)
    {
      if (!db)
        return 0;
      if (!stmp_sql)
        return 0;
      return (t_int)sqlite3_column_int (stmp_sql, (int)col);
    }
  t_double
  sql_well::get_sql_real (t_int col)
    {
      if (!db)
        return 0;
      if (!stmp_sql)
        return 0;
      return (t_double)sqlite3_column_double (stmp_sql, (int)col);
    }
  bool
  sql_well::get_sql_bool (t_int col)
    {
      if (!db)
        return 0;
      if (!stmp_sql)
        return 0;
      return (bool)sqlite3_column_int (stmp_sql, (int)col);
    }
  std::string
  sql_well::get_sql_str (t_int col)
    {
      if (!db)
        return 0;
      if (!stmp_sql)
        return 0;
      return std::string ((const char *)sqlite3_column_text (stmp_sql, (int)col));
    }
  bool
  sql_well::get_sql_exist (t_int col)
  {
	  if (!db)
		  return false;
	  if (!stmp_sql)
		  return false;
	  int n = sqlite3_column_bytes(stmp_sql, (int)col);
	  return (bool)n;
  }

  int
  sql_well::exec_sql (const std::string &sql)
    {
      if (!db)
        return 0;
      if (stmp_sql)
        finalize_sql ();

      int rc = 0;
      char *zErrMsg = 0;

      rc = sqlite3_exec (db, sql.c_str (), NULL, 0, &zErrMsg);
      if( rc != SQLITE_OK )
        {
          fprintf (stderr, "SQL error: %s, when query \"%s\"\n", zErrMsg, sql.c_str());
          sqlite3_free (zErrMsg);
          return -4;
        }
      return 0;
    }

  int
  sql_well::exec_sql_and_return_rowid (const std::string &sql)
    {
      if (!db)
        return 0;
      if (stmp_sql)
        finalize_sql ();

      int rc = 0;
      char *zErrMsg = 0;
      int rowid = -1;

      rc = sqlite3_exec (db, sql.c_str (), NULL, 0, &zErrMsg);
      rowid = sqlite3_last_insert_rowid(db);
      if( rc != SQLITE_OK )
        {
          fprintf (stderr, "SQL error: %s\n", zErrMsg);
          sqlite3_free (zErrMsg);
          return -4;
        }
      return rowid;
    }

#if 0
  int clear_table (void *pData, int nColumns,
                   char **values, char ** /*columns*/)
    {
      int rc = 0;
      char *zErrMsg = 0;
      if (nColumns != 1)
        return 1; // Error

      sqlite3* db = (sqlite3*)pData;

      char *stmt = sqlite3_mprintf("DELETE FROM backup.%q",
                                   values[0]);
      rc = sqlite3_exec (db, stmt, NULL, NULL, &zErrMsg);
      sqlite3_free (stmt);
      if (rc != SQLITE_OK)
        {
          fprintf (stderr, "SQL error: %s\n", zErrMsg);
          sqlite3_free (zErrMsg);
        }

      return 0;
    }
  int main_to_backup (void *pData, int nColumns,
                      char **values, char ** /*columns*/)
    {
      int rc = 0;
      char *zErrMsg = 0;
      if (nColumns != 1)
        return 1; // Error

      sqlite3* db = (sqlite3*)pData;

      char *stmt = sqlite3_mprintf("insert into backup.%q select * from main.%q",
                                   values[0], values[0]);
      rc = sqlite3_exec (db, stmt, NULL, NULL, &zErrMsg);
      sqlite3_free (stmt);
      if (rc != SQLITE_OK)
        {
          fprintf (stderr, "SQL error: %s\n", zErrMsg);
          sqlite3_free (zErrMsg);
        }

      return 0;
    }

  int backup_to_main (void *pData, int nColumns,
                      char **values, char ** /*columns*/)
    {
      int rc = 0;
      char *zErrMsg = 0;
      if (nColumns != 1)
        return 1; // Error

      sqlite3* db = (sqlite3*)pData;

      char *stmt = sqlite3_mprintf("insert into main.%q select * from backup.%q",
                                   values[0], values[0]);
      rc = sqlite3_exec (db, stmt, NULL, NULL, &zErrMsg);
      sqlite3_free (stmt);
      if (rc != SQLITE_OK)
        {
          fprintf (stderr, "SQL error: %s\n", zErrMsg);
          sqlite3_free (zErrMsg);
        }

      return 0;
    }
#endif

} /* namespace blue_sky */

