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

#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <string>
#include <fstream>
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/fstream.hpp"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>



#include "bs_kernel.h"
#include "sql_well.h"
#include "frac_comp_ident.h"
#include "well_path_ident.h"
#include "i_cant_link_2_mesh.h"

using namespace boost;

#ifdef BSPY_EXPORTING_PLUGIN

#include <boost/python.hpp>
using namespace boost::python;

#endif //BSPY_EXPORTING_PLUGIN


namespace blue_sky
{

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
    }

  int
  sql_well::open_db (const std::string &file)
    {
      int rc = 0;
      char *zErrMsg = 0;

      if (db)
        close_db ();
      printf ("SQL open_db %s\n", file.c_str ());
      if (!strcmp(file.c_str(),":memory:") || !boost::filesystem::exists (file))
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
          if (rc)
            {
              fprintf (stderr, "Can't open database: %s (%s) - 2\n", sqlite3_errmsg (db), file.c_str());
              sqlite3_close (db);
              db = 0;
              return -1;
            }
        }


#if 0
      char buf[2048];
      rc = sqlite3_open (":memory:", &db);
      if (rc)
        {
          fprintf (stderr, "Can't open database: %s\n", sqlite3_errmsg (db));
          sqlite3_close (db);
          db = 0;
          return -1;
        }
      file_name = file;
      // load from file to memory
      rc = create_db (db);
      if (rc)
        return rc;
      sprintf (buf, "ATTACH DATABASE '%s' as backup; BEGIN", file.c_str ());
      rc = sqlite3_exec (db, buf, NULL, NULL, &zErrMsg);
      if (rc != SQLITE_OK)
        {
          fprintf (stderr, "SQL error: %s\n", zErrMsg);
          sqlite3_free (zErrMsg);
        }

      rc = sqlite3_exec(db, "SELECT name FROM backup.sqlite_master WHERE type='table'",
                        &backup_to_main, db, &zErrMsg);
      if (rc != SQLITE_OK)
        {
          fprintf (stderr, "SQL error: %s\n", zErrMsg);
          sqlite3_free (zErrMsg);
        }
      sqlite3_exec(db, "COMMIT; DETACH DATABASE backup", NULL, NULL, NULL);
#else //0
#endif //0
      return 0;
    }

  void
  sql_well::backup_to_file (const std::string &filename)
    {
      if (!db)
        return;
      sqlite3 *ddb = 0;
      int rc = sqlite3_open (filename.c_str (), &ddb);
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
          sqlite3_backup_step(b, -1);
          sqlite3_backup_finish(b);
        }
      rc = sqlite3_errcode(ddb);
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

      int rc = 0;
      char *zErrMsg = 0;

      std::string sql = "attach '" + dbname + "' as tomerge;    insert or ignore into groups select * from tomerge.groups; \
                                                                insert or ignore into wells select * from tomerge.wells; \
                                                                insert or ignore into completions select * from tomerge.completions; \
                                                                insert or ignore into dates select * from tomerge.dates; \
                                                                insert or ignore into fractures select * from tomerge.fractures; \
                                                                insert or ignore into permfrac select * from tomerge.permfrac; \
                                                                insert or ignore into well_hist select * from tomerge.well_hist; \
                                                                insert or ignore into wells_in_group select * from tomerge.wells_in_group; \
                                                                insert or replace into branches select * from tomerge.branches; \
                                                                detach database tomerge";

      rc = sqlite3_exec (db, sql.c_str (), NULL, 0, &zErrMsg);
      if( rc != SQLITE_OK )
      {
        fprintf (stderr, "SQL error: %s\n", zErrMsg);
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
				    horiz INTEGER DEFAULT 0);\
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
					   well_log BLOB,\
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
      return create_db (db);
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
                            sp_gis_t g)
    {
      if (!db)
        return -1;
      if (stmp_sql)
        finalize_sql ();

      int rc = 0;
      //char *zErrMsg = 0;
      const char *ttt;
      sqlite3_stmt *stmp;
      char buf[4096];
      sprintf (buf, "UPDATE branches SET well_log = @q WHERE well_name = '%s' AND branch_name = '%s'",
               wname.c_str (), branch.c_str ());
      rc = sqlite3_prepare_v2 (db, buf, strlen (buf), &stmp, &ttt);
      if (rc)
        {
          fprintf (stderr, "Can't make select: %s\n", sqlite3_errmsg (db));
          return -1;
        }

      std::ostringstream oss;
      boost::archive::text_oarchive oar(oss);
      g->save (oar);
      std::string s = oss.str ();
      std::vector<char> ch (s.begin (), s.end ());
      //printf ("GIS\n%s\n", s.c_str ());
      //printf ("GIS INT %d %d\n", (int)strlen (s.c_str ()), (int)s.length ());
      rc = sqlite3_bind_blob (stmp, 1, &ch[0], ch.size (), SQLITE_STATIC);
      if (rc)
        {
          fprintf (stderr, "Can't make select: %s\n", sqlite3_errmsg (db));
          return -3;
        }
      sqlite3_step (stmp); // UPDATE
      sqlite3_finalize (stmp);

      return 0;
    }

  sql_well::sp_gis_t
  sql_well::get_branch_gis (const std::string &wname, const std::string &branch) const
    {
      sp_gis_t sp_gis = BS_KERNEL.create_object ("gis");
      if (!db)
        return sp_gis_t ();
      if (stmp_sql)
        return sp_gis_t ();

      int rc = 0;
      //char *zErrMsg = 0;
      const char *ttt;
      sqlite3_stmt *stmp;
      char buf[4096];
      sprintf (buf, "SELECT well_log FROM branches WHERE well_name = '%s' AND branch_name = '%s'",
               wname.c_str (), branch.c_str ());
      rc = sqlite3_prepare_v2 (db, buf, strlen (buf) + 1, &stmp, &ttt);
      if (rc)
        {
          fprintf (stderr, "Can't make select: %s\n", sqlite3_errmsg (db));
          sqlite3_finalize (stmp);
          return sp_gis_t ();
        }
      if (sqlite3_step (stmp) == SQLITE_ROW) // UPDATE
        {
          int n = sqlite3_column_bytes (stmp, 0);
          if (n < 1)
            {
              sqlite3_finalize (stmp);
              return sp_gis_t ();
            }
          const char *b = (const char *)sqlite3_column_blob (stmp, 0);
          //std::string s = (const char *)sqlite3_column_text (stmp, 0);
          std::string s;
          s.assign (b, n);
          printf ("READ GIS %d\n", (int)s.length ());
          std::istringstream iss;
          iss.str (s);
          boost::archive::text_iarchive iar(iss);
          sp_gis->load (iar);

        }
	  else
	    {
		  return sp_gis_t ();
	    }
      sqlite3_finalize (stmp);
      return sp_gis;
    }

   int
   sql_well::add_branch_traj (const std::string &wname, const std::string &branch,
                              sp_traj_t t)
     {
      if (!db)
        return -1;
      if (stmp_sql)
        finalize_sql ();

      int rc = 0;
      //char *zErrMsg = 0;
      const char *ttt;
      sqlite3_stmt *stmp;
      char buf[4096];
      sprintf (buf, "UPDATE branches SET traj = @w WHERE well_name = '%s' AND branch_name = '%s'",
               wname.c_str (), branch.c_str ());
      rc = sqlite3_prepare_v2 (db, buf, strlen (buf) + 1, &stmp, &ttt);
      if (rc)
        {
          fprintf (stderr, "Can't make select: %s\n", sqlite3_errmsg (db));
          return -1;
        }

      std::ostringstream oss;
      boost::archive::text_oarchive oar(oss);
      t->save (oar);
      std::string s = oss.str ();
      std::vector<char> ch (s.begin (), s.end ());
      rc = sqlite3_bind_blob (stmp, 1, &ch[0], ch.size (), SQLITE_STATIC);
      //printf ("TRAJ\n%s\n", s.c_str ());
      //printf ("TRAJ INT %d %d\n", (int)strlen (s.c_str ()), (int)s.length ());
      if (rc)
        {
          fprintf (stderr, "Can't make select: %s\n", sqlite3_errmsg (db));
          return -3;
        }

      sqlite3_step (stmp); // UPDATE
      sqlite3_finalize (stmp);

      return 0;
     }

   sql_well::sp_traj_t
   sql_well::get_branch_traj (const std::string &wname, const std::string &branch) const
     {
      sp_traj_t sp_traj = BS_KERNEL.create_object ("traj");
      if (!db)
        return sp_traj_t ();
      if (stmp_sql)
        return sp_traj_t ();

      int rc = 0;
      //char *zErrMsg = 0;
      const char *ttt;
      sqlite3_stmt *stmp;
      char buf[4096];
      sprintf (buf, "SELECT traj FROM branches WHERE well_name = '%s' AND branch_name = '%s'",
               wname.c_str (), branch.c_str ());
      rc = sqlite3_prepare_v2 (db, buf, strlen (buf) + 1, &stmp, &ttt);
      if (rc)
        {
          fprintf (stderr, "Can't make select: %s\n", sqlite3_errmsg (db));
          sqlite3_finalize (stmp);
          return sp_traj_t ();
        }
      if (sqlite3_step (stmp) == SQLITE_ROW) // UPDATE
        {
          int n = sqlite3_column_bytes (stmp, 0);
          if (!n)
            {
              sqlite3_finalize (stmp);
              return sp_traj_t ();
            }
          const char *b = (const char *)sqlite3_column_blob (stmp, 0);
          //std::string s = (const char *)sqlite3_column_text (stmp, 0);
          std::string s;
          s.assign (b, n);
          printf ("READ TRAJ %d\n", (int)s.length ());
          std::istringstream iss;
          iss.str (s);
          //printf ("hkdjhkf: %s\n", iss.str ().c_str ());
          boost::archive::text_iarchive iar(iss);
          sp_traj->load (iar);

        }
	  else
	    {
		  return sp_traj_t ();
	    }
      sqlite3_finalize (stmp);
      return sp_traj;
     }

  int
   sql_well::update_branch_traj (const std::string &wname, const std::string &branch,
                                 sp_traj_t t)
     {
      if (!db)
        return -1;
      if (stmp_sql)
        finalize_sql ();

      int rc = 0;
      //char *zErrMsg = 0;
      const char *ttt;
      sqlite3_stmt *stmp;
      char buf[4096];
      sprintf (buf, "SELECT traj FROM branches WHERE well_name = '%s' AND branch_name = '%s'",
               wname.c_str (), branch.c_str ());
      rc = sqlite3_prepare_v2 (db, buf, strlen (buf) + 1, &stmp, &ttt);
      if (rc)
        {
          fprintf (stderr, "Can't make select: %s\n", sqlite3_errmsg (db));
          if(stmp) sqlite3_finalize (stmp);
          return -1;
        }

      std::ostringstream oss;
      boost::archive::text_oarchive oar(oss);
      t->save (oar);
      std::string s = oss.str ();
      std::vector<char> ch (s.begin (), s.end ());
      rc = sqlite3_bind_blob (stmp, 1, &ch[0], ch.size (), SQLITE_STATIC);
      //printf ("TRAJ\n%s\n", s.c_str ());
      //printf ("TRAJ INT %d %d\n", (int)strlen (s.c_str ()), (int)s.length ());
      if (rc)
        {
          fprintf (stderr, "Can't make select: %s\n", sqlite3_errmsg (db));
          if(stmp) sqlite3_finalize (stmp);
          return -3;
        }

      sqlite3_step (stmp); // UPDATE
      sqlite3_finalize (stmp);

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
          fprintf (stderr, "SQL error: %s\n", zErrMsg);
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

  int
  sql_well::read_from_ascii_file (const std::string &fname, double starting_date)
    {
      //if (fr_file)
      //  delete fr_file;
      //fr_file = new FRead (fname.c_str (), fname.c_str ());
      fr_file = BS_KERNEL.create_object ("bos_reader");
      fr_file->open (fname.c_str (), fname.c_str ());

      int rc = 0;
      char buf[4096];
      char *id = 0;
      char *other = 0;
      double d = 0;

      for (;;)
        {
          rc = fr_file->read_line (buf, 4096, FREAD_DONT_CONVERT_CASE);
          //printf ("%s\n", buf);
          if (rc <= 0)
            break;
          d = starting_date;
          // read date and time
          rc = read_date_and_time (buf, &id, &d);
          if (rc)
            return -4;
          // read id
          fr_file->trim_left (&id);
          other = id;
          rc |= fr_file->get_phrase (&other);
          if (rc)
            return -1;
          //trim_right_s (id);
          fr_file->locale_ucase (id);

          if (!strcmp (id, "W_SPEC"))
            {
              rc = read_w_spec (other);
              if (rc)
                return rc;
            }
          else if (!strcmp (id, "W_BRANCH_F"))
            {
              rc = read_w_branch_f (other);
              if (rc)
                return rc;
            }
          else if (!strcmp (id, "W_COMP"))
            {
              rc = read_w_comp (other, d);
              if (rc)
                return rc;
            }
          else if (!strcmp (id, "W_FRAC"))
            {
              rc = read_w_frac (other, d);
              if (rc)
                return rc;
            }
          else if (!strcmp (id, "W_PROD"))
            {
              rc = read_w_prod (other, d);
              if (rc)
                return rc;
            }
          else if (!strcmp (id, "W_INJ"))
            {
              rc = read_w_inj (other, d);
              if (rc)
                return rc;
            }
          //printf ("date %lf\n", d);
          //printf ("next: %s\n", next);
        }
      return 0;
    }

  int
  sql_well::read_w_branch_f (char *buf)
    {
      int rc = 0;
      char *nx = 0;
      char wname[1024];
      char branch[1024];
      char parent[1024];
      char fname[1024];
      char fname2[1024];

      char sql[1024];
      double md = -1;
      // read well name
      wname[0] = '\0';
      nx = buf;
      // read well name
      rc = fr_file->get_phrase_str (&nx, wname);
      if (rc || wname[0] == '\0')
        {
          fprintf (stderr, "Error: well name in keyword W_BRANCH_F can not be set by default\n");
          return -1;
        }
      // read branch
      strcpy (branch, "main");
      rc = fr_file->get_phrase_str (&nx, branch);
      if (rc)
        {
          fprintf (stderr, "Error: W_BRANCH_F can not read branch name\n");
          return -1;
        }
      // read parent
      parent[0] = '\0';
      rc = fr_file->get_phrase_str (&nx, parent);
      if (rc)
        {
          fprintf (stderr, "Error: W_BRANCH_F can not read parent name\n");
          return -1;
        }
      //read md
      rc = fr_file->get_phrase_double (&nx, &md);
      if (rc)
        {
          fprintf (stderr, "Error: W_BRANCH_F\n");
          return -1;
        }
      // read file name
      fname[0] = '\0';
      rc = fr_file->get_phrase_filepath (&nx, fname);
      if (rc)
        {
          fprintf (stderr, "Error: W_BRANCH_F can not read file name\n");
          return -1;
        }
      // read well_log file name
      fname2[0] = '\0';
      rc = fr_file->get_phrase_filepath (&nx, fname2);
      if (rc)
        {
          fprintf (stderr, "Error: W_BRANCH_F can not read file name\n");
          return -1;
        }
      // add to data base
      sprintf (sql, "INSERT OR REPLACE INTO branches(well_name, branch_name, md, parent) VALUES ('%s', '%s', %lf, '%s')",
               wname, branch, md, parent);
      printf ("SQL: %s\n", sql);
      if (exec_sql (sql))
        return -1;
      if (fname[0] != '\0')
        {
          sp_traj_t sp_traj = BS_KERNEL.create_object ("traj");
          rc = sp_traj->read_from_dev_file (fr_file->get_incdir() + std::string(fname));
          if (rc)
            {
              return rc;
            }
          if (add_branch_traj (wname, branch, sp_traj))
            return -6;
          printf ("TRAJ\n %s\n", sp_traj->py_str ().c_str ());
        }
      if (fname2[0] != '\0')
        {
          sp_gis_t sp_gis = BS_KERNEL.create_object ("gis");
          rc = sp_gis->read_from_las_file (fr_file->get_incdir() + std::string(fname2));
          if (rc)
            {
              return rc;
            }
          if (add_branch_gis (wname, branch, sp_gis))
            return -6;
        }

      printf ("W_BRANCH_F %s %s %s %lf\n", wname, branch, parent, md);

      return 0;
    }
  int
  sql_well::read_w_spec (char *buf)
    {
      int rc = 0;
      char *nx = 0;
      char wname[1024];
      char sql[1024];
      double x = -1, y = -1;
      // read well name
      wname[0] = '\0';
      nx = buf;
      rc = fr_file->get_phrase_str (&nx, wname);
      //printf ("WNAME: %s\n", wname);
      if (rc || wname[0] == '\0')
        {
          fprintf (stderr, "Error: well name in keyword W_SPEC can not be set by default\n");
          return -1;
        }
      rc = fr_file->get_phrase_double (&nx, &x);
      //printf ("rc: %lf\n", x);
      rc |= fr_file->get_phrase_double (&nx, &y);
      //printf ("rc: %lf\n", y);
      if (rc)
        {
          fprintf (stderr, "Error: W_SPEC\n");
          return -1;
        }
      printf ("W_SPEC %s %lf %lf\n", wname, x, y);
      // add to data base
      sprintf (sql, "INSERT INTO wells(name, x, y) VALUES ('%s', %lf, %lf)",
               wname, x, y);
      return exec_sql (sql);
    }
  int
  sql_well::read_w_frac (char *buf, double d)
    {
      int rc = 0;
      char *nx = 0;
      char wname[1024];
      char branch[1024];
      char status[1024];
      int i_status;
      char sql[1024];
      double md = -1;
      double angle  = 0;
      double half_length_1 = 50.0;
      double half_length_2 = 50.0;
      double half_up = 5.0;
      double half_down = 5.0;
      double perm = -1;
      double half_thin = 0.005;
      // read well name
      wname[0] = '\0';
      nx = buf;
      rc = fr_file->get_phrase_str (&nx, wname);
      //printf ("WNAME: %s\n", wname);
      if (rc || wname[0] == '\0')
        {
          fprintf (stderr, "Error: well name in keyword W_FRAC can not be set by default\n");
          return -1;
        }
      strcpy (branch, "main");
      rc = fr_file->get_phrase_str (&nx, branch);
      //printf ("WNAME: %s\n", wname);
      if (rc)
        {
          fprintf (stderr, "Error: well branch name in keyword W_FRAC\n");
          return -1;
        }
      strcpy (status, "SHUT");
      rc = fr_file->get_phrase_str (&nx, status);
      //printf ("WNAME: %s\n", wname);
      if (rc)
        {
          fprintf (stderr, "Error: well status in keyword W_FRAC\n");
          return -1;
        }

      rc = fr_file->get_phrase_double (&nx, &md);
      rc |= fr_file->get_phrase_double (&nx, &angle);
      rc |= fr_file->get_phrase_double (&nx, &half_length_1);
      rc |= fr_file->get_phrase_double (&nx, &half_length_2);
      rc |= fr_file->get_phrase_double (&nx, &half_up);
      rc |= fr_file->get_phrase_double (&nx, &half_down);
      rc |= fr_file->get_phrase_double (&nx, &perm);
      rc |= fr_file->get_phrase_double (&nx, &half_thin);
      if (rc)
        {
          fprintf (stderr, "Error: W_SPEC\n");
          return -1;
        }
      // check input data
      fr_file->locale_ucase (status);
      if (status[0] == 'S')
        i_status = STATUS_CON_SHUT;
      //else if (status[0] == 'C') // close
      //  i_status = 1;
      else if (status[0] == 'O') // OPEN
        i_status = STATUS_CON_OPEN;
      else
        {
          fprintf (stderr, "Error: status %s for W_FRAC not aloowed\n", status);
          return -9;
        }
      if (md < 0)
        {
          fprintf (stderr, "Error: you should specify md for W_FRAC \n");
          return -9;
        }
      if (half_length_1 <= 0)
        {
          fprintf (stderr, "Error: HALF_LENGTH_1 should be > 0 for W_FRAC \n");
          return -9;
        }
      if (half_length_2 <= 0)
        {
          fprintf (stderr, "Error: HALF_LENGTH_2 should be > 0 for W_FRAC \n");
          return -9;
        }
      if (half_up <= 0)
        {
          fprintf (stderr, "Error: HALF_UP should be > 0 for W_FRAC \n");
          return -9;
        }
      if (half_down <= 0)
        {
          fprintf (stderr, "Error: HALF_DOWN should be > 0 for W_FRAC \n");
          return -9;
        }
      if (half_thin <= 0)
        {
          fprintf (stderr, "Error: HALF_THIN should be > 0 for W_FRAC \n");
          return -9;
        }

      printf ("W_FRAC %s %s %s %lf %lf %lf %lf %lf %lf %lf %lf\n",
              wname, branch, status, md, angle, half_length_1, half_length_2,
              half_up, half_down, perm, half_thin);
      // add to data base
      sprintf (sql, "INSERT INTO fractures(well_name, branch_name, md, d, status, \
half_up, half_down, angle, half_length_1, half_length_2, perm, half_thin) \
VALUES ('%s', '%s', %lf, %lf, %d, %lf, %lf, %lf, %lf, %lf, %lf, %lf)",
               wname, branch, md, d, i_status, half_up, half_down, angle,
               half_length_1, half_length_2, perm, half_thin);
      return exec_sql (sql);
      return 0;
    }
  int
  sql_well::read_w_comp (char *buf, double d)
    {
      int rc = 0;
      char *nx = 0;
      char wname[1024];
      char branch[1024];
      char status[1024];
      int i_status;
      char sql[1024];
      double md = -1, length  = 1, rw = 0.08, skin = 0, khmult = 1.0;
      // read well name
      wname[0] = '\0';
      nx = buf;
      rc = fr_file->get_phrase_str (&nx, wname);
      //printf ("WNAME: %s\n", wname);
      if (rc || wname[0] == '\0')
        {
          fprintf (stderr, "Error: well name in keyword W_COMP can not be set by default\n");
          return -1;
        }
      strcpy (branch, "main");
      rc = fr_file->get_phrase_str (&nx, branch);
      //printf ("WNAME: %s\n", wname);
      if (rc)
        {
          fprintf (stderr, "Error: well branch name in keyword W_COMP\n");
          return -1;
        }
      strcpy (status, "SHUT");
      rc = fr_file->get_phrase_str (&nx, status);
      //printf ("WNAME: %s\n", wname);
      if (rc)
        {
          fprintf (stderr, "Error: well status in keyword W_COMP\n");
          return -1;
        }

      rc = fr_file->get_phrase_double (&nx, &md);
      rc |= fr_file->get_phrase_double (&nx, &length);
      rc |= fr_file->get_phrase_double (&nx, &rw);
      rc |= fr_file->get_phrase_double (&nx, &skin);
      rc |= fr_file->get_phrase_double (&nx, &khmult);
      if (rc)
        {
          fprintf (stderr, "Error: W_SPEC\n");
          return -1;
        }
      // check input data
      fr_file->locale_ucase (status);
      if (status[0] == 'S')
        i_status = STATUS_CON_SHUT;
      //else if (status[0] == 'C') // close
      //  i_status = 1;
      else if (status[0] == 'O') // OPEN
        i_status = STATUS_CON_OPEN;
      else
        {
          fprintf (stderr, "Error: status %s for W_COMP not aloowed\n", status);
          return -9;
        }
      if (md < 0)
        {
          fprintf (stderr, "Error: you should specify md for W_COMP \n");
          return -9;
        }
      if (length <= 0)
        {
          fprintf (stderr, "Error: length should be > 0 for W_COMP \n");
          return -9;
        }
      if (rw <= 0)
        {
          fprintf (stderr, "Error: rw should be > 0 for W_COMP \n");
          return -9;
        }
      if (khmult <= 0)
        {
          fprintf (stderr, "Error: khmult should be > 0 for W_COMP \n");
          return -9;
        }

      printf ("W_COMP %s %s %lf %lf %lf %lf %lf\n",
              wname, branch, md, length, rw, skin, khmult);
      // add to data base
      sprintf (sql, "INSERT INTO completions(well_name, branch_name, md, d, length, status, rw, kh_mult) VALUES ('%s', '%s', %lf, %lf, %lf, %d, %lf, %lf)",
               wname, branch, md, d, length, i_status, rw, khmult);
      return exec_sql (sql);
      return 0;
    }

  int
  sql_well::read_w_inj (char *buf, double d)
    {
      int rc = 0;
      char *nx = 0;
      char wname[1024];
      char status[1024];
      char ctrl[1024];
      char fluid[1024];
      int i_status = 0;
      int i_ctrl = 0;
      char sql[1024];
      double bhp = -1;
      double rate = -1;
      double orate = -1;
      double wrate = -1;
      double grate = -1;
      double lim_bhp = -1;
      double lim_rate = -1;
      double lim_orate = -1;
      double lim_wrate = -1;
      double lim_grate = -1;
      double wefac = 1.0;

      // read well name
      wname[0] = '\0';
      nx = buf;
      rc = fr_file->get_phrase_str (&nx, wname);
      //printf ("WNAME: %s\n", wname);
      if (rc || wname[0] == '\0')
        {
          fprintf (stderr, "Error: well name in keyword W_INJ can not be set by default\n");
          return -1;
        }
      strcpy (status, "SHUT");
      rc = fr_file->get_phrase_str (&nx, status);
      //printf ("WNAME: %s\n", wname);
      if (rc)
        {
          fprintf (stderr, "Error: well status in keyword W_INJ\n");
          return -1;
        }
      strcpy (ctrl, "BHP");
      rc = fr_file->get_phrase_str (&nx, ctrl);
      if (rc)
        {
          fprintf (stderr, "Error: well control in keyword W_INJ\n");
          return -1;
        }
      strcpy (fluid, "WATER");
      rc = fr_file->get_phrase_str (&nx, ctrl);
      if (rc)
        {
          fprintf (stderr, "Error: well control in keyword W_INJ\n");
          return -1;
        }

      rc = fr_file->get_phrase_double (&nx, &bhp);
      rc |= fr_file->get_phrase_double (&nx, &rate);
      rc = fr_file->get_phrase_double  (&nx, &lim_bhp);
      rc |= fr_file->get_phrase_double (&nx, &lim_rate);
      rc |= fr_file->get_phrase_double (&nx, &wefac);
      if (rc)
        {
          fprintf (stderr, "Error: W_PROD\n");
          return -1;
        }
      // check input data
      fr_file->locale_ucase (status);
      fr_file->locale_ucase (ctrl);
      fr_file->locale_ucase (fluid);
      if (status[0] == 'S')
        i_status = STATUS_SHUT;
      else if (status[0] == 'C') // close
        i_status = STATUS_CLOSE;
      else if (status[0] == 'O') // OPEN
        i_status = STATUS_OPEN;
      else
        {
          fprintf (stderr, "Error: status %s for W_COMP not aloowed\n", status);
          return -9;
        }
      if (ctrl[0] == 'B')
        {
          i_ctrl = CTRL_I_BHP;
          if (i_status == 2 && bhp <= 0)
            {
              fprintf (stderr, "Error: BHP = %lf for CONTROL %s in keyword W_INJ",
                       bhp, ctrl);
              return -1;
            }
        }
      else if (ctrl[0] == 'R') // rate
        {
          if (i_status == STATUS_OPEN && rate <= 0)
            {
              fprintf (stderr, "Error: RATE = %lf for CONTROL %s in keyword W_INJ",
                       wrate, ctrl);
              return -1;
            }
          if (fluid[0] == 'W')
            {
              i_ctrl = CTRL_I_WRATE;
              wrate = rate;
              lim_wrate = lim_rate;
            }
          else if (fluid[0] == 'O')
            {
              i_ctrl = -CTRL_I_ORATE;
              orate = rate;
              lim_orate = lim_rate;
            }
          else if (fluid[0] == 'G')
            {
              i_ctrl = CTRL_I_GRATE;
              grate = rate;
              lim_grate = lim_rate;
            }
          else
            {
              fprintf (stderr, "Error: FLUID = %s not allowed in keyword W_INJ",
                       fluid);
              return -1;
            }
        }
      else
        {
          fprintf (stderr, "Error: control %s for W_INJ not aloowed\n", ctrl);
          return -9;
        }
      if (i_status == STATUS_OPEN && wefac <= 0)
        {
          fprintf (stderr, "Error: WEFAC = %lf in keyword W_INJ",
                   wefac);
          return -1;
        }

      printf ("W_INJ %s %s %s %s %lf %lf %lf %lf %lf\n",
              wname, status, ctrl, fluid, bhp, rate, lim_bhp, lim_rate, wefac);
      // add to data base
      sprintf (sql, "INSERT INTO well_hist(well_name, d, i_or, i_wr, i_gr, \
i_bhp, wefac, ctrl, status, lim_i_or, lim_i_wr, lim_i_gr, lim_i_bhp) \
VALUES ('%s', %lf, %lf, %lf, %lf, %lf, %lf, %d, %d, %lf, %lf, %lf, %lf)",
               wname, d, orate, wrate, grate, bhp, wefac, i_ctrl, i_status,
               lim_orate, lim_wrate, lim_grate, lim_bhp);
      return exec_sql (sql);
      return 0;
    }

  int
  sql_well::read_w_prod (char *buf, double d)
    {
      int rc = 0;
      char *nx = 0;
      char wname[1024];
      char status[1024];
      char ctrl[1024];
      int i_status = 0;
      int i_ctrl = 0;
      char sql[1024];
      double bhp = -1;
      double orate = -1;
      double wrate = -1;
      double grate = -1;
      double lrate = -1;
      double lim_bhp = -1;
      double lim_orate = -1;
      double lim_wrate = -1;
      double lim_grate = -1;
      double lim_lrate = -1;
      double wefac = 1.0;

      // read well name
      wname[0] = '\0';
      nx = buf;
      rc = fr_file->get_phrase_str (&nx, wname);
      //printf ("WNAME: %s\n", wname);
      if (rc || wname[0] == '\0')
        {
          fprintf (stderr, "Error: well name in keyword W_PROD can not be set by default\n");
          return -1;
        }
      strcpy (status, "SHUT");
      rc = fr_file->get_phrase_str (&nx, status);
      //printf ("WNAME: %s\n", wname);
      if (rc)
        {
          fprintf (stderr, "Error: well status in keyword W_PROD\n");
          return -1;
        }
      strcpy (ctrl, "BHP");
      rc = fr_file->get_phrase_str (&nx, ctrl);
      if (rc)
        {
          fprintf (stderr, "Error: well control in keyword W_PROD\n");
          return -1;
        }

      rc = fr_file->get_phrase_double (&nx, &bhp);
      rc |= fr_file->get_phrase_double (&nx, &wrate);
      rc |= fr_file->get_phrase_double (&nx, &orate);
      rc |= fr_file->get_phrase_double (&nx, &grate);
      rc |= fr_file->get_phrase_double (&nx, &lrate);
      rc = fr_file->get_phrase_double  (&nx, &lim_bhp);
      rc |= fr_file->get_phrase_double (&nx, &lim_wrate);
      rc |= fr_file->get_phrase_double (&nx, &lim_orate);
      rc |= fr_file->get_phrase_double (&nx, &lim_grate);
      rc |= fr_file->get_phrase_double (&nx, &lim_lrate);

      rc |= fr_file->get_phrase_double (&nx, &wefac);
      if (rc)
        {
          fprintf (stderr, "Error: W_PROD\n");
          return -1;
        }
      // check input data
      fr_file->locale_ucase (status);
      fr_file->locale_ucase (ctrl);
      if (status[0] == 'S')
        i_status = STATUS_SHUT;
      else if (status[0] == 'C') // close
        i_status = STATUS_CLOSE;
      else if (status[0] == 'O') // OPEN
        i_status = STATUS_OPEN;
      else
        {
          fprintf (stderr, "Error: status %s for W_COMP not aloowed\n", status);
          return -9;
        }
      if (ctrl[0] == 'B')
        {
          i_ctrl = CTRL_P_BHP;
          if (i_status == 2 && bhp <= 0)
            {
              fprintf (stderr, "Error: BHP = %lf for CONTROL %s in keyword W_PROD",
                       bhp, ctrl);
              return -1;
            }
        }
      else if (ctrl[0] == 'W') // close
        {
          i_ctrl = CTRL_P_WRATE;
          if (i_status == 2 && wrate <= 0)
            {
              fprintf (stderr, "Error: WRATE = %lf for CONTROL %s in keyword W_PROD",
                       wrate, ctrl);
              return -1;
            }
        }
      else if (ctrl[0] == 'O') // close
        {
          i_ctrl = CTRL_P_ORATE;
          if (i_status == 2 && orate <= 0)
            {
              fprintf (stderr, "Error: ORATE = %lf for CONTROL %s in keyword W_PROD",
                       orate, ctrl);
              return -1;
            }
        }
      else if (ctrl[0] == 'G') // close
        {
          i_ctrl = CTRL_P_GRATE;
          if (i_status == 2 && grate <= 0)
            {
              fprintf (stderr, "Error: GRATE = %lf for CONTROL %s in keyword W_PROD",
                       grate, ctrl);
              return -1;
            }
        }
      else if (ctrl[0] == 'L') // close
        {
          i_ctrl = CTRL_P_LRATE;
          if (i_status == 2 && lrate <= 0)
            {
              fprintf (stderr, "Error: LRATE = %lf for CONTROL %s in keyword W_PROD",
                       lrate, ctrl);
              return -1;
            }
        }
      else
        {
          fprintf (stderr, "Error: control %s for W_PROD not aloowed\n", status);
          return -9;
        }
      if (i_status == STATUS_OPEN && wefac <= 0)
        {
          fprintf (stderr, "Error: WEFAC = %lf in keyword W_PROD",
                   wefac);
          return -1;
        }

      printf ("W_PROD %s %s %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
              wname, status, ctrl, bhp, wrate, orate, grate, lrate, lim_bhp,
              lim_wrate, lim_orate, lim_grate, lim_lrate, wefac);
      // add to data base
      sprintf (sql, "INSERT INTO well_hist(well_name, d, p_or, p_wr, p_gr, p_lr, \
p_bhp, wefac, ctrl, status, lim_p_or, lim_p_wr, lim_p_gr, lim_p_lr, lim_p_bhp) \
VALUES ('%s', %lf, %lf, %lf, %lf, %lf, %lf, %lf, %d, %d, %lf, %lf, %lf, %lf, %lf)",
               wname, d, orate, wrate, grate, lrate, bhp, wefac, i_ctrl, i_status,
               lim_orate, lim_wrate, lim_grate, lim_lrate, lim_bhp);
      return exec_sql (sql);
      return 0;
    }

  int
  sql_well::read_date_and_time (char *buf, char **next_start, double *dd)
    {
      int rc = 0;
      double days;
      double t;
      char *nx;

      *next_start = buf;
      fr_file->trim_left (next_start);
      nx = *next_start;
      rc |= fr_file->get_phrase (&nx);
      if (**next_start == '*')
        days = 0;
      else
        {
          rc |= fr_file->get_dt ()->cstr2d (*next_start, days);
          *dd = (double)days;
        }

      *next_start = nx;
      fr_file->trim_left (next_start);
      nx = *next_start;
      rc |= fr_file->get_phrase (&nx);
      if (**next_start == '*')
        t = 0;
      else
        {
          rc |= fr_file->get_dt ()->cstr2t (*next_start, t);
        }
      *next_start = nx;

      *dd += (double)t;
      return rc;
    }

  int
  sql_well::save_to_bos_ascii_file (const std::string &fname, sp_pool_t pool, sp_prop_t prop)
    {
      FILE *fp = fopen (fname.c_str (), "w");
      std::list<double> dates;
      std::list<double>::iterator di, de;
      boost::python::list dims;
      fci::cd_storage cd;
      fci::cd_storage::iterator cdi, cde;
      fci::frac_storage ft;
      fci::frac_storage::iterator fti, fte;
      int w_spec_flag = 0;
      char s_buf[2048];
      int nx, ny, nz;
      // interface to mesh
      sp_himesh himesh = BS_KERNEL.create_object("handy_mesh_iface");
      BS_ASSERT(himesh);
      //const double eps = 1e-10;
      sp_dt_t dt_t = BS_KERNEL.create_object ("dt_tools");

      if (!fp)
        {
          printf ("Error: Can not open destination file %s\n", fname.c_str ());
          return -1;
        }
      // get dates list
      std::string sql = "SELECT d FROM dates ORDER BY d ASC";

      if (prepare_sql (sql.c_str ()))
        return -1;
      for (; !step_sql ();)
        {
          double d = get_sql_real (0);
          dates.push_back (d);
        }
      finalize_sql ();
      de = dates.end ();

      dims = pool->py_get_pool_dims();
      nx = boost::python::extract<int>(dims[0]);
      ny = boost::python::extract<int>(dims[1]);
      nz = boost::python::extract<int>(dims[2]);
      unsigned long nx_ny = nx * ny;
      BS_SP (well_pool_iface) sp_wp = this;
#if 0
      fci::compdat_builder compdats (nx, ny, pool->get_fp_data("COORD"), pool->get_fp_data("ZCORN"), sp_wp);
      fci::fracture_builder fractures (nx, ny, pool->get_fp_data("COORD"), pool->get_fp_data("ZCORN"), sp_wp);
#else
      fci::compl_n_frac_builder compl_n_frac (nx, ny, pool->get_fp_data("COORD"), pool->get_fp_data("ZCORN"), sp_wp);
#endif
      for (di = dates.begin (); di != de; ++di)
        {
          char d_buf[1024];
          if (*di - int(*di) == 0) // if *di is integer
            {
              dt_t->d2ecl (*di, d_buf);
              printf ("DATE %s\n", d_buf);
              fprintf (fp, "DATES\n%s\n/\n\n", d_buf);
            }
          else
            {
        	  std::list<double>::iterator di_prev = di;
        	  di_prev--;
              double dt = *di - *di_prev;
              printf ("TSTEP %lf\n", dt);
              fprintf (fp, "TSTEP\n%lf\n/\n\n", dt);
            }
          // Building compdat to put first completion I, J to WELLSPECS
          compl_n_frac.clear ();
          cd = compl_n_frac.compl_build (*di);
          cde = cd.end();

          // wells in mesh come here
          std::set< std::string > good_wells;

          if (!w_spec_flag)
            {
              w_spec_flag = 1;
              spv_double point = BS_KERNEL.create_object (v_double::bs_type ());
              prepare_sql("SELECT COUNT(*) FROM wells");
              if (!step_sql ())
                {
                  BSOUT << "Wells count " << get_sql_int(0) << bs_end;
                  point->resize (3 * get_sql_int(0));
                }
              else
                return -1;

              finalize_sql();
              t_double *point_ptr = &(*point)[0];
              // well specs
              fprintf (fp, "WELSPECS\n");

              if (prepare_sql ("SELECT name, x, y FROM wells ORDER BY name ASC"))
                return -1;
              for (; !step_sql ();)
                {
                  point_ptr[0] = get_sql_real (1);
                  point_ptr[1] = get_sql_real (2);
                  point_ptr[2] = prop->get_f ("min_z");
                  point_ptr += 3;
                }
              finalize_sql ();

              spv_uint cells = himesh->where_is_points_2d(nx, ny, pool->get_fp_data("COORD"), pool->get_fp_data("ZCORN"), point);

              if (prepare_sql ("SELECT name FROM wells ORDER BY name ASC"))
                return -1;

              for (int i = 0; !step_sql (); i++)
                {
                  std::string s = get_sql_str (0);
                  BSOUT << "Writing well " << s << bs_end;
                  unsigned long cell = (*cells)[i];
                  if (cell >= nx_ny * (unsigned long)nz) {
                    // don't write out of mesh wells
                    BSERR << std::string("Well ") + s + "is out of mesh! Omitting from WELLSPEC section" << bs_end;
                    continue;
                    //throw bs_exception ("", "Well's X Y is out of mesh!");
                  }
                  unsigned long k1 = cell / nx_ny;
                  unsigned long j1 = (cell - k1 * nx_ny) / nx;
                  unsigned long i1 = cell - k1 * nx_ny - j1 * nx;
                  fprintf (fp, "\'%s\' \'FIELD\' %lu %lu /\n", s.c_str (), i1 + 1, j1 + 1);
                  // remember well's name for filtering COMPDATS
                  good_wells.insert(s);
                }
              fprintf (fp, "/\n\n");
              finalize_sql ();
            }


          // COMPDAT
          const double eps = 1.0e-5;
          const std::set< std::string >::const_iterator good_wells_end = good_wells.end();
          if (!cd.empty ())
          {
            fprintf (fp, "COMPDAT\n");
            cde = cd.end();
            for (cdi = cd.begin(); cdi != cde; ++cdi)
            {
              // skip out of mesh wells
              if (fabs(cdi->kh_mult) > eps && good_wells.find(cdi->well_name) != good_wells_end)
                {
                  if (cdi->status)
                    fprintf (fp, "\'%s\' %lu %lu %lu %lu \'OPEN\' 2* %lf 1* %lf 1* \'%c\' /\n", cdi->well_name.c_str(), cdi->cell_pos[0] + 1, cdi->cell_pos[1] + 1, cdi->cell_pos[2] + 1, cdi->cell_pos[3] + 1, cdi->diam, cdi->skin, cdi->dir);
                  else
                    fprintf (fp, "\'%s\' %lu %lu %lu %lu \'SHUT\' 2* %lf 1* %lf 1* \'%c\' /\n", cdi->well_name.c_str(), cdi->cell_pos[0] + 1, cdi->cell_pos[1] + 1, cdi->cell_pos[2] + 1, cdi->cell_pos[3] + 1, cdi->diam, cdi->skin, cdi->dir);
                }
            }
            fprintf (fp, "/\n\n");
          }


          // WPIMULT

          if (!cd.empty ())
          {

            cde = cd.end();
            int wpimult_exist = 0;
            for (cdi = cd.begin(); cdi != cde; ++cdi)
            {
              if (fabs(cdi->kh_mult - 1.0) > 1e-6 && std::abs(cdi->kh_mult) > eps && good_wells.find(cdi->well_name) != good_wells_end)
                {
                  if (!wpimult_exist)
                    fprintf (fp, "WPIMULT\n");
                  wpimult_exist = 1;
                  fprintf (fp, "\'%s\' %lf %lu %lu %lu /\n", cdi->well_name.c_str(), cdi->kh_mult, cdi->cell_pos[0] + 1, cdi->cell_pos[1] + 1, cdi->cell_pos[2] + 1);
                }
            }
            if (wpimult_exist)
              fprintf (fp, "/\n\n");
          }


          // FRACTURES
          ft = compl_n_frac.frac_build(*di);
          if (!ft.empty ())
          {
            char main_k_str[1024];
            char perm_str[1024];
            fprintf (fp, "FRACTURE\n");
            fte = ft.end ();
            for (fti = ft.begin (); fti != fte; ++fti)
            {
              // skip out-of-mesh wells
              if(good_wells.find(fti->well_name) == good_wells_end)
                continue;

              if (fti->frac_perm > 0)
                sprintf (perm_str, "%lf", fti->frac_perm);
              else
                sprintf (perm_str, " *   * ");

              int horiz = 0;
              std::string sql = "SELECT name, horiz FROM wells WHERE name = ";
              sql += fti->well_name.c_str ();
              if (prepare_sql (sql.c_str ()))
                return -1;
              if (!step_sql ())
                {
                  horiz = get_sql_int (1);
                }
              finalize_sql ();

              if (horiz)
                sprintf (main_k_str, "%d", fti->frac_main_k + 1);
              else
                sprintf (main_k_str, " * ");

              fprintf (fp, "\'%s\' %lu %lu %lu %lu %lf %lf %lf \'%s\' %lf %s %s %s %s /\n", fti->well_name.c_str (), fti->cell_pos[0] + 1, fti->cell_pos[1] + 1, fti->cell_pos[2] + 1,
                fti->cell_pos[3] + 1, fti->frac_half_length_1, fti->frac_angle - 90, fti->frac_skin, fti->frac_status == 0 ? "SHUT" : "OPEN", fti->frac_half_thin, perm_str, " * ", " * ", main_k_str);
            }
            fprintf (fp, "/\n\n");
          }

          // WCONINJE
          sprintf (s_buf, "SELECT * FROM well_hist WHERE d=%lf AND ctrl < 0 ORDER BY well_name ASC", *di);
          if (prepare_sql (s_buf))
            return -1;
          t_uint wconinje_flag = 0;
          for (; !step_sql ();)
            {
              if (wconinje_flag == 0)
                {
                  fprintf (fp, "WCONINJE\n");
                  wconinje_flag++;
                }
              std::string s = get_sql_str (0);

              // skip out-of-mesh wells
              if(good_wells.find(s) == good_wells_end)
                continue;

              int status = get_sql_int (14);
              int ctrl = get_sql_int (13);
              double i_or = get_sql_real (8);
              double i_wr = get_sql_real (9);
              double i_gr = get_sql_real (10);
              double i_bhp = get_sql_real (11);
              double rate = i_wr;
              char s_status[1024];
              char s_ctrl[1024];
              char s_phase[1024];
              char s_params[1024];
              if (status == STATUS_OPEN)
                sprintf (s_status, "OPEN");
              else
                sprintf (s_status, "SHUT");
              if (ctrl == CTRL_I_BHP)
                {
                  sprintf (s_ctrl, "BHP");
                  // TODO: add injection phase into DB
                  sprintf (s_phase, "WATER");
                  sprintf (s_params, "2* %lf", i_bhp);
                }
              else
                {
                  sprintf (s_ctrl, "RATE");
                  if (ctrl == CTRL_I_WRATE)
                    sprintf (s_phase, "WATER");
                  else if (ctrl == CTRL_I_ORATE)
                    {
                      sprintf (s_phase, "OIL");
                      rate = i_or;
                    }
                  else if (ctrl == CTRL_I_GRATE)
                    {
                      sprintf (s_phase, "GAS");
                      rate = i_gr;
                    }
                  else
                    sprintf (s_phase, "WATER");

                  sprintf (s_params, "%lf 1* %lf", rate, i_bhp);
                }
              fprintf (fp, "\'%s\' \'%s\' \'%s\' \'%s\' %s ", s.c_str (), s_phase, s_status, s_ctrl, s_params);
              /*
              if (rate < 0)
                fprintf (fp, "2* ");
              else
                fprintf (fp, "%lf 1* ", rate);
              if (i_bhp < 0)
                fprintf (fp, "1* ");
              else
                fprintf (fp, "%lf ", i_bhp);
              */
              fprintf (fp, "/\n");
            }
          if (wconinje_flag)
            fprintf (fp, "/\n\n");
          finalize_sql ();

          // WCONPROD
          sprintf (s_buf, "SELECT * FROM well_hist WHERE d=%lf AND ctrl > 0 ORDER BY well_name ASC", *di);
          if (prepare_sql (s_buf))
            return -1;

          t_uint wconprod_flag = 0;
          for (; !step_sql ();)
            {
              if (wconprod_flag == 0)
                {
                  fprintf (fp, "WCONPROD\n");
                  wconprod_flag++;
                }
              std::string s = get_sql_str (0);

              // skip out-of-mesh wells
              if(good_wells.find(s) == good_wells_end)
                continue;

              int status = get_sql_int (14);
              int ctrl = get_sql_int (13);
              double p_or = get_sql_real (2);
              double p_wr = get_sql_real (3);
              double p_gr = get_sql_real (4);
              double p_lr = get_sql_real (5);
              double p_bhp = get_sql_real (6);
              char s_status[1024];
              char s_ctrl[1024];
              char s_params[1024];
              if (status == STATUS_OPEN)
                sprintf (s_status, "OPEN");
              else
                sprintf (s_status, "SHUT");
              if (ctrl == CTRL_P_LRATE)
                {
                  sprintf (s_ctrl, "LRAT");
                  sprintf (s_params, "3* %lf 1* %lf", p_lr, p_bhp);
                }
              else if (ctrl == CTRL_P_BHP)
                {
                  sprintf (s_ctrl, "BHP");
                  sprintf (s_params, "5* %lf", p_bhp);
                }
              else if (ctrl == CTRL_P_ORATE)
                {
                  sprintf (s_ctrl, "ORAT");
                  sprintf (s_params, "%lf 4* %lf", p_or, p_bhp);
                }
              else if (ctrl == CTRL_P_WRATE)
                {
                  sprintf (s_ctrl, "WRAT");
                  sprintf (s_params, "1* %lf 3* %lf", p_wr, p_bhp);
                }
              else if (ctrl == CTRL_P_GRATE)
                {
                  sprintf (s_ctrl, "GRAT");
                  sprintf (s_params, "2* %lf 2* %lf", p_gr, p_bhp);
                }

              fprintf (fp, "\'%s\' \'%s\' \'%s\' %s ", s.c_str (), s_status, s_ctrl, s_params);
              /*
              if (p_or < 0)
                fprintf (fp, "1* ");
              else
                fprintf (fp, "%lf ", p_or);
              if (p_wr < 0)
                fprintf (fp, "1* ");
              else
                fprintf (fp, "%lf ", p_wr);
              if (p_gr < 0)
                fprintf (fp, "1* ");
              else
                fprintf (fp, "%lf ", p_gr);
              if (p_lr < 0)
                fprintf (fp, "2* ");
              else
                fprintf (fp, "%lf * ", p_lr);
              if (p_bhp < 0)
                fprintf (fp, "1* ");
              else
                fprintf (fp, "%lf ", p_bhp);
              */
              fprintf (fp, "/\n");
            }


          if (wconprod_flag)
            fprintf (fp, "/\n\n");
          finalize_sql ();


          // WEFAC
          sprintf (s_buf, "SELECT well_name, wefac FROM well_hist WHERE d=%lf ORDER BY well_name ASC", *di);
          if (prepare_sql (s_buf))
            return -1;

          t_uint wefac_flag = 0;
          for (; !step_sql ();)
            {
              double wefac = get_sql_real (1);
              if (wefac == 1.0)
                continue;
              if (wefac_flag == 0)
                {
                  fprintf (fp, "WEFAC\n");
                  wefac_flag++;
                }
              std::string s = get_sql_str (0);

              // skip out-of-mesh wells
              if(good_wells.find(s) == good_wells_end)
                continue;

              fprintf (fp, "\'%s\' %lf", s.c_str (), wefac);
              fprintf (fp, "/\n");
            }


          if (wefac_flag)
            fprintf (fp, "/\n\n");
          finalize_sql ();
        }

      fclose (fp);
      return 0;
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
