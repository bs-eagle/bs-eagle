/** 
 * @file sql_well.cpp
 * @brief implementation of frac storage
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-29
 */

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
        : bs_refcounter ()
    {
      *this = rhs;
      fr_file = 0;
    }
  sql_well::~sql_well ()
    {
      if (fr_file)
        delete fr_file;
      fr_file = 0;
    }

  int 
  sql_well::open_db (const std::string &file)
    {
      int rc = 0;
      char *zErrMsg = 0;

      if (db)
        close_db ();
      printf ("SQL open_db %s\n", file.c_str ());
      if (file == ":memory:" || !boost::filesystem::exists (file))
        {
          rc = sqlite3_open (file.c_str (), &db);
          if (rc)
            {
              fprintf (stderr, "Can't open database: %s\n", sqlite3_errmsg (db));
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
              fprintf (stderr, "Can't open database: %s\n", sqlite3_errmsg (db));
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
  

  int 
  sql_well::create_db (sqlite3 *db_in)
    {
      int rc = 0;
      char *zErrMsg = 0;
      const char *sql = \
"\
BEGIN;\
CREATE TABLE wells(name TEXT UNIQUE PRIMARY KEY, \
				    x REAL DEFAULT -1, \
				    y REAL DEFAULT -1);\
CREATE TABLE groups(name TEXT UNIQUE PRIMARY KEY);\
COMMIT;\
BEGIN;\
CREATE TABLE wells_in_group(gr_name TEXT REFERENCES groups(name) ON UPDATE CASCADE ON DELETE CASCADE,\
						    well_name TEXT REFERENCES wells(name) ON UPDATE CASCADE ON DELETE CASCADE);\
CREATE INDEX i4 ON wells_in_group (gr_name ASC);\
CREATE INDEX i5 ON wells_in_group (well_name ASC);\
CREATE TABLE well_hist(well_name TEXT REFERENCES wells(name) ON UPDATE CASCADE ON DELETE CASCADE,\
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
CREATE TABLE well_res(well_name TEXT REFERENCES wells(name) ON UPDATE CASCADE ON DELETE CASCADE,\
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
COMMIT;\
BEGIN;\
CREATE TABLE branches(well_name TEXT REFERENCES wells(name) ON UPDATE CASCADE ON DELETE CASCADE,\
					   branch_name TEXT DEFAULT 'main', \
                       md REAL DEFAULT -1,\
                       parent TEXT DEFAULT '',\
					   traj BLOB, \
					   well_log BLOB,\
					   PRIMARY KEY (well_name, branch_name));\
CREATE INDEX i9 ON branches (well_name ASC);\
CREATE UNIQUE INDEX i10 ON branches (well_name, branch_name ASC);\
COMMIT;\
CREATE TRIGGER tr1 AFTER INSERT ON wells\
	BEGIN\
		INSERT INTO branches(well_name, branch_name) VALUES(new.name, 'main');\
	END;\
CREATE TRIGGER tr2 AFTER INSERT ON wells\
	BEGIN\
		INSERT INTO wells_in_group(gr_name, well_name) VALUES ('field', new.name);\
	END;\
BEGIN;\
COMMIT;\
BEGIN;\
  CREATE TABLE fractures(well_name TEXT NOT NULL,\
					     branch_name TEXT DEFAULT 'main', \
					     md REAL NOT NULL, \
					     d REAL NOT NULL,\
					     status  INTEGER DEFAULT 0,\
					     half_up REAL DEFAULT 5,\
					     half_down REAL DEFAULT 5,\
					     angle REAL DEFAULT 0,\
					     half_length_1 REAL DEFAULT 50,\
					     half_length_2 REAL DEFAULT 50,\
					     perm REAL DEFAULT -1,\
					     half_thin REAL DEFAULT 0.005,\
					     FOREIGN KEY (well_name, branch_name) REFERENCES branches(well_name, branch_name) ON UPDATE CASCADE ON DELETE CASCADE\
					     );\
CREATE INDEX i11 ON fractures (well_name ASC);\
CREATE INDEX i12 ON fractures (well_name, branch_name ASC);\
CREATE TABLE completions(well_name TEXT NOT NULL, \
					     branch_name TEXT DEFAULT 'main', \
					     md REAL NOT NULL, \
					     d REAL NOT NULL,\
					     status  INTEGER DEFAULT 0,\
					     length REAL DEFAULT 1,\
					     rw REAL DEFAULT 0.08,\
					     r0 REAL DEFAULT -1,\
					     kh REAL DEFAULT -1,\
					     kh_mult REAL DEFAULT 1,\
					     FOREIGN KEY (well_name, branch_name) REFERENCES branches(well_name, branch_name) ON UPDATE CASCADE ON DELETE CASCADE\
					     );					     \
CREATE INDEX i13 ON completions (well_name ASC);\
CREATE INDEX i14 ON completions (well_name, branch_name ASC);\
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
      
      if (stmp)
        fprintf(stderr, "stmp not null1\n");
      else
         fprintf(stderr, "stmp null1\n");
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
        return sp_gis;
      if (stmp_sql)
        return sp_gis;

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
          return sp_gis;
        }
      if (sqlite3_step (stmp) == SQLITE_ROW) // UPDATE
        {
          int n = sqlite3_column_bytes (stmp, 0);  
          if (n < 1)
            {
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
        return sp_traj;
      if (stmp_sql)
        return sp_traj;

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
          return sp_traj;
        }
      if (sqlite3_step (stmp) == SQLITE_ROW) // UPDATE
        {
          
          int n = sqlite3_column_bytes (stmp, 0);  
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
      sqlite3_finalize (stmp);
      return sp_traj;
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
          fprintf (stderr, "Can't make select: %s\n", sqlite3_errmsg (db));
          return -1;
        }
      return 0;
    }

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
      return 0;
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
  sql_well::read_from_ascii_file (const std::string &fname, double starting_date)
    {
      if (fr_file)
        delete fr_file;
      fr_file = new FRead (fname.c_str (), fname.c_str ());
      
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
          rc = trim_left (&id);
          other = id;
          rc |= get_phrase (&other);
          if (rc)
            return -1;
          //trim_right_s (id);
          locale_ucase (id);

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
      rc = get_phrase_str (&nx, wname);
      if (rc || wname[0] == '\0')
        {
          fprintf (stderr, "Error: well name in keyword W_BRANCH_F can not be set by default\n");
          return -1;
        }
      // read branch
      strcpy (branch, "main");
      rc = get_phrase_str (&nx, branch);
      if (rc)
        {
          fprintf (stderr, "Error: W_BRANCH_F can not read branch name\n");
          return -1;
        }
      // read parent
      parent[0] = '\0';
      rc = get_phrase_str (&nx, parent);
      if (rc)
        {
          fprintf (stderr, "Error: W_BRANCH_F can not read parent name\n");
          return -1;
        }
      //read md
      rc = get_phrase_double (&nx, &md);
      if (rc)
        {
          fprintf (stderr, "Error: W_BRANCH_F\n");
          return -1;
        }
      // read file name
      fname[0] = '\0';
      rc = get_phrase_filepath (&nx, fname);
      if (rc)
        {
          fprintf (stderr, "Error: W_BRANCH_F can not read file name\n");
          return -1;
        }
      // read well_log file name
      fname2[0] = '\0';
      rc = get_phrase_filepath (&nx, fname2);
      if (rc)
        {
          fprintf (stderr, "Error: W_BRANCH_F can not read file name\n");
          return -1;
        }
      // add to data base
      sprintf (sql, "INSERT OR REPLACE INTO branches(well_name, branch_name, md, parent) VALUES ('%s', '%s', %lf, '%s')",
               wname, branch, md, parent);
      if (exec_sql (sql))
        return -1;
      if (fname[0] != '\0')
        {
          sp_traj_t sp_traj = BS_KERNEL.create_object ("traj"); 
          rc = sp_traj->read_from_dev_file (std::string(fr_file->get_incdir()) + std::string(fname));
          if (rc)
            {
              return rc;
            }
          if (add_branch_traj (wname, branch, sp_traj))
            return -6;
        }
      if (fname2[0] != '\0')
        {
          sp_gis_t sp_gis = BS_KERNEL.create_object ("gis"); 
          rc = sp_gis->read_from_las_file (std::string(fr_file->get_incdir()) + std::string(fname2));
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
      rc = get_phrase_str (&nx, wname);
      //printf ("WNAME: %s\n", wname);
      if (rc || wname[0] == '\0')
        {
          fprintf (stderr, "Error: well name in keyword W_SPEC can not be set by default\n");
          return -1;
        }
      rc = get_phrase_double (&nx, &x);
      //printf ("rc: %lf\n", x); 
      rc |= get_phrase_double (&nx, &y);
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
      rc = get_phrase_str (&nx, wname);
      //printf ("WNAME: %s\n", wname);
      if (rc || wname[0] == '\0')
        {
          fprintf (stderr, "Error: well name in keyword W_FRAC can not be set by default\n");
          return -1;
        }
      strcpy (branch, "main");
      rc = get_phrase_str (&nx, branch);
      //printf ("WNAME: %s\n", wname);
      if (rc)
        {
          fprintf (stderr, "Error: well branch name in keyword W_FRAC\n");
          return -1;
        }
      strcpy (status, "SHUT");
      rc = get_phrase_str (&nx, status);
      //printf ("WNAME: %s\n", wname);
      if (rc)
        {
          fprintf (stderr, "Error: well status in keyword W_FRAC\n");
          return -1;
        }
      
      rc = get_phrase_double (&nx, &md);
      rc |= get_phrase_double (&nx, &angle);
      rc |= get_phrase_double (&nx, &half_length_1);
      rc |= get_phrase_double (&nx, &half_length_2);
      rc |= get_phrase_double (&nx, &half_up);
      rc |= get_phrase_double (&nx, &half_down);
      rc |= get_phrase_double (&nx, &perm);
      rc |= get_phrase_double (&nx, &half_thin);
      if (rc)
        {
          fprintf (stderr, "Error: W_SPEC\n");
          return -1;
        }
      // check input data
      locale_ucase (status);
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
      rc = get_phrase_str (&nx, wname);
      //printf ("WNAME: %s\n", wname);
      if (rc || wname[0] == '\0')
        {
          fprintf (stderr, "Error: well name in keyword W_COMP can not be set by default\n");
          return -1;
        }
      strcpy (branch, "main");
      rc = get_phrase_str (&nx, branch);
      //printf ("WNAME: %s\n", wname);
      if (rc)
        {
          fprintf (stderr, "Error: well branch name in keyword W_COMP\n");
          return -1;
        }
      strcpy (status, "SHUT");
      rc = get_phrase_str (&nx, status);
      //printf ("WNAME: %s\n", wname);
      if (rc)
        {
          fprintf (stderr, "Error: well status in keyword W_COMP\n");
          return -1;
        }
      
      rc = get_phrase_double (&nx, &md);
      rc |= get_phrase_double (&nx, &length);
      rc |= get_phrase_double (&nx, &rw);
      rc |= get_phrase_double (&nx, &skin);
      rc |= get_phrase_double (&nx, &khmult);
      if (rc)
        {
          fprintf (stderr, "Error: W_SPEC\n");
          return -1;
        }
      // check input data
      locale_ucase (status);
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
      rc = get_phrase_str (&nx, wname);
      //printf ("WNAME: %s\n", wname);
      if (rc || wname[0] == '\0')
        {
          fprintf (stderr, "Error: well name in keyword W_INJ can not be set by default\n");
          return -1;
        }
      strcpy (status, "SHUT");
      rc = get_phrase_str (&nx, status);
      //printf ("WNAME: %s\n", wname);
      if (rc)
        {
          fprintf (stderr, "Error: well status in keyword W_INJ\n");
          return -1;
        }
      strcpy (ctrl, "BHP");
      rc = get_phrase_str (&nx, ctrl);
      if (rc)
        {
          fprintf (stderr, "Error: well control in keyword W_INJ\n");
          return -1;
        }
      strcpy (fluid, "WATER");
      rc = get_phrase_str (&nx, ctrl);
      if (rc)
        {
          fprintf (stderr, "Error: well control in keyword W_INJ\n");
          return -1;
        }
      
      rc = get_phrase_double (&nx, &bhp);
      rc |= get_phrase_double (&nx, &rate);
      rc = get_phrase_double  (&nx, &lim_bhp);
      rc |= get_phrase_double (&nx, &lim_rate);
      rc |= get_phrase_double (&nx, &wefac);
      if (rc)
        {
          fprintf (stderr, "Error: W_PROD\n");
          return -1;
        }
      // check input data
      locale_ucase (status);
      locale_ucase (ctrl);
      locale_ucase (fluid);
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
      rc = get_phrase_str (&nx, wname);
      //printf ("WNAME: %s\n", wname);
      if (rc || wname[0] == '\0')
        {
          fprintf (stderr, "Error: well name in keyword W_PROD can not be set by default\n");
          return -1;
        }
      strcpy (status, "SHUT");
      rc = get_phrase_str (&nx, status);
      //printf ("WNAME: %s\n", wname);
      if (rc)
        {
          fprintf (stderr, "Error: well status in keyword W_PROD\n");
          return -1;
        }
      strcpy (ctrl, "BHP");
      rc = get_phrase_str (&nx, ctrl);
      if (rc)
        {
          fprintf (stderr, "Error: well control in keyword W_PROD\n");
          return -1;
        }
      
      rc = get_phrase_double (&nx, &bhp);
      rc |= get_phrase_double (&nx, &wrate);
      rc |= get_phrase_double (&nx, &orate);
      rc |= get_phrase_double (&nx, &grate);
      rc |= get_phrase_double (&nx, &lrate);
      rc = get_phrase_double  (&nx, &lim_bhp);
      rc |= get_phrase_double (&nx, &lim_wrate);
      rc |= get_phrase_double (&nx, &lim_orate);
      rc |= get_phrase_double (&nx, &lim_grate);
      rc |= get_phrase_double (&nx, &lim_lrate);

      rc |= get_phrase_double (&nx, &wefac);
      if (rc)
        {
          fprintf (stderr, "Error: W_PROD\n");
          return -1;
        }
      // check input data
      locale_ucase (status);
      locale_ucase (ctrl);
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
      rc = trim_left (next_start);
      nx = *next_start;
      rc |= get_phrase (&nx);
      if (**next_start == '*')
        days = 0;
      else
        {
          rc |= key_read_date (*next_start, days);
          *dd = (double)days;
        }

      *next_start = nx;
      rc = trim_left (next_start);
      nx = *next_start;
      rc |= get_phrase (&nx);
      if (**next_start == '*')
        t = 0;
      else
        {
          rc |= key_read_time (*next_start, t);
        }
      *next_start = nx;

      *dd += (double)t;
      return rc;
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
