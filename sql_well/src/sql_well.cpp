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
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/fstream.hpp"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>



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
      
    }
  sql_well::sql_well (const sql_well& rhs) 
        : bs_refcounter ()
    {
      *this = rhs;
    }

  int 
  sql_well::open_db (const std::string &file)
    {
      int rc = 0;
      char *zErrMsg = 0;
      char buf[2048];

      if (db)
        close_db ();
      if (!boost::filesystem::exists (file))
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
          sqlite3_close (db);
        }

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
      
      return 0;
    }
  

  int 
  sql_well::create_db (sqlite3 *db_in)
    {
      int rc = 0;
      char *zErrMsg = 0;
      const char *sql = \
"CREATE TABLE well_dynamic (well_name TEXT,\
                           date REAL NOT NULL,\
                           h_orate REAL DEFAULT -1,\
                           h_wrate REAL DEFAULT -1,\
                           h_grate REAL DEFAULT -1,\
                           h_bhp   REAL DEFAULT -1,\
                           h_wefac REAL DEFAULT -1,\
                           l_orate REAL DEFAULT -1,\
                           l_wrate REAL DEFAULT -1,\
                           l_grate REAL DEFAULT -1,\
                           l_bhp   REAL DEFAULT -1,\
                           c_orate REAL DEFAULT -1,\
                           c_wrate REAL DEFAULT -1,\
                           c_grate REAL DEFAULT -1,\
                           c_bhp   REAL DEFAULT -1,\
                           c_wefac REAL DEFAULT -1,\
                           ct_orate REAL DEFAULT -1,\
                           ct_wrate REAL DEFAULT -1,\
                           ct_grate REAL DEFAULT -1,\
                           h_status INTEGER DEFAULT -1,\
                           h_ctrl   INTEGER DEFAULT -1,\
                           c_status INTEGER DEFAULT -1,\
                           c_ctrl   INTEGER DEFAULT -1);\
CREATE INDEX well_name ON well_dynamic (well_name ASC);\
CREATE INDEX date ON well_dynamic (date ASC);\
CREATE UNIQUE INDEX iname_date ON well_dynamic (well_name, date ASC);\
CREATE TABLE wells_in_group (group_name TEXT, well_name TEXT);\
CREATE INDEX well_name2 ON wells_in_group (well_name ASC);\
CREATE INDEX group_name ON wells_in_group (group_name ASC);\
CREATE TABLE fractures (well_name TEXT, \
                        parent TEXT DEFAULT 'main', \
                        date REAL DEFAULT -1,\
                        is_vertical BOOLEAN DEFAULT 1,\
                        is_symmetric BOOLEAN DEFAULT 1,\
                        perm REAL DEFAULT -1,\
                        wf REAL DEFAULT 0.005,\
                        half_length REAL DEFAULT 50,\
                        up_half_height REAL DEFAULT 10,\
                        down_half_height REAL DEFAULT 10,\
                        angle REAL DEFAULT 0,\
                        hor_main_radius REAL DEFAULT 50,\
                        hor_sec_radius REAL DEFAULT 50,\
                        parent_md REAL DEFAULT 10);\
CREATE INDEX well_name3 ON fractures (well_name ASC);\
CREATE INDEX by_date ON fractures (well_name, date);\
CREATE INDEX duo_name ON fractures (well_name, parent);\
CREATE TABLE branches (well_name TEXT,\
                       branch_name TEXT,\
                       parent TEXT DEFAULT '',\
                       parent_md REAL DEFAULT -1.0,\
                       prop TEXT,\
                       traj TEXT,\
                       well_log TEXT);\
CREATE INDEX well_name4 ON branches (well_name ASC);\
CREATE UNIQUE INDEX duo_name2 ON branches (well_name ASC, branch_name ASC);\
CREATE INDEX parent on branches (well_name ASC, parent);\
CREATE TABLE wells (well_name TEXT,\
                    x REAL DEFAULT -1,\
                    y REAL DEFAULT -1);\
CREATE UNIQUE INDEX names ON wells (well_name ASC);\
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
          int rc = 0;
          char *zErrMsg = 0;
          char buf[2048];
          
          
          if (stmp_sql)
            finalize_sql ();
            
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

      sprintf (buf, "INSERT INTO wells (well_name) VALUES('%s'); \
INSERT INTO wells_in_group (group_name, well_name) VALUES('FIELD', '%s');\
INSERT INTO branches (well_name, branch_name) VALUES('%s', 'main')", 
               well_name.c_str (),
               well_name.c_str (),
               well_name.c_str ());
      rc = sqlite3_exec (db, buf, NULL, 0, &zErrMsg);
      if( rc != SQLITE_OK )
        {
          fprintf (stderr, "SQL error: %s\n", zErrMsg);
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
      std::string sql = "SELECT well_name FROM wells";
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
      rc = sqlite3_prepare_v2 (db, buf, strlen (buf) + 1, &stmp, &ttt);
      if (rc)
        {
          fprintf (stderr, "Can't make select: %s\n", sqlite3_errmsg (db));
          return -1;
        }

      std::ostringstream oss (std::ios_base::binary | std::ios_base::out | std::ios_base::in);
      boost::archive::binary_oarchive oar(oss);
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
          const char *b = (const char *)sqlite3_column_blob (stmp, 0);
          //std::string s = (const char *)sqlite3_column_text (stmp, 0);
          std::string s;
          s.assign (b, n);
          printf ("READ GIS %d\n", (int)s.length ());
          std::istringstream iss;
          iss.str (s);
          boost::archive::binary_iarchive iar(iss);
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
      boost::archive::binary_oarchive oar(oss);
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
          boost::archive::binary_iarchive iar(iss);
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
