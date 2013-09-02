# coding=utf-8
import sys
import bs


if __name__ == '__main__':
  db = bs.sql_well.sql_well ()
  db.open_db (sys.argv[1])
  #db.create_db_struct ()
  db.fill_db ()

  #test add_well
  db.add_well ("hello_1")
  db.add_well ("42")
  print db.get_well_names ()

  db.close_db ()
  
