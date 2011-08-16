# coding=utf-8


if __name__ == '__main__':
  import sys
  import bs
  ws = bs.sql_well.sql_well ();
  ws.open_db (sys.argv[1])

  ws.read_from_ascii_file (sys.argv[2], 10)
  ws.close_db ()
  
