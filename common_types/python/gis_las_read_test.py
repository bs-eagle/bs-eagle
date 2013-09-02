# coding=utf-8

import bs
import sys

if __name__ == '__main__':
  g = bs.comm.gis ()
  g.read_from_las_file (sys.argv[1])
  print g
