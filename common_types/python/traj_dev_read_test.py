# coding=utf-8

import bs
import sys

if __name__ == '__main__':
  g = bs.comm.traj ()
  g.read_from_dev_file (sys.argv[1])
  print g
