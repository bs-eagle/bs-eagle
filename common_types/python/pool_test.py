# coding=utf-8
import sys
import numpy as np
import bs
import time as tm
import random as rnd



def test_pool (fname):
  n = 100

  p = bs.comm.h5_pool ()
  p.open_file (fname, "/pool")
  for i in xrange (n):
    a = np.random.random ((100, 100, 10))
    p.set_fp_data (str(i), a)

  b = np.ones ((10, 10), 'int32')
  p.set_fp_data ("pressure", a)
  p.set_i_data ('fipnum', b)
  print p
  p.flush ()
  p.close_file ()
  p.open_file (fname, "/pool")
  t_start = tm.clock ()
  for i in xrange (10):
    c = p.get_i_data ('fipnum')
    d = p.get_fp_data ('pressure')
  print "Time: ", tm.clock () - t_start

  t_start = tm.clock ()
  for i in xrange (1000):
    s_idx = str (rnd.randint (0, n - 1))
    d = p.get_fp_data (s_idx)
    #print d[49,87,9]
  print "Time: ", tm.clock () - t_start
  print p
  #print c
  #print d


if __name__ == '__main__':
  test_pool (sys.argv[1])
