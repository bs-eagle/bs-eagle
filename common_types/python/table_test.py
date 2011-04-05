# coding=utf-8

import numpy as np
import bs

t = bs.common_types.table ()

n_rows = 10
n_cols = 5

t.init (n_rows, n_cols);

for i in xrange (n_cols):
  t.set_col_name (i, "Col " + str (i))
  a = np.linspace (float (i), float (i + 1), n_rows)
  t.set_col_values (i, a)
print t

