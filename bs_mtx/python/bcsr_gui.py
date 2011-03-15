# coding=utf-8

import numpy as np
from enthought.traits.api import HasTraits, Property, Array, Long, Instance, Float

from enthought.traits.ui.api import View, Item, TabularEditor, Group, HSplit

from enthought.traits.ui.tabular_adapter import TabularAdapter

from enthought.traits.ui.ui_editors.array_view_editor import ArrayViewEditor


import bs

class bcsr_gui (HasTraits):
  n_rows = Long (0)
  n_cols = Long (0)
  n_mem  = Float (0)
  n_bs   = Long (0)
  rows   = Array
  cols   = Array
  values   = Array
  #matrix = Instance (tt)

  view = View (Group (Item ('n_rows', style = 'readonly', label = 'Number of rows'),
                      Item ('n_cols', style = 'readonly', label = 'Number of columns'),
                      Item ('n_bs', style = 'readonly', label = 'Block size'),
                      Item ('n_mem', style = 'readonly', label = 'Total allocated memory (MB)'),
                      HSplit (Item ('rows', show_label = False, style = 'readonly', 
                                    editor = ArrayViewEditor (titles = ['rows_ptr'],
                                                              format = '%d',
                                                              font = 'Arial 8')),
                              Item ('cols', show_label = False, style = 'readonly', 
                                    editor = ArrayViewEditor (titles = ['cols_ind'],
                                                              format = '%d',
                                                              font = 'Arial 8')),
                              Item ('values', show_label = False, style = 'readonly', 
                                    editor = ArrayViewEditor (titles = ['values'],
                                                              format = '%lf',
                                                              font = 'Arial 8'))),
                      label = 'Matrix properties:',
                      show_border = True ),

               title    = 'Block CSR matrix viewer',
               width     = 0.5,
               height    = 0.5,
               resizable = True

              )



if __name__ == '__main__':
  m = bs.mx.bcsr_matrix ()
  mt = bs.mx.bcsr_matrix_tools ()
  mt.random_init (m, 1000000, 1, 2, 27)
  print m

  b = bcsr_gui (n_rows = m.n_rows, 
                n_cols = m.n_cols, 
                n_mem = m.allocated_memory_in_mbytes, 
                n_bs = m.n_block_size, 
                #matrix = t
                rows = m.get_rows_ptr (),
                cols = m.get_cols_ind (),
                values = m.get_values ())
  b.configure_traits ()
