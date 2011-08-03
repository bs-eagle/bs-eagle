/** 
 * @file py_table.h
 * @brief Python #table and #table_iface interface
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-02-26
 */
#ifndef PY_TABLE_PZHJ8KVY

#define PY_TABLE_PZHJ8KVY

#include <string>
#include "table_iface.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "dummy_base.h"
#include "construct_python_object.h"
#include "make_me_happy.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include "export_python_wrapper.h"

#ifdef BSPY_EXPORTING_PLUGIN
namespace blue_sky
  {
  namespace python
    {

  PY_EXPORTER (py_table_exporter, default_exporter)
    .def ("init",                               &T::init, 
        args ("n_rows", "n_cols"), "Initialize table n_rows * n_cols")
    .def ("clear",                               &T::clear, 
        args (""), "Clear all")
    .def ("set_col_name",                       &T::set_col_name, 
        args ("col_idx", "name"), "Set name for the col_idx column")
    .def ("get_col_name",                       &T::get_col_name, 
        args ("col_idx"), "Return name of the col_idx column")
    .def ("get_n_rows",                         &T::get_n_rows, 
        args (""), "Return number of rows in table")
    .def ("get_n_cols",                         &T::get_n_cols, 
        args (""), "Return number of columns in table")
    .def ("set_col_values",                     &T::set_col_values, 
        args ("col_idx", "column"), "Set values of the column col_idx")
    .def ("get_col_values",                     &T::get_col_values, 
        args ("col_idx"), "Return values of the column col_idx")
    .def ("add_row",                            &T::add_row, 
        args ("row_idx"), "Insert row at the position row")
    .def ("get_value",                            &T::get_value, 
        args ("row_idx", "column_idx"), "Return value at given row and column")
    .def ("set_value",                            &T::set_value, 
        args ("row_idx", "column_idx", "value"), "Set value at given row and column")
    .def ("check_serial",                         &T::check_serial, 
        args (""), "should return copy of the table")
    .def ("__str__",                            &T::py_str)
  PY_EXPORTER_END;                               

  //PY_EXPORTER (py_dummy_exporter, default_exporter)
                                                 
  //PY_EXPORTER_END;                               
                                                 
  //! export matrices to python                  
  void py_export_table ();                    

    
  } // namespace python
} // namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN


#endif /* end of include guard: PY_TABLE_PZHJ8KVY */
