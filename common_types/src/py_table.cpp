/** 
 * @file py_table.cpp
 * @brief Python #table and #table_iface interface
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-02-26
 */
#ifdef BSPY_EXPORTING_PLUGIN
#include <boost/python.hpp>
#endif

#include "table_iface.h"
#include "table.h"
#include "export_python_wrapper.h"

using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {
  //////////////////////////////////////////////////////////////////////////
  t_long (table_iface::*add_col_vector1)(const std::wstring&) = &table_iface::add_col_vector;
  void   (table_iface::*add_col_vector2)(t_long, std::wstring const&, spv_double) =
    &table_iface::add_col_vector;

  PY_EXPORTER (py_table_exporter, default_exporter)
    .def ("init",                               &T::init, 
        args ("n_rows", "n_cols"), "Initialize table n_rows * n_cols")
    .def ("clear",                               &T::clear, 
        args (""), "Clear all")
    .def ("copy",                               &T::copy, 
        args ("table"), "Copy from another table")
    .def ("set_col_name",                       &T::set_col_name, 
        args ("col_idx", "name"), "Set name for the col_idx column")
    .def ("get_col_name",                       &T::get_col_name, 
        args ("col_idx"), "Return name of the col_idx column")
    .def ("get_col_names",                      &T::get_col_names, 
        args (""), "Return list of column names")
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
    .def ("remove_row",                            &T::remove_row, 
        args ("row_idx"), "Remove row at the position row")
    .def ("get_value",                            &T::get_value, 
        args ("row_idx", "column_idx"), "Return value at given row and column")
    .def ("set_value",                            &T::set_value, 
        args ("row_idx", "column_idx", "value"), "Set value at given row and column")
    .def ("check_serial",                         &T::check_serial, 
        args (""), "should return copy of the table")
    .def ("to_str",                               &T::to_str, 
        args (""), "serialize class content to string")
    .def ("from_str",                             &T::from_str, 
        args ("string"), "restore class content from string")
    .def ("__str__",                            &T::py_str)
// uentity: newly exported methods
    .def("get_col_vector", &T::get_col_vector,
        boost::python::return_internal_reference<>())
    .def("add_col_vector", add_col_vector1)
    .def("add_col_vector", add_col_vector2)
  PY_EXPORTER_END;                               

  //////////////////////////////////////////////////////////////////////////
  //! export matrices to python
  void py_export_table ()
  {
    using namespace boost::python;

    //base_exporter <matrix_iface <seq_vector<float>, seq_vector<int> >,   default_exporter>::export_class ("matrix_iface_fi"); 
    //matrix_exporter::export_base <matrix_iface, py_matrix_iface_exporter> ("matrix_iface");
    base_exporter <table_iface, py_table_exporter>::export_class ("table_iface");

    class_exporter <table, table_iface, py_table_exporter>::export_class ("table");

  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
