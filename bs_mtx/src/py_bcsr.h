/**
 * @file py_bcsr.h
 * @brief python interface for BCSR matrix
 * @author
 * @date 2009-12-07
 */
#ifndef __PY_BCSR_H
#define __PY_BCSR_H

#include "bcsr_matrix_iface.h"
#include "bcsr_amg_matrix_iface.h"
#include "bcsr.h"
#include "bcsr_matrix_tools.h"
#include "py_matrix_iface.h"

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

  PY_EXPORTER (py_bcsr_matrix_iface_exporter, py_matrix_iface_exporter)
    .def ("init_by_matrix",
          &T::init_by_matrix,
          args ("matrix"), "Initialize matrix by matrix :)")
    .def ("init",
          &T::init,
          args ("n_rows", "n_cols", "n_block_size", "n_non_zeros"), "Initialize matrix")
    .def ("init_struct",
          &T::init_struct,
          args ("n_rows", "n_cols", "n_non_zeros"), "Initialize matrix structure with out values")
    .def ("alloc_rows_ptr",
          &T::alloc_rows_ptr,
          args ("n_rows"), "Initialize rows_ptr vector")
    .def ("alloc_cols_ind",
          &T::alloc_cols_ind,
          args ("n_non_zeros"), "Initialize cols_ind vector")
    .def ("alloc_values",
          &T::alloc_values,
          args ("n_non_zeros", "n_block_size"), "Initialize values vector")
    .def ("alloc_cols_ind_and_values",
          &T::alloc_cols_ind_and_values,
          args ("n_non_zeros", "n_block_size"), "Initialize values and cols_ind vectors")
    .def ("copy",
          &T::copy,
          args ("matrix"), "Initialize matrix by matrix and copy all content")
    .def ("get_rows_ptr",
          &T::get_rows_ptr,
          args (""), "Return reference to the rows_ptr vector")
    .def ("get_cols_ind",
          &T::get_cols_ind,
          args (""), "Return reference to the cols_ind vector")
    .def ("get_values",
          &T::get_values,
          args (""), "Return reference to values vector")
    .def ("internal_check",
          &T::internal_check,
          args (""), "Return 0 if OK")
    .def ("get_n_rows",
          &T::get_n_rows,
          args (""), "Return number of matrix rows")
  PY_EXPORTER_END;

  PY_EXPORTER (py_bcsr_amg_matrix_iface_exporter, py_bcsr_matrix_iface_exporter)
    .def ("build_transpose",
          &T::build_transpose,
          args ("matrix", "rows_offset", "cols_offset", "new_n_rows"), "Build transpose matrix")
    .def ("build_transpose_struct",
          &T::build_transpose_struct,
          args ("matrix", "rows_offset", "cols_offset", "new_n_rows"), "Build transpose matrix (only structure)")
    .def ("triple_matrix_product",
          &T::triple_matrix_product,
          args ("R", "A", "P", "update_flag"), "calculate M = RAP (if update_flag == true calculate only values do not rebuild structure)")
  PY_EXPORTER_END;

  PY_EXPORTER (py_matrix_bcsr_tools_exporter, default_exporter)
    .def ("ascii_read_from_csr_format",
          &T::ascii_read_from_csr_format,
          args ("matrix", "filename"), "Read given matrix from file in ascii format")
    .def ("ascii_write_to_csr_format",
          &T::ascii_write_to_csr_format,
          args ("matrix", "filename", "sort cols flag"), "Write given matrix to file in ascii format")
    .def ("random_init",
          &T::random_init,
          args ("matrix", "n_rows", "n_block_size", "value_dispertion", "elems_in_row"),
          "Initialize matrix by random values")
    .def ("gen_2d_laplas",
          &T::gen_2d_laplas,
          args ("matrix", "dimension of problem"),
          "Initialize matrix for 2d laplas problem")
    .def ("dense_init",
          &T::dense_init,
          args ("matrix", "n_rows", "n_block_size", "value_dispertion"),
          "Initialize matrix by random values (make dense matrix)")

  PY_EXPORTER_END;


  //PY_EXPORTER (py_dummy_exporter, default_exporter)

  //PY_EXPORTER_END;

  //! export matrices to python
  void py_export_bcsr_matrices ();


  } // namespace python
} // namespace blue_sky
#endif //
#endif //__PY_BCSR_H

