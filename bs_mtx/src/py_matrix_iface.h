#ifndef PY_MATRIX_IFACE_H_
#define PY_MATRIX_IFACE_H_
#include "matrix_iface.h"
#include "bdiag_matrix.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "dummy_base.h"
#include "construct_python_object.h"
#include "make_me_happy.h"
#include "python/py_bs_object_base.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include "export_python_wrapper.h"

#ifdef BSPY_EXPORTING_PLUGIN
namespace blue_sky
  {
  namespace python
    {

  PY_EXPORTER (py_matrix_iface_exporter, default_exporter)
    .add_property ("n_block_size",
                   &T::get_n_block_size,
                   "Size of block for storing block matricies")
    .add_property ("n_rows",
                   &T::get_n_rows,
                   "Number of rows in matrix")
    .add_property ("n_cols",
                   &T::get_n_cols,
                   "Number of columns in matrix")
    .add_property ("allocated_memory_in_mbytes",
                   &T::get_allocated_memory_in_mbytes,
                   "Total amount of allocated memory in MBytes")
    .def ("matrix_vector_product",
          &T::matrix_vector_product,
          args ("v", "r"), "Matrix vector product (r += Av)")
    .def ("matrix_vector_product_t",
          &T::matrix_vector_product_t,
          args ("v", "r"), "Transpose matrix vector product (r += A^(T)v)")
    .def ("calc_lin_comb",
          &T::calc_lin_comb,
          args ("alpha", "beta", "u", "v", "r"), "r = alpha * Au + beta * v")
    .def ("is_square",
          &T::is_square,
          args (""), "Return true if matrix is square")
    .def ("init_vector",
          &T::init_vector,
          args ("v"), "Initialize vector (for non MPI: call v.resize (n_rows, 0.0)")
    .def ("__str__",
          &T::py_str)

  PY_EXPORTER_END;

  PY_EXPORTER (py_bdiag_matrix_iface_exporter, py_matrix_iface_exporter)
    .def ("init_by_matrix",
          &T::init_by_matrix,
          args ("matrix"), "Initialize matrix by matrix :)")
    .def ("init",
          &T::init,
          args ("n_rows", "n_block_size"), "Initialize matrix")
    .def ("copy",
          &T::copy,
          args ("matrix"), "Initialize matrix by matrix and copy all content")
    .def ("get_diag",
          &T::get_diag,
          args (""), "Return reference to the diagonal vector")

  PY_EXPORTER_END;


  PY_EXPORTER (py_dummy_exporter, default_exporter)

  PY_EXPORTER_END;

  //! export matrices to python
  void py_export_matrices ();


  } // namespace python
} // namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
#endif // #ifndef PY_MATRIX_BASE_H_
