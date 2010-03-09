/** 
 * \file py_bcsr_matrix.cpp
 * \brief
 * \author Sergey Miryanov
 * \date 06.05.2009
 * */

#include "bs_matrix_stdafx.h"
#include "py_bcsr_matrix.h"

#ifdef BSPY_EXPORTING_PLUGIN
#include BS_FORCE_PLUGIN_IMPORT ()
#include "construct_python_object.h"
#include BS_STOP_PLUGIN_IMPORT ()

using namespace boost::python;

namespace blue_sky {
namespace python {

  template<class matrix_t> void
  py_export_(const char * name)
  {
    int    (matrix_t::*init1)(const matrix_t &)= &matrix_t::init;
    int    (matrix_t::*init2)(const typename matrix_t::i_type_t, const typename matrix_t::i_type_t, const typename matrix_t::i_type_t, const typename matrix_t::i_type_t)= &matrix_t::init;

    typename matrix_t::fp_vector_type & (matrix_t::*get_values) () = &matrix_t::get_values;
    typename matrix_t::i_vector_type &  (matrix_t::*get_rows_ptr) () = &matrix_t::get_rows_ptr;
    typename matrix_t::i_vector_type &  (matrix_t::*get_cols_ind) () = &matrix_t::get_cols_ind;
    typename matrix_t::i_vector_type &  (matrix_t::*get_diag_ind) () = &matrix_t::get_diag_ind;

    //int (matrix_t::*mvp_df) (const seq_vector <double> &, const seq_vector <float> &)  const = reinterpret_cast <int (matrix_t::*)(const seq_vector <double> &, const seq_vector <float> &)> (&matrix_t::matrix_vector_product);
    //int (matrix_t::*mvp_ff) (const seq_vector <float> &,  const seq_vector <float> &)  const = &matrix_t::matrix_vector_product;
    //int (matrix_t::*mvp_dd) (const seq_vector <double> &, const seq_vector <double> &) const = &matrix_t::matrix_vector_product;

    class_ <matrix_t, bases <typename matrix_t::base_t>, boost::noncopyable> (name, no_init)
      .def ("__cons__",                     make_constructor (construct_python_object <matrix_t>))
      .def ("__init__",                     make_function (init_python_object <matrix_t>))
      .def("init",                          init1)
      .def("init",                          init2)
      .def("init_struct",                   &matrix_t::init_struct)
      .def("init_diag_ind",                 &matrix_t::init_diag_ind)

      .def("alloc_rows_ptr",                &matrix_t::alloc_rows_ptr)
      .def("alloc_cols_ind",                &matrix_t::alloc_cols_ind)
      .def("alloc_values",                  &matrix_t::alloc_values)
      .def("alloc_cols_ind_and_values",     &matrix_t::alloc_cols_ind_and_values)

      .def("rand_init",                     &matrix_t::rand_init)
      .def("gen_2d_laplas",                 &matrix_t::gen_2d_laplas)

      .def("build_transpose",               &matrix_t::build_transpose)
      .def("build_transpose_struct",        &matrix_t::build_transpose_struct)

      .def("internal_check",                &matrix_t::internal_check)

      .def("get_allocated_memory_in_bytes", &matrix_t::get_allocated_memory_in_bytes)
      .def("ascii_write_in_csr_format",     &matrix_t::ascii_write_in_csr_format)
      .def("ascii_write_in_ij_format",      &matrix_t::ascii_write_in_ij_format)
      .def("ascii_read_from_csr_format",    &matrix_t::ascii_read_from_csr_format)

      .def ("merge",                        &matrix_t::merge)

      .add_property ("values",              make_function (get_values, return_internal_reference <> ()))//, make_function (&matrix_t::set_values))
      .add_property ("rows",                make_function (get_rows_ptr, return_internal_reference <> ()))//, make_function (&matrix_t::set_rows))
      .add_property ("cols",                make_function (get_cols_ind, return_internal_reference <> ()))//, make_function (&matrix_t::set_cols))
      .add_property ("diags",               make_function (get_diag_ind, return_internal_reference <> ()))//, make_function (&matrix_t::set_diags))
      .add_property ("nnz",                 make_function (&matrix_t::get_n_non_zeros))
      ;
  }

  void
  py_export_base_matrices()
  {
    class_ <base_matrix_di, boost::noncopyable> ("base_matrix_di", no_init)
      .add_property ("n_block_size",        make_function (&base_matrix_di::get_n_block_size), make_function (&base_matrix_di::set_n_block_size))
      .add_property ("n_rows",              make_function (&base_matrix_di::get_n_rows), make_function (&base_matrix_di::set_n_rows))
      .add_property ("n_cols",              make_function (&base_matrix_di::get_n_cols), make_function (&base_matrix_di::set_n_cols))
      ;
    class_ <base_matrix_fi, boost::noncopyable> ("base_matrix_fi", no_init)
      .add_property ("n_block_size",        make_function (&base_matrix_fi::get_n_block_size), make_function (&base_matrix_fi::set_n_block_size))
      .add_property ("n_rows",              make_function (&base_matrix_fi::get_n_rows), make_function (&base_matrix_fi::set_n_rows))
      .add_property ("n_cols",              make_function (&base_matrix_fi::get_n_cols), make_function (&base_matrix_fi::set_n_cols))
      ;

    py_export_<bcsr_matrix_di> ("bcsr_matrix_di");
    py_export_<bcsr_matrix_fi> ("bcsr_matrix_fi");

    register_ptr_to_python <matrix_base <shared_vector <double>, shared_vector <int> > *> ();
    register_ptr_to_python <base_matrix_fi *> ();
    register_ptr_to_python <smart_ptr <base_matrix_di, true> > ();
    register_ptr_to_python <smart_ptr <base_matrix_fi, true> > ();
    register_ptr_to_python <smart_ptr <bcsr_matrix_di, true> > ();
    register_ptr_to_python <smart_ptr <bcsr_matrix_fi, true> > ();
  }

} // namespace python
} // namespace blue_sky
#endif