/**
 * \file py_jacobian_matrix.cpp
 * \brief impl of python wrapper of
 * \author Sergey Miryanov
 * \date 16.06.2008
 * */
#include "bs_matrix_stdafx.h"

#include "py_jacobian_matrix.h"

#ifdef BSPY_EXPORTING_PLUGIN

#include "jacobian_matrix.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "construct_python_object.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {
  namespace python
    {

    //template <class matrix_t>
    //py_jacobian_matrix<matrix_t>::py_jacobian_matrix ()
    //    : base_t (matrix_t::bs_type ())
    //{
    //}

    //template <class matrix_t> void
    //py_jacobian_matrix<matrix_t>::read_regular_from_csr_matrix (const char *filename)
    //{
    //  this->get_lspx (this)->read_regular_from_csr_matrix (filename);
    //}
    //template <class matrix_t> void
    //py_jacobian_matrix<matrix_t>::read_regular_from_b_matrix (const char *filename, bool summ_diag)
    //{
    //  this->get_lspx (this)->read_regular_from_b_matrix (filename, summ_diag);
    //}
    //template <class matrix_t> void
    //py_jacobian_matrix<matrix_t>::read_rhs_vector (const char *filename)
    //{
    //  this->get_lspx (this)->read_rhs_vector (filename);
    //}
    //template <class matrix_t> void
    //py_jacobian_matrix<matrix_t>::read_rhs_flux_vector (const char *filename)
    //{
    //  this->get_lspx (this)->read_rhs_flux_vector (filename);
    //}
    //template <class matrix_t> void
    //py_jacobian_matrix<matrix_t>::read_solution_vector (const char *filename)
    //{
    //  this->get_lspx (this)->read_solution_vector (filename);
    //}
    //template <class matrix_t> void
    //py_jacobian_matrix<matrix_t>::read_acc_diag_from_b_matrix (const char *filename)
    //{
    //  this->get_lspx (this)->read_acc_diag_from_b_matrix (filename);
    //}
    //template <class matrix_t> void
    //py_jacobian_matrix<matrix_t>::read_irregular_csr_matrix (const char *filename)
    //{
    //  this->get_lspx (this)->read_irregular_csr_matrix (filename);
    //}

    //template <class matrix_t> typename py_jacobian_matrix <matrix_t>::item_array_t &
    //py_jacobian_matrix<matrix_t>::get_rhs ()
    //{
    //  return this->get_lspx (this)->get_rhs ();
    //}
    //template <class matrix_t> typename py_jacobian_matrix <matrix_t>::item_array_t &
    //py_jacobian_matrix<matrix_t>::get_rhs_flux ()
    //{
    //  return this->get_lspx (this)->get_rhs_flux ();
    //}
    //template <class matrix_t>
    //typename py_jacobian_matrix <matrix_t>::item_array_t &
    //py_jacobian_matrix <matrix_t>::get_sol ()
    //{
    //  return this->get_lspx (this)->get_solution ();
    //}
    //template <class matrix_t>
    //typename py_jacobian_matrix <matrix_t>::rhs_item_array_t &
    //py_jacobian_matrix <matrix_t>::get_reg_acc_diag ()
    //{
    //  return this->get_lspx (this)->get_regular_acc_diag ();
    //}

    //template <class matrix_t>
    //typename py_jacobian_matrix <matrix_t>::item_array_t &
    //py_jacobian_matrix <matrix_t>::get_sec_rhs ()
    //{
    //  return this->get_lspx (this)->get_sec_rhs ();
    //}
    //template <class matrix_t>
    //typename py_jacobian_matrix <matrix_t>::item_array_t &
    //py_jacobian_matrix <matrix_t>::get_sec_sol ()
    //{
    //  return this->get_lspx (this)->get_sec_solution ();
    //}
    //template <class matrix_t>
    //typename py_jacobian_matrix <matrix_t>::item_array_t &
    //py_jacobian_matrix <matrix_t>::get_ss_diag ()
    //{
    //  return this->get_lspx (this)->get_ss_diagonal ();
    //}
    //template <class matrix_t>
    //typename py_jacobian_matrix <matrix_t>::item_array_t &
    //py_jacobian_matrix <matrix_t>::get_sp_diag ()
    //{
    //  return this->get_lspx (this)->get_sp_diagonal ();
    //}



    //template <class matrix_t>
    //typename py_jacobian_matrix <matrix_t>::py_bcsr_matrix_t
    //py_jacobian_matrix <matrix_t>::get_regular_matrix () const
    //{
    //  return py_bcsr_matrix_t (this->template get_spx <wrapped_t> ()->get_regular_matrix ());
    //}
    //template <class matrix_t>
    //typename py_jacobian_matrix <matrix_t>::py_bcsr_matrix_t
    //py_jacobian_matrix <matrix_t>::get_irregular_matrix () const
    //{
    //  return py_bcsr_matrix_t (this->template get_spx <wrapped_t> ()->get_irregular_matrix ());
    //}
    //template <class matrix_t>
    //typename py_jacobian_matrix <matrix_t>::py_bcsr_matrix_t
    //py_jacobian_matrix <matrix_t>::get_prepared_matrix () const
    //{
    //  return py_bcsr_matrix_t (this->template get_spx <wrapped_t> ()->get_prepared_matrix ());
    //}

    //template <typename matrix_t>
    //void
    //py_jacobian_matrix <matrix_t>::prepare_matrix ()
    //{
    //  this->get_lspx (this)->prepare_matrix ();
    //}
    //////////////////////////////////////////////////////////////////////////

    template <class matrix_t> void
    export_jacobian_matrix (const char *name)
    {
      using namespace boost::python;

      const typename matrix_t::sp_csr_matrix_t &  (matrix_t::*get_regular_matrix) () const    = &matrix_t::get_regular_matrix;
      const typename matrix_t::sp_csr_matrix_t &  (matrix_t::*get_irregular_matrix) () const  = &matrix_t::get_irregular_matrix;
      typename matrix_t::rhs_item_array_t & (matrix_t::*get_reg_acc_diag) ()                  = &matrix_t::get_regular_acc_diag;
      typename matrix_t::rhs_item_array_t & (matrix_t::*get_rhs) ()                           = &matrix_t::get_rhs;
      typename matrix_t::rhs_item_array_t & (matrix_t::*get_rhs_flux) ()                      = &matrix_t::get_rhs_flux;
      typename matrix_t::item_array_t & (matrix_t::*get_sol) ()                               = &matrix_t::get_solution;
      typename matrix_t::rhs_item_array_t & (matrix_t::*get_sec_rhs) ()                       = &matrix_t::get_sec_rhs;
      typename matrix_t::item_array_t & (matrix_t::*get_sec_sol) ()                           = &matrix_t::get_sec_solution;
      typename matrix_t::rhs_item_array_t & (matrix_t::*get_ss_diag) ()                       = &matrix_t::get_ss_diagonal;
      typename matrix_t::rhs_item_array_t & (matrix_t::*get_sp_diag) ()                       = &matrix_t::get_sp_diagonal;

      class_ <matrix_t, bases <typename matrix_t::base_t>, boost::noncopyable> (name, no_init)
        .def ("__cons__",                     make_constructor (construct_python_object <matrix_t>))
        .def ("__init__",                     make_function (init_python_object <matrix_t>))
        .def ("read_regular_from_csr_matrix", &matrix_t::read_regular_from_csr_matrix)
        .def ("read_regular_from_b_matrix",   &matrix_t::read_regular_from_b_matrix)
        .def ("read_rhs_vector",              &matrix_t::read_rhs_vector)
        .def ("read_rhs_flux_vector",         &matrix_t::read_rhs_flux_vector)
        .def ("read_solution_vector",         &matrix_t::read_solution_vector)
        .def ("read_acc_diag_from_b_matrix",  &matrix_t::read_acc_diag_from_b_matrix)
        .def ("read_irregular_csr_matrix",    &matrix_t::read_irregular_csr_matrix)
        .def ("prepare_matrix",               &matrix_t::prepare_matrix)
        .add_property ("regular_matrix",      make_function (get_regular_matrix, boost::python::return_value_policy <copy_const_reference>()))
        .add_property ("irregular_matrix",    make_function (get_irregular_matrix, boost::python::return_value_policy <copy_const_reference>()))
        .add_property ("reg_acc_diag",        make_function (get_reg_acc_diag, boost::python::return_value_policy <copy_non_const_reference>()))
        .add_property ("rhs",                 make_function (get_rhs, boost::python::return_value_policy <copy_non_const_reference>()))
        .add_property ("rhs_flux",            make_function (get_rhs_flux, boost::python::return_value_policy <copy_non_const_reference>()))
        .add_property ("sol",                 make_function (get_sol, boost::python::return_value_policy <copy_non_const_reference>()))
        .add_property ("sec_rhs",             make_function (get_sec_rhs, boost::python::return_value_policy <copy_non_const_reference>()))
        .add_property ("sec_sol",             make_function (get_sec_sol, boost::python::return_value_policy <copy_non_const_reference>()))
        .add_property ("ss_diag",             make_function (get_ss_diag, boost::python::return_value_policy <copy_non_const_reference>()))
        .add_property ("sp_diag",             make_function (get_sp_diag, boost::python::return_value_policy <copy_non_const_reference>()))
        .add_property ("prepared_mx",         make_function (&matrix_t::get_prepared_matrix))
        ;

      register_ptr_to_python <matrix_t *> ();
      register_ptr_to_python <smart_ptr <matrix_t, true> > ();
    }

    void
    py_export_jacobian_matrix ()
    {
      export_jacobian_matrix <jacobian_matrix <base_strategy_fi> >    ("jacobian_matrix_seq_fi");
      export_jacobian_matrix <jacobian_matrix <base_strategy_di> >    ("jacobian_matrix_seq_di");
      export_jacobian_matrix <jacobian_matrix <base_strategy_mixi> >  ("jacobian_matrix_seq_mixi");
    }

  } // namespace python
} // namespace blue_sky

#endif  // #ifdef BSPY_EXPORTING_PLUGIN
