/**
 * \file py_jacobian_matrix.h
 * \brief python wrapper for jacobian matrix
 * \author Sergey Miryanov
 * \date 16.06.2008
 * */
#ifndef PY_BS_JACOBIAN_MATRIX_H_
#define PY_BS_JACOBIAN_MATRIX_H_

#ifdef BSPY_EXPORTING_PLUGIN

#include "py_matrix_base.h"
#include "py_bcsr_matrix.h"

namespace blue_sky
  {
  namespace python
    {

    //template <class matrix_type>
    //class py_jacobian_matrix : public py_matrix_base <typename matrix_type::rhs_item_array_t, typename matrix_type::index_array_t>
    //  {
    //  public:

    //    typedef matrix_type                                   matrix_t;
    //    typedef typename matrix_t::item_array_t               item_array_t;
    //    typedef typename matrix_t::rhs_item_array_t           rhs_item_array_t;
    //    typedef typename matrix_t::index_array_t              index_array_t;
    //    typedef py_matrix_base <rhs_item_array_t, index_array_t>  base_t;
    //    typedef matrix_t                                      wrapped_t;

    //    typedef py_bcsr_matrix<rhs_item_array_t, index_array_t>   py_bcsr_matrix_t;


    //    py_jacobian_matrix ();

    //    py_jacobian_matrix (sp_obj sp_obj_) : base_t (sp_obj_) {}

    //    //! read regular matrix from csr file
    //    void read_regular_from_csr_matrix (const char *filename);

    //    //! read regular matrix from b_matrix file
    //    void read_regular_from_b_matrix (const char *filename, bool summ_diag);

    //    //! read rhs vector from section
    //    void read_rhs_vector (const char *filename);

    //    //! read rhs_flux vector from section
    //    void read_rhs_flux_vector (const char *filename);

    //    //! read solution vector form file
    //    void read_solution_vector (const char *filename);

    //    //! read accumulative diag of regular matrix from b_matrix file
    //    void read_acc_diag_from_b_matrix (const char *filename);

    //    //! read irregular matrix from csr file
    //    void read_irregular_csr_matrix (const char *filename);

    //    py_bcsr_matrix_t get_regular_matrix () const;
    //    py_bcsr_matrix_t get_irregular_matrix () const;

    //    py_bcsr_matrix_t get_prepared_matrix () const;

    //    void
    //    prepare_matrix ();

    //    rhs_item_array_t &get_reg_acc_diag ();
    //    item_array_t &get_rhs ();
    //    item_array_t &get_rhs_flux ();
    //    item_array_t &get_sol ();

    //    item_array_t &get_sec_rhs ();
    //    item_array_t &get_sec_sol ();
    //    item_array_t &get_ss_diag ();
    //    item_array_t &get_sp_diag ();

    //  };

    void
    py_export_jacobian_matrix ();


  } // namespace python
} // namespace blue_sky

#endif  // #ifdef BSPY_EXPORTING_PLUGIN

#endif  // #ifndef PY_BS_JACOBIAN_MATRIX_H_
