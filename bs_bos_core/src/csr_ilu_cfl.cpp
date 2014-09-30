#if 0
/**
 *       \file  csr_ilu_cfl.cpp
 *      \brief  Builds matrix based on CFL for csr_ilu preconditioner
 *              and solves this matrix
 *     \author  Elmira Salimgareeva
 *       \date  11.05.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "csr_ilu_cfl.h"

namespace blue_sky
{
  /*!
  brief constructor
  */
  template <class strategy_t>
  csr_ilu_cfl_prec<strategy_t>::csr_ilu_cfl_prec (bs_type_ctor_param param /* = NULL */)
    : csr_ilu_prec<strategy_t>(param),
    sp_ilu_cfl_(BS_KERNEL.create_object(bcsr_matrix_t::bs_type()))
    {
    }

  template <class strategy_t>
  csr_ilu_cfl_prec<strategy_t>::csr_ilu_cfl_prec(const csr_ilu_cfl_prec& solver)
    : bs_refcounter (solver), 
    csr_ilu_prec <strategy_t> (solver)
  {
    if (&solver != this)
      *this = solver;
  }

  /*!
  * \brief destructor
  */
  template <class strategy_t>
  csr_ilu_cfl_prec<strategy_t>::~csr_ilu_cfl_prec ()
  {
  }

  template <class strategy_t> int
  csr_ilu_cfl_prec<strategy_t>::solve (matrix_t * /*matrix*/, rhs_item_array_t &rhs, item_array_t &sol)
  {
    return this->csr_ilu_prec_t::solve (sp_ilu_cfl_.lock (), rhs, sol);
  }

  template <class strategy_t> int
  csr_ilu_cfl_prec<strategy_t>::solve_prec (matrix_t * /*matrix*/, item_array_t &rhs, item_array_t &sol)
  {
    return this->csr_ilu_prec_t::solve_prec (sp_ilu_cfl_.lock (), rhs, sol);
  }


  template <class strategy_t>
  typename csr_ilu_cfl_prec<strategy_t>::sp_bcsr_matrix_t
  csr_ilu_cfl_prec<strategy_t>::create_merged_matrix_for_ilu (jmatrix_t *matrix, const rhs_item_array_t &cfl_vector)
  {
   if (!matrix)
      {
        throw bs_exception ("csr_ilu_cfl: create_merge_matrix_for_ilu", "Passed matrix is null");
      }


   BS_ASSERT (sp_ilu_cfl_);
   sp_bcsr_matrix_t locked_regular (matrix->get_regular_matrix ());
   sp_bcsr_matrix_t locked_irregular (matrix->get_irregular_matrix ());

   if (!locked_regular)
     throw bs_exception ("csr_ilu_cfl: create_merge_matrix_for_ilu", "Regular matrix is null");

   if (!locked_irregular)
     throw bs_exception ("csr_ilu_cfl: create_merge_matrix_for_ilu", "Irregular matrix is null");

   //----------------------------------------------------------------------------

   bcsr_matrix_t &regular_matrix   = *locked_regular;
   bcsr_matrix_t &irregular_matrix = *locked_irregular;
   const rhs_item_array_t &regular_acc_diag = matrix->get_regular_acc_diag ();

   irregular_matrix.init_diag_ind ();

   index_array_t &reg_rows      = regular_matrix.get_rows_ptr ();
   index_array_t &irreg_rows    = irregular_matrix.get_rows_ptr ();
   //---------------------------------------------------------------------------

   // if !condition then print message about error and values of variables
   BS_ASSERT (regular_matrix.n_rows == irregular_matrix.n_rows) (regular_matrix.n_rows) (irregular_matrix.n_rows);

   if (regular_matrix.n_rows != irregular_matrix.n_rows)
     throw bs_exception ("csr_ilu_cfl: create_merge_matrix_for_ilu", "Number of rows in regular matrix is not the same as in irregular");

   if (regular_matrix.n_block_size != irregular_matrix.n_block_size)
     throw bs_exception ("csr_ilu_cfl: create_merge_matrix_for_ilu", "Size of block in regular matrix is not the same as in irregular");

   //----------------------------------------------------------------------------


   bcsr_matrix_t &new_matrix = *sp_ilu_cfl_;

   if (new_matrix.init (1, 1, regular_matrix.n_block_size, 1))
     throw bs_exception ("csr_ilu_cfl: create_merge_matrix_for_ilu", "can't init merged_reg_irreg_matrix");


   //bool ok = false;
   if (reg_rows.empty () || irreg_rows.empty ())
     throw bs_exception ("csr_ilu_cfl: create_merge_matrix_for_ilu", "Reg rows or irreg rows is empty");


   if (new_matrix.merge (&regular_matrix, &irregular_matrix, 1, cfl_vector))
     throw bs_exception ("csr_ilu_cfl: create_merge_matrix_for_ilu", "can't merge regular and irregular matricies");

   if (new_matrix.init_diag_ind ())
     throw bs_exception ("csr_ilu_cfl: create_merge_matrix_for_ilu", "can't init merged reg and irreg matricies diag ind's");

   index_t n_block_size = regular_matrix.n_block_size; // size of block in matrix
   index_t b_sqr = n_block_size*n_block_size;

   rhs_item_array_t &new_vals = new_matrix.get_values ();

   index_array_t &new_rows = new_matrix.get_rows_ptr ();
   index_array_t &new_cols = new_matrix.get_cols_ind ();

   for (index_t i = 0, n_rows = regular_matrix.n_rows; i < n_rows; ++i)
     {
       index_t l = new_rows[i];
       BS_ASSERT (new_cols[l] == i) (new_cols[l]) (i);
       for (index_t k = 0; k < b_sqr; ++k)
         {
           new_vals[l * b_sqr + k] += regular_acc_diag[i * b_sqr + k];
         }
     }

    new_matrix.n_rows = regular_matrix.n_rows;
    new_matrix.n_cols = regular_matrix.n_cols;
    new_matrix.n_block_size = regular_matrix.n_block_size;

    return sp_ilu_cfl_;//locked_new_matrix;
  }


  template <class strategy_t> int
  csr_ilu_cfl_prec<strategy_t>::setup (matrix_t *matrix_)
  {
    BS_ASSERT (matrix_);
    if (!matrix_)
      throw bs_exception ("csr_ilu_cfl::setup", "Passed matrix is null");

    jacobian_matrix <strategy_t> *jmx = dynamic_cast <jacobian_matrix <strategy_t> *> (matrix_);
    BS_ASSERT (jmx) (bs::type_name (*matrix_));
    if (jmx)
      {
        sp_ilu_cfl_ = create_merged_matrix_for_ilu (jmx, jmx->get_cfl_vector ());
        return csr_ilu_prec_t::setup_internal (sp_ilu_cfl_);
      }
    else
      {
        return csr_ilu_prec_t::setup_internal (dynamic_cast <bcsr_matrix_t *> (matrix_));
      }
  }


//////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE_T_DEF(csr_ilu_cfl_prec, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF(csr_ilu_cfl_prec, (class));

  BLUE_SKY_TYPE_IMPL_T_EXT(1, (csr_ilu_cfl_prec<base_strategy_fi>) , 1, (linear_solver_base<base_strategy_fi>), "csr_ilu_cfl_prec_seq_fi", "CSR ILU Preconditioner", "CSR ILU Preconditioner", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (csr_ilu_cfl_prec<base_strategy_di>) , 1, (linear_solver_base<base_strategy_di>), "csr_ilu_cfl_prec_seq_di", "CSR ILU Preconditioner", "CSR ILU Preconditioner", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (csr_ilu_cfl_prec<base_strategy_mixi>) , 1, (linear_solver_base<base_strategy_mixi>), "csr_ilu_cfl_prec_seq_mixi", "CSR ILU Preconditioner", "CSR ILU Preconditioner", false);

//////////////////////////////////////////////////////////////////////////


}
#endif // #if 0
