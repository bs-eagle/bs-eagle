/**
 *       \file  two_stage_preconditioner.cpp
 *      \brief  Class declaration for two stage preconditioner
 *     \author  Borschuk Oleg
 *       \date  04.08.2006
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"
#include "two_stage_preconditioner.h"
#include "solve_helper.h"

namespace blue_sky
  {
  /**
   * \brief  'default' ctor for two_stage_preconditioner
   * \param  param Additional params for ctor
   * */
  two_stage_preconditioner::two_stage_preconditioner (bs_type_ctor_param param)
      : linear_solver_base (param)
  {
    //sp_link prec_2_link = bs_link::create (bs_node::create_node (), "prec_2");

    //bs_node::insert (prec_2_link, false);
  }

  /**
   * \brief  copy-ctor for two_stage_preconditioner
   * \param  prec Instance of two_stage_preconditioner to be copied
   * */
  two_stage_preconditioner::two_stage_preconditioner (const two_stage_preconditioner &prec)
      : bs_refcounter (prec), linear_solver_base (prec)
  {
    if (this != &prec)
      *this = prec;
  }

  //! destructor
  two_stage_preconditioner::~two_stage_preconditioner ()
  {
  }

  int
  two_stage_preconditioner::solve(matrix_t *matrix, rhs_item_array_t &rhs, item_array_t &sol)
  {
    return templ_solve (matrix, rhs, sol);
  }

  int
  two_stage_preconditioner::solve_prec(matrix_t *matrix, item_array_t &rhs, item_array_t &sol)
  {
    return templ_solve (matrix, rhs, sol);
  }


  //! solve preconditioner
  template <class rhs_t> int
  two_stage_preconditioner::templ_solve (matrix_t *matrix, rhs_t &rhs, item_array_t &sol)
  {
    BS_ASSERT (matrix);
    BS_ASSERT (rhs.size ()) (rhs.size ());
    BS_ASSERT (sol.size ()) (sol.size ());
    BS_ASSERT (rhs.size () == sol.size ()) (rhs.size ()) (sol.size ());
    BS_ASSERT (base_t::prec);
    BS_ASSERT (prec_2);

    index_t n = matrix->n_rows * matrix->n_block_size;
    BS_ASSERT ((size_t)n == rhs.size ()) (n) (rhs.size ());

    if (solve_helper (&(*base_t::prec), matrix, rhs, sol))
      {
        bs_throw_exception ("TWO_STAGE: PREC 1 failed");
      }

    if ((index_t)r_array.size () != n)
      {
        r_array.assign (n, 0);
        w_array.assign (n, 0);
      }
    BS_ASSERT (r_array.size () == w_array.size ()) (r_array.size ()) (w_array.size ());

    //r_array = rhs - sol - new rhs - step 5 in Olegs book
    matrix->calc_lin_comb (-1.0, 1.0, sol, rhs, r_array);

    // solve system with ILU - step 6
    if (prec_2->solve_prec (matrix, r_array, w_array))
      {
        bs_throw_exception ("TWO_STAGE: PREC 2 failed");
      }

    // Then we make a correction of our solution x = x + w - step 7
    index_t i = 0;
    index_t n2 = n - (n % 4);
    for (; i < n2; i+=4)
      {
        sol[i] += w_array[i];
        sol[i+1] += w_array[i+1];
        sol[i+2] += w_array[i+2];
        sol[i+3] += w_array[i+3];
      }

    for (; i < n; i++)
      {
        sol[i] += w_array[i];
      }

    return 0;
  }

  //! setup preconditioner
  int
  two_stage_preconditioner::setup (matrix_t *matrix)
  {
    if (!matrix)
      {
        bs_throw_exception ("TWO_STAGE: Passed matrix is null");
      }

    BS_ASSERT (base_t::prec);
    BS_ASSERT (prec_2);

    base_t::prec->setup (matrix);
    prec_2->setup (matrix);

    return 0;
  }

  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE (two_stage_preconditioner);
  BLUE_SKY_TYPE_STD_COPY (two_stage_preconditioner);
  BLUE_SKY_TYPE_IMPL (two_stage_preconditioner, linear_solver_base, "two_stage_prec", "two_stage_prec", "two_stage_prec");

  bool two_stage_prec_register_type (const blue_sky::plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, two_stage_preconditioner::bs_type ());

    return res;
  }



} // namespace blue_sky
