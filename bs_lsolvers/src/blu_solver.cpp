/**
 * @file blu_solver.cpp
 * @brief block solver for linear systems using LU decomposition
 * @date 2009-12-09
 */

#include "bs_lsolvers_stdafx.h"

#include "blu_solver.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "strategies.h"
#include "matrix_iface.h"
#include "matrix_macroses.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {
  /*!
   * \brief constructor
   */

  blu_solver::blu_solver (bs_type_ctor_param /*param*/)
  {
      prop = BS_KERNEL.create_object ("prop");
      if (!prop)
        {
          bs_throw_exception ("Type (prop) not registered");
        }
      init_prop ();
  }


  blu_solver::blu_solver(const blu_solver& solver)
      : bs_refcounter (solver) //, lsolver_iface  ()
  {
    if (&solver != this)
      *this = solver;
  }

  //! set solver's properties

  void blu_solver::set_prop (sp_prop_t prop_)
  {
    prop = prop_;

    init_prop ();
  }

  void
  blu_solver::init_prop ()
    {
      prop->add_property_i (true, L"block_size",
                            L"Size of calculation block");
    }


  int blu_solver::solve (sp_matrix_t matrix, spv_double sp_rhs, spv_double sp_sol)
  {
    BS_ASSERT (matrix);
    BS_ASSERT (sp_rhs->size ());
    BS_ASSERT (sp_sol->size ());
    BS_ASSERT (sp_rhs->size () == sp_sol->size ()) (sp_rhs->size ()) (sp_sol->size ());

    t_long n;
    t_long nb;

    sp_dens_matrix_t lu;

    if (dynamic_cast<dens_matrix_iface_t *> (matrix.lock ()))
      {
        lu = matrix;
        BS_ASSERT (lu);
      }
    else
      {
        bs_throw_exception ("LU solve: Passed matrix is not dense");
      }

    spv_float spvalues = lu->get_values ();
    t_float *values = &(*spvalues)[0];
    t_double *sol = &(*sp_sol)[0];
    t_double *rhs = &(*sp_rhs)[0];


    n         = lu->get_n_rows ();
    nb        = lu->get_calc_block_size ();
    if (nb < 2)
      nb = n;

    BS_ASSERT (n == (t_long)sp_sol->size ());

    memcpy (sol, rhs, n * sizeof (t_double));

    return impl.lu_find_L_roots (n, values, sol, nb)
           || impl.lu_find_U_roots (n, values, sol, nb);
  }


  int blu_solver::solve_prec (sp_matrix_t matrix, spv_double rhs, spv_double sol)
  {
     return solve (matrix, rhs, sol);
  }

  /**
   * \brief setup preconditioner (merge matrices if needed)
   *
   * \param matrix Various matrix
   * \return 0 if success
   */
  int
  blu_solver::setup (sp_matrix_t matrix)
  {
    BS_ASSERT (matrix);

    t_long n;
    t_long nb;

    sp_dens_matrix_t lu;
    if (dynamic_cast<dens_matrix_iface_t *> (matrix.lock ()))
      {
        lu = matrix;
        BS_ASSERT (lu);
      }
    else
      {
        bs_throw_exception ("LU setup: Passed matrix is not dense");
      }

    spv_float values = lu->get_values ();

    n         = lu->get_n_rows ();
    nb        = lu->get_calc_block_size ();
    if (nb < 2)
      nb = n;

    return impl.lu_decomposition (n, &(*values)[0], nb);
  }


  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE (blu_solver);
  BLUE_SKY_TYPE_STD_COPY (blu_solver);

  BLUE_SKY_TYPE_IMPL (blu_solver, lsolver_iface, "blu_solver", "BLU solver", "Block LU solver for dense matricies");
} // namespace blue_sky

