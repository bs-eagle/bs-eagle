/**
 * @file gs_solver.cpp
 * @brief Gause - Zeidel iterative solver
 * @date 2009-12-16
 */

#include "bs_lsolvers_stdafx.h"

#include "gs_solver.h"
#include "mv_functions.h"


#include BS_FORCE_PLUGIN_IMPORT ()
#include "bos_report.h"
#include "matrix_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {

    //////////////////////////////////////////////////////////////////////////
    //  gs_solver

    //! constructor

    gs_solver::gs_solver (bs_type_ctor_param /*param*/)
      : amg_smoother_iface ()
    {
      prop = BS_KERNEL.create_object ("prop");
      if (!prop)
        {
          bs_throw_exception ("Type (prop) not registered");
        }
      init_prop ();
      sp_r = BS_KERNEL.create_object (v_double::bs_type ());
    }

    //! copy constructor

    gs_solver::gs_solver(const gs_solver &solver)
      : bs_refcounter (), amg_smoother_iface ()
    {
      if (&solver != this)
        *this = solver;
    }

    //! destructor

    gs_solver::~gs_solver ()
    {}

    //! set solver's properties

    void gs_solver::set_prop (sp_prop_t prop_)
    {
      prop = prop_;

      init_prop ();
    }

     void
    gs_solver::init_prop ()
      {
        prop->add_property_f (1.0e-4, tol_idx,
                              L"Target tolerance for linear solver");
        prop->add_property_f (1, final_res_idx,
                              L"Solution residual");
        prop->add_property_i (20, max_iters_idx,
                              L"Maximum allowed number of iterations");
        prop->add_property_i (0, iters_idx,
                              L"Total number of used solver iterations");
        prop->add_property_b (false, success_idx,
                              L"True if solver successfully convergent");
        prop->add_property_b (false, invers_idx,
                              L"If True use inverse main loop");
        prop->add_property_i (0, cf_type_idx,
                              L"0 -- use all points, -1 -- only negative points, ");
      }

    int gs_solver::smooth (sp_bcsr_t matrix, spv_long cf_markers, const t_long iter_number,
                                    spv_double sp_rhs, spv_double sp_sol)
      {
        BS_ASSERT (matrix);
        BS_ASSERT (sp_rhs->size ());
        BS_ASSERT (sp_sol->size ());
        BS_ASSERT (sp_rhs->size () == sp_sol->size ()) (sp_rhs->size ()) (sp_sol->size ());
        BS_ASSERT (matrix->get_n_block_size () >= 1) (matrix->get_n_block_size ());

        t_long i, j1, j2, j, cl;
        t_double d, b;
        static const double eps = 1.0e-24;
        bool inv = prop->get_b (invers_idx);
        int  cf  = prop->get_i (cf_type_idx);

        t_long n                      = matrix->get_n_rows ();
        spv_float sp_data         = matrix->get_values ();
        spv_long          sp_rows         = matrix->get_rows_ptr ();
        spv_long          sp_cols         = matrix->get_cols_ind ();

        t_float     *data           = &(*sp_data)[0];
        t_long              *rows           = &(*sp_rows)[0];
        t_long              *cols           = &(*sp_cols)[0];

        t_double *rhs = &(*sp_rhs)[0];
        t_double *sol = &(*sp_sol)[0];
        t_long *cf_m = 0;
        if (cf)
          cf_m                  = &(*cf_markers)[0];



#define GAUSS_LOOP                              \
            j1 = rows[i];                           \
            j2 = rows[i + 1];                       \
            d = data[j1];                           \
            if (d > eps || d < -eps)                \
              {                                     \
                b = rhs[i];                         \
                for (j = j1 + 1; j < j2; ++j)       \
                  {                                 \
                    cl = cols[j];                   \
                    b -= sol[cl] * data[j];         \
                  }                                 \
                sol[i] = b / d;                     \
              }

        // TODO: OPENMP BUG
        // main loop
#ifdef OTHER_NON_IMPORTANT_PARALLEL
#pragma omp parallel for private (j1, j2, j, b, d, cl)
#endif //OTHER_NON_IMPORTANT_PARALLEL

        if ((!inv) && (cf!=0))
          for (t_long ii = 0; ii < iter_number; ++ii)
            for (i = 0; i < n; ++i)
              {
                if ((cf * cf_m[i]) >= 0)
                  {
                    GAUSS_LOOP
                  }
              }
        if ((inv) && (cf!=0))
          for (t_long ii = 0; ii < iter_number; ++ii)
            for (i = n-1; i >= 0; --i)
              {
                if ((cf * cf_m[i]) >= 0)
                  {
                    GAUSS_LOOP
                  }
              }
        if ((!inv) && (cf==0))
          for (t_long ii = 0; ii < iter_number; ++ii)
            for (i = 0; i < n; ++i)
              {
                GAUSS_LOOP
              }
        if ((inv) && (cf==0))
          for (t_long ii = 0; ii < iter_number; ++ii)
            for (i = n-1; i >= 0; --i)
              {
                GAUSS_LOOP
              }
        return 0;
      }


    int gs_solver::solve (sp_matrix_t matrix, spv_double sp_rhs, spv_double sp_sol)
    {
      BS_ERROR (matrix, "gs_solve");
      BS_ERROR (sp_rhs->size (), "gs_solve");
      BS_ERROR (sp_sol->size (), "gs_solve");
      BS_ERROR (prop, "gs_solve");

      int iter;
      const double epsmac = 1e-24;
      t_double r_norm, b_norm, den_norm;

      t_double *rhs = &(*sp_rhs)[0];
      //t_double *sol = &(*sp_sol)[0];

      sp_bcsr_t bcsr;
      if (!dynamic_cast<bcsr_t *> (matrix.get ()))
        {
          bcsr = matrix;
          BS_ASSERT (bcsr);
        }

      t_long n = matrix->get_n_rows () * matrix->get_n_block_size ();
      spv_long flags;

      t_double tol = prop->get_f (tol_idx);
      tol *= tol;

      int max_iter  = prop->get_i (max_iters_idx);
      prop->reset_i (cf_type_idx);

      prop->set_b (success_idx, false);

      matrix->init_vector (sp_r);
      t_double *r = &(*sp_r)[0];

      matrix->calc_lin_comb (-1.0, 1.0, sp_sol, sp_rhs, sp_r);
      r_norm = mv_vector_inner_product_n (r, r, n);

      b_norm = sqrt (mv_vector_inner_product_n (rhs, rhs, n));

      if (b_norm > epsmac) // choose convergence criterion
        {
          // |r_i|/|b| <= eps if |b| > 0
          tol *= b_norm;
          den_norm = b_norm;
        }
      else // (r_norm > epsmac)
        {
          // |r_i|/|r0| <= eps if |b| = 0
          tol *= r_norm;
          den_norm = r_norm;
        }

      // main loop
      for (iter = 0; iter < max_iter; ++iter)
        {

          smooth (bcsr, flags, 1, sp_rhs, sp_sol);
          matrix->calc_lin_comb (-1.0, 1.0, sp_sol, sp_rhs, sp_r);
          r_norm = mv_vector_inner_product_n (r, r, n);
          if (r_norm <= tol) // initial guess quite good
            break;
        } // end of main loop

      if (iter < max_iter)
        {
          prop->set_i (iters_idx, iter + 1);
          prop->set_b (success_idx, true);
        }
      else
        {
          prop->set_i (iters_idx, iter + 1);
          prop->set_b (success_idx, false);
        }

      if (den_norm > epsmac)
        prop->set_f(final_res_idx, r_norm / den_norm);
      else
        prop->set_f (final_res_idx, r_norm);

      BOSOUT (section::solvers, level::low) << "r_norm = " << r_norm << " r_norm / den_norm = " << r_norm / den_norm << " iter = " << iter << bs_end;

      return 0;
    }


    int gs_solver::solve_prec(sp_matrix_t matrix, spv_double rhs, spv_double solution)
    {
      return solve (matrix, rhs, solution);
    }

    /**
    * @brief setup for CGS
    *
    * @param matrix -- input matrix
    *
    * @return 0 if success
    */
     int
    gs_solver::setup (sp_matrix_t matrix)
    {
      if (!matrix)
        {
          bs_throw_exception ("gs: Passed matrix is null");
        }

      return 0;
    }

    //////////////////////////////////////////////////////////////////////////
    BLUE_SKY_TYPE_STD_CREATE (gs_solver);
    BLUE_SKY_TYPE_STD_COPY (gs_solver);

    BLUE_SKY_TYPE_IMPL (gs_solver, lsolver_iface, "gs_solver", "gs linear solver", "gs linear solver");
  }

