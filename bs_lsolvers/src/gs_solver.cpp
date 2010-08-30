/** 
 * @file gs_solver.cpp
 * @brief Gause - Zeidel iterative solver 
 * @date 2009-12-16
 */

#include "bs_kernel.h"
#include "bs_assert.h"

#include "gs_solver.h"
#include "mv_functions.h"


#include BS_FORCE_PLUGIN_IMPORT ()
#include "bos_report.h"
#include "save_seq_vector.h"
#include "strategies.h"
#include "matrix_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {

    //////////////////////////////////////////////////////////////////////////
    //  gs_solver

    //! constructor
    template <class strat_t>
    gs_solver<strat_t>::gs_solver (bs_type_ctor_param /*param*/)
      : amg_smoother_iface<strat_t> ()
    {
      prop = BS_KERNEL.create_object ("prop");
      if (!prop)
        {
          bs_throw_exception ("Type (prop) not registered");
        }
      init_prop ();
      sp_r = BS_KERNEL.create_object (fp_array_t::bs_type ());
    }

    //! copy constructor
    template <class strat_t>
    gs_solver<strat_t>::gs_solver(const gs_solver &solver)
      : bs_refcounter (), amg_smoother_iface<strat_t> ()
    {
      if (&solver != this)
        *this = solver;
    }

    //! destructor
    template <class strat_t>
    gs_solver<strat_t>::~gs_solver ()
    {}

    //! set solver's properties
    template <class strat_t>
    void gs_solver<strat_t>::set_prop (sp_prop_t prop_)
    {
      prop = prop_;

      init_prop ();
    }

    template <class strat_t> void
    gs_solver<strat_t>::init_prop ()
      {
        tol_idx = prop->get_index_f (std::string ("tolerance"));
        if (tol_idx < 0)
          tol_idx = prop->add_property_f (1.0e-4, std::string ("tolerance"), std::string ("Target tolerance for linear solver"));

        final_res_idx = prop->get_index_f (std::string ("final_residual"));
        if (final_res_idx < 0)
          final_res_idx = prop->add_property_f (1, std::string ("final_residual"), std::string ("Solution residual"));

        max_iters_idx = prop->get_index_i (std::string ("maxiters"));
        if (max_iters_idx < 0)
          max_iters_idx = prop->add_property_i (20, std::string ("maxiters"), std::string ("Maximum allowed number of iterations"));

        iters_idx = prop->get_index_i (std::string ("iters"));
        if (iters_idx < 0)
          iters_idx = prop->add_property_i (0, std::string ("iters"), std::string ("Total number of used solver iterations"));

        success_idx = prop->get_index_b (std::string ("is_success"));
        if (success_idx < 0)
          success_idx = prop->add_property_b (false, std::string ("is_success"), std::string ("True if solver successfully convergent"));

        invers_idx = prop->get_index_b (std::string ("inverse"));
        if (invers_idx < 0)
          invers_idx = prop->add_property_b (false, std::string ("inverse"), std::string ("If True use inverse main loop"));

        cf_type_idx = prop->get_index_i (std::string ("cf_type"));
        if (cf_type_idx < 0)
          cf_type_idx = prop->add_property_i (0, std::string ("cf_type"), std::string ("0 -- use all points, -1 -- only negative points, "));
        
        if (tol_idx < 0
            || final_res_idx < 0
            || max_iters_idx < 0
            || iters_idx < 0
            || success_idx < 0)
          {
            bs_throw_exception ("Can not regidter some properties");
          }
      }
    template <class strat_t>
    int gs_solver<strat_t>::smooth (sp_bcsr_t matrix, sp_i_array_t cf_markers, const i_type_t iter_number, 
                                    sp_fp_array_t sp_rhs, sp_fp_array_t sp_sol)
      {
        BS_ASSERT (matrix);
        BS_ASSERT (sp_rhs->size ());
        BS_ASSERT (sp_sol->size ());
        BS_ASSERT (sp_rhs->size () == sp_sol->size ()) (sp_rhs->size ()) (sp_sol->size ());
        BS_ASSERT (matrix->get_n_block_size () >= 1) (matrix->get_n_block_size ());

        i_type_t i, j1, j2, j, cl;
        fp_type_t d, b;
        static const double eps = 1.0e-24;
        bool inv = prop->get_b (invers_idx);
        int  cf  = prop->get_i (cf_type_idx);

        i_type_t n                      = matrix->get_n_rows ();
        sp_fp_storage_array_t sp_data         = matrix->get_values ();
        sp_i_array_t          sp_rows         = matrix->get_rows_ptr ();
        sp_i_array_t          sp_cols         = matrix->get_cols_ind ();

        fp_storage_type_t     *data           = &(*sp_data)[0];
        i_type_t              *rows           = &(*sp_rows)[0];
        i_type_t              *cols           = &(*sp_cols)[0];

        fp_type_t *rhs = &(*sp_rhs)[0];
        fp_type_t *sol = &(*sp_sol)[0];
        i_type_t *cf_m = 0;
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
          for (i_type_t ii = 0; ii < iter_number; ++ii)
            for (i = 0; i < n; ++i)
              {
                if ((cf * cf_m[i]) > 0)
                  {
                    GAUSS_LOOP
                  }
              }
        if ((inv) && (cf!=0))
          for (i_type_t ii = 0; ii < iter_number; ++ii)
            for (i = n-1; i >= 0; --i)
              {
                if ((cf * cf_m[i]) > 0)
                  {
                    GAUSS_LOOP
                  }
              }
        if ((!inv) && (cf==0))
          for (i_type_t ii = 0; ii < iter_number; ++ii)
            for (i = 0; i < n; ++i)
              {
                GAUSS_LOOP
              }
        if ((inv) && (cf==0))
          for (i_type_t ii = 0; ii < iter_number; ++ii)
            for (i = n-1; i >= 0; --i)
              {
                GAUSS_LOOP
              }
        return 0;
      }

    template <class strat_t>
    int gs_solver<strat_t>::solve (sp_matrix_t matrix, sp_fp_array_t sp_rhs, sp_fp_array_t sp_sol)
    {
      BS_ERROR (matrix, "gs_solve");
      BS_ERROR (sp_rhs->size (), "gs_solve");
      BS_ERROR (sp_sol->size (), "gs_solve");
      BS_ERROR (prop, "gs_solve");

      int iter;
      const double epsmac = 1e-24;
      fp_type_t r_norm, b_norm, den_norm;

      fp_type_t *rhs = &(*sp_rhs)[0];
      //fp_type_t *sol = &(*sp_sol)[0];

      sp_bcsr_t bcsr;
      if (!dynamic_cast<bcsr_t *> (matrix.lock ()))
        {
          bcsr = matrix;
          BS_ASSERT (bcsr);
        }

      i_type_t n = matrix->get_n_rows () * matrix->get_n_block_size ();
      sp_i_array_t flags;

      fp_type_t tol = prop->get_f (tol_idx);
      tol *= tol;

      int max_iter  = prop->get_i (max_iters_idx);
      prop->reset_i (cf_type_idx);

      prop->set_b (success_idx, false);

      matrix->init_vector (sp_r);
      fp_type_t *r = &(*sp_r)[0];

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

    template <class strat_t>
    int gs_solver<strat_t>::solve_prec(sp_matrix_t matrix, sp_fp_array_t rhs, sp_fp_array_t solution)
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
    template <class strat_t> int
    gs_solver<strat_t>::setup (sp_matrix_t matrix)
    {
      if (!matrix)
        {
          bs_throw_exception ("gs: Passed matrix is null");
        }

      return 0;
    }

    //////////////////////////////////////////////////////////////////////////
    BLUE_SKY_TYPE_STD_CREATE_T_DEF(gs_solver, (class));
    BLUE_SKY_TYPE_STD_COPY_T_DEF(gs_solver, (class));

    BLUE_SKY_TYPE_IMPL_T_EXT(1, (gs_solver<base_strategy_fif>) , 1, (lsolver_iface<base_strategy_fif>), "gs_solver_base_fif", "gs linear solver", "gs linear solver", false);
    BLUE_SKY_TYPE_IMPL_T_EXT(1, (gs_solver<base_strategy_did>) , 1, (lsolver_iface<base_strategy_did>), "gs_solver_base_did", "gs linear solver", "gs linear solver", false);
    BLUE_SKY_TYPE_IMPL_T_EXT(1, (gs_solver<base_strategy_dif>) , 1, (lsolver_iface<base_strategy_dif>), "gs_solver_base_dif", "gs linear solver", "gs linear solver", false);

    BLUE_SKY_TYPE_IMPL_T_EXT(1, (gs_solver<base_strategy_flf>) , 1, (lsolver_iface<base_strategy_flf>), "gs_solver_base_flf", "gs linear solver", "gs linear solver", false);
    BLUE_SKY_TYPE_IMPL_T_EXT(1, (gs_solver<base_strategy_dld>) , 1, (lsolver_iface<base_strategy_dld>), "gs_solver_base_dld", "gs linear solver", "gs linear solver", false);
    BLUE_SKY_TYPE_IMPL_T_EXT(1, (gs_solver<base_strategy_dlf>) , 1, (lsolver_iface<base_strategy_dlf>), "gs_solver_base_dlf", "gs linear solver", "gs linear solver", false);
  }

