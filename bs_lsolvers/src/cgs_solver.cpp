/**
 * \file cgs_solver.cpp
 * \brief CGS linear solver impl
 * \author Salimgareeva E.M. (algorithm from "A fast Lanczos-type solver for nonsymmetric linear systems", P.Sonneveld)
 * \date 11.02.2009
 * */

#include "bs_kernel.h"
#include "bs_assert.h"

#include "cgs_solver.h"
#include "mv_functions.h"


#include BS_FORCE_PLUGIN_IMPORT ()
#include "save_seq_vector.h"
#include "strategies.h"
#include "matrix_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {

    //////////////////////////////////////////////////////////////////////////
    //  cgs_solver

    //! constructor
    template <class strat_t>
    cgs_solver<strat_t>::cgs_solver (bs_type_ctor_param /*param*/)
      : lsolver_iface<strat_t> ()
    {
      prop = BS_KERNEL.create_object ("prop");
      if (!prop)
        {
          bs_throw_exception ("Type (prop) not registered");
        }
      init_prop ();
      sp_p                 = BS_KERNEL.create_object (fp_array_t::bs_type ()); 
      sp_phat              = BS_KERNEL.create_object (fp_array_t::bs_type ()); 
      sp_v                 = BS_KERNEL.create_object (fp_array_t::bs_type ()); 
      sp_tmp               = BS_KERNEL.create_object (fp_array_t::bs_type ()); 
      sp_q                 = BS_KERNEL.create_object (fp_array_t::bs_type ()); 
      sp_u                 = BS_KERNEL.create_object (fp_array_t::bs_type ()); 
      sp_d                 = BS_KERNEL.create_object (fp_array_t::bs_type ()); 
      sp_dhat              = BS_KERNEL.create_object (fp_array_t::bs_type ()); 
      sp_r                 = BS_KERNEL.create_object (fp_array_t::bs_type ()); 
      sp_rtilde            = BS_KERNEL.create_object (fp_array_t::bs_type ()); 
      sp_r_old             = BS_KERNEL.create_object (fp_array_t::bs_type ()); 
    }

    //! copy constructor
    template <class strat_t>
    cgs_solver<strat_t>::cgs_solver(const cgs_solver &solver)
      : bs_refcounter (), lsolver_iface<strat_t> ()
    {
      if (&solver != this)
        *this = solver;
    }

    //! destructor
    template <class strat_t>
    cgs_solver<strat_t>::~cgs_solver ()
    {}

    //! set solver's properties
    template <class strat_t>
    void cgs_solver<strat_t>::set_prop (sp_prop_t prop_)
    {
      prop = prop_;

      init_prop ();
    }
    template <class strat_t> void
    cgs_solver<strat_t>::init_prop ()
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
    int cgs_solver<strat_t>::solve (sp_matrix_t matrix, sp_fp_array_t sp_rhs, sp_fp_array_t sp_sol)
    {
#ifdef _DEBUG
      BOSOUT (section::solvers, level::debug) << "CGS\n" << bs_end;
#endif
      BS_ERROR (sp_rhs->size (), "cgs_solve");
      BS_ERROR (sp_sol->size (), "cgs_solve");
      BS_ERROR (prop, "cgs_solve");

      fp_type_t rho_1, rho_2 = 1, alpha = 1, beta, sigma;
      int iter;
      const double epsmac = 1e-24;
      fp_type_t r_norm, b_norm, den_norm;
      //fp_type *x = solution;
      fp_type_t *rhs = &(*sp_rhs)[0];
      fp_type_t *sol = &(*sp_sol)[0];

      const fp_type_t one = 1.0;
      //OMP_TIME_MEASURE_START (cgs_solve_timer);

      fp_type_t tol = prop->get_f (tol_idx);
      tol *= tol;
      //resid = prop->get_residuals ();
      //convergence_rate = prop->get_convergence_rate ();

      int max_iter  = prop->get_i (max_iters_idx);
      i_type_t n    = matrix->get_n_rows () * matrix->get_n_block_size ();
      BS_ASSERT (n == (i_type_t)sp_sol->size ());

      fp_type_t *p               = &(*sp_p)[0];
      fp_type_t *phat            = &(*sp_phat)[0];
      fp_type_t *v               = &(*sp_v)[0];
      fp_type_t *tmp             = &(*sp_tmp)[0];
      fp_type_t *q               = &(*sp_q)[0];
      fp_type_t *u               = &(*sp_u)[0];
      fp_type_t *d               = &(*sp_d)[0];
      fp_type_t *dhat            = &(*sp_dhat)[0];
      fp_type_t *r               = &(*sp_r)[0];
      fp_type_t *rtilde          = &(*sp_rtilde)[0];
      //fp_type_t *r_old           = &(*sp_r_old)[0];

      prop->set_b (success_idx, false);

      // solution = {0}
      //assign (solution, n, 0);
      memset (sol, 0, sizeof (fp_type_t) * n);
      //solution.assign (n, 0);

      sp_p->resize (n);
      sp_phat->resize (n);
      sp_v->resize (n);
      sp_tmp->resize (n);
      sp_q->resize (n);
      sp_u->resize (n);
      sp_d->resize (n);
      sp_dhat->resize (n);
      sp_r->resize (n);
      sp_rtilde->resize (n);
      sp_r_old->resize (n);

      // r = {0}
      //r.assign (n, 0);
      memset (r, 0, sizeof (fp_type_t) * n);
      memset (tmp, 0, sizeof (fp_type_t) * n);
      memset (p, 0, sizeof (fp_type_t) * n);
      memset (v, 0, sizeof (fp_type_t) * n);
      memset (q, 0, sizeof (fp_type_t) * n);

      
       // TODO:paste
      //tmp.assign (n, 0);
      //p.assign (n, 0);
      //v.assign (n, 0);
      //q.assign (n, 0);

      // p0 = u0 = r0;
      //u.assign (r.begin (), r.end ());
      memcpy (u, r, sizeof (fp_type_t) * n);
      // TODO:end

      // r = Ax0 - b
      matrix->calc_lin_comb (-1.0, 1.0, sp_sol, sp_rhs, sp_r);
      //rtilde.assign (r.begin (), r.end ());
      memcpy (rtilde, r, sizeof (fp_type_t) * n);

      //tools::save_seq_vector (tools::string_formater ("1_well_bhp.%s.txt", it->first).str).save (it->second);

      r_norm = mv_vector_inner_product_n (r, r, n);

      if (r_norm <= tol) // initial guess quite good
        return 0;

      rho_1 = r_norm;
      b_norm = sqrt (mv_vector_inner_product_n (rhs, rhs, n));

      // TODO:delete
      //p.assign      (r.begin (), r.end ());
      //rtilde.assign (r.begin (), r.end ());
      //v.assign      (n, 0);
      // TODO:end

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

      // set up initial norm and convergense factor
      //prop->set_relative_factor (den_norm);

      // main loop
      for (iter = 0; iter < max_iter; ++iter)
        {
          //printf ("CGS iteration: %d, resid = %le\n", iter, r_norm);
          //fflush (stdout);
          // TODO: paste
          if (iter)
            {
              //rho_1 = (r,rtilde)
              rho_1 = mv_vector_inner_product_n (r, rtilde, n); //in first iter equals to r_norm
              if (rho_1 == 0) // failure
                {
                  if (den_norm > epsmac)
                    prop->set_f (final_res_idx, r_norm / den_norm);
                  else
                    prop->set_f (final_res_idx, r_norm);
                  bs_throw_exception ("CGS: Failure - rho_1 == 0");
                }
            }

           beta = rho_1/rho_2; // beta = rho_n/rho_n-1
           rho_2 = rho_1;

           // u = r + beta*q
           sum_vector_n (r, one, q, beta, u, n);
           // tmp = q+beta*p_old
           sum_vector_n (q, one, p, beta, tmp, n);
           // p_new = u + beta*tmp
           sum_vector_n (u, one, tmp, beta, p, n);

           //temp_p.assign (p.begin (), p.end ());
           if (prec)
             {
               if (prec->solve_prec (matrix, sp_p, sp_phat))
                 {
                   bs_throw_exception ("CGS: Preconditioner failed");
                 }
             }
           else // no precondition (preconditioner=identity_matrix)
             {
               memcpy (phat, p, sizeof (fp_type_t) * n);
               //phat.assign (p.begin (), p.end ());
             }

          // v = A * phat = A * p, if no precondition;
          //v.assign (n, 0);
          memset (v, 0, sizeof (fp_type_t) * n);

          matrix->matrix_vector_product (sp_phat, sp_v);
          // sigma = (v,rtilde)
          sigma = mv_vector_inner_product_n (rtilde, v, n);

          if (sigma > epsmac || sigma < -epsmac)
          // alpha = rho_1/sigma
            alpha = rho_1 / sigma;
          else // failure
            {
              if (den_norm > epsmac)
                prop->set_f (final_res_idx, r_norm / den_norm);
              else
                prop->set_f (final_res_idx, r_norm);
              bs_throw_exception ("CGS: Failure - sigma == 0");
            }

          // q = u - alpha*v
          sum_vector_n (u, one, v, -alpha, q, n);
          // d = u + q
          sum_vector_n (u, one, q, one, d, n);

          // dhat = M^(-1) * d;
          //temp_d.assign (d.begin (), d.end ());
          if (prec)
            {
              if(prec->solve_prec (matrix, sp_d, sp_dhat))
                {
                  bs_throw_exception ("CGS: Preconditioner failed");
                }
            }
          else // no precondition (preconditioner=identity_matrix)
            {
              //dhat.assign (d.begin (), d.end ());
              memcpy (dhat, d, sizeof (fp_type_t) * n);
            }

          //tmp.assign (n, 0);
          memset (tmp, 0, sizeof (fp_type_t) * n);
          // tmp = A*d
          matrix->matrix_vector_product (sp_dhat, sp_tmp);

          // r = r - alpha*tmp
          sum_vector_n (r, one, tmp, -alpha, r, n);
          // x = x + alpha*dhat
          sum_vector_n (sol, one, dhat, alpha, sol, n);

          r_norm = mv_vector_inner_product_n (r, r, n);


          if (r_norm <= tol) // && check_resid_for_matbalance (n_rows, nb, r, matb_tol))
            break;
     }

     //tools::save_seq_vector ("solution.txt").save(solution);

     //TODO: end
     prop->set_i (iters_idx, iter + 1);
     prop->set_b (success_idx, true);

      /*
      //additional checking convergence
      mv_calc_lin_comb (matrix, -1.0, 1.0, solution, rhs, r);
      r_norm = mv_vector_inner_product (r, r, n);
      */
      if (den_norm > epsmac)
        prop->set_f (final_res_idx, r_norm / den_norm);
      else
        prop->set_f (final_res_idx, r_norm);

      //printf ("CGS OK! iters = %d, resid = %le\n", lprop->iters, lprop->final_resid);
      //OMP_TIME_MEASURE_END (bicgstab_solve_timer);

      return 0;
    }

    template <class strat_t>
    int cgs_solver<strat_t>::solve_prec (sp_matrix_t matrix, sp_fp_array_t rhs, sp_fp_array_t sol)
    {
      return solve (matrix, rhs, sol);
    }

    /**
    * @brief setup for CGS
    *
    * @param matrix -- input matrix
    *
    * @return 0 if success
    */
    template <class strat_t> int
    cgs_solver<strat_t>::setup (sp_matrix_t matrix)
    {
      if (!matrix)
        {
          bs_throw_exception ("CGS: Passed matrix is null");
        }

      BS_ASSERT (prop);
      if (prec)
        {
          return prec->setup (matrix);
        }

      return 0;
    }

    //////////////////////////////////////////////////////////////////////////
    BLUE_SKY_TYPE_STD_CREATE_T_DEF(cgs_solver, (class));
    BLUE_SKY_TYPE_STD_COPY_T_DEF(cgs_solver, (class));

    BLUE_SKY_TYPE_IMPL_T_EXT(1, (cgs_solver<base_strategy_fif>) , 1, (lsolver_iface<base_strategy_fif>), "cgs_solver_fif", "CGS linear solver", "CGS linear solver", false);
    BLUE_SKY_TYPE_IMPL_T_EXT(1, (cgs_solver<base_strategy_did>) , 1, (lsolver_iface<base_strategy_did>), "cgs_solver_did", "CGS linear solver", "CGS linear solver", false);
    BLUE_SKY_TYPE_IMPL_T_EXT(1, (cgs_solver<base_strategy_dif>) , 1, (lsolver_iface<base_strategy_dif>), "cgs_solver_dif", "CGS linear solver", "CGS linear solver", false);

    BLUE_SKY_TYPE_IMPL_T_EXT(1, (cgs_solver<base_strategy_flf>) , 1, (lsolver_iface<base_strategy_flf>), "cgs_solver_flf", "CGS linear solver", "CGS linear solver", false);
    BLUE_SKY_TYPE_IMPL_T_EXT(1, (cgs_solver<base_strategy_dld>) , 1, (lsolver_iface<base_strategy_dld>), "cgs_solver_dld", "CGS linear solver", "CGS linear solver", false);
    BLUE_SKY_TYPE_IMPL_T_EXT(1, (cgs_solver<base_strategy_dlf>) , 1, (lsolver_iface<base_strategy_dlf>), "cgs_solver_dlf", "CGS linear solver", "CGS linear solver", false);
  }
