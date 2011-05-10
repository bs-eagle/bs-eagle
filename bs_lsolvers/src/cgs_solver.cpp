/**
 * \file cgs_solver.cpp
 * \brief CGS linear solver impl
 * \author Salimgareeva E.M. (algorithm from "A fast Lanczos-type solver for nonsymmetric linear systems", P.Sonneveld)
 * \date 11.02.2009
 * */

#include "lsolvers_stdafx.h"

#include "cgs_solver.h"
#include "mv_functions.h"


#include BS_FORCE_PLUGIN_IMPORT ()
#include "matrix_iface.h"
#include "bos_report.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {

    //////////////////////////////////////////////////////////////////////////
    //  cgs_solver

    //! constructor
    
    cgs_solver::cgs_solver (bs_type_ctor_param /*param*/)
      : lsolver_iface ()
    {
      prop = BS_KERNEL.create_object ("prop");
      if (!prop)
        {
          bs_throw_exception ("Type (prop) not registered");
        }
      init_prop ();
      sp_p                 = BS_KERNEL.create_object (v_double::bs_type ()); 
      sp_phat              = BS_KERNEL.create_object (v_double::bs_type ()); 
      sp_v                 = BS_KERNEL.create_object (v_double::bs_type ()); 
      sp_tmp               = BS_KERNEL.create_object (v_double::bs_type ()); 
      sp_q                 = BS_KERNEL.create_object (v_double::bs_type ()); 
      sp_u                 = BS_KERNEL.create_object (v_double::bs_type ()); 
      sp_d                 = BS_KERNEL.create_object (v_double::bs_type ()); 
      sp_dhat              = BS_KERNEL.create_object (v_double::bs_type ()); 
      sp_r                 = BS_KERNEL.create_object (v_double::bs_type ()); 
      sp_rtilde            = BS_KERNEL.create_object (v_double::bs_type ()); 
      sp_r_old             = BS_KERNEL.create_object (v_double::bs_type ()); 
    }

    //! copy constructor
    
    cgs_solver::cgs_solver(const cgs_solver &solver)
      : bs_refcounter (), lsolver_iface ()
    {
      if (&solver != this)
        *this = solver;
    }

    //! destructor
    
    cgs_solver::~cgs_solver ()
    {}

    //! set solver's properties
    
    void cgs_solver::set_prop (sp_prop_t prop_)
    {
      prop = prop_;

      init_prop ();
    }
     void
    cgs_solver::init_prop ()
      {
        prop->add_property_f (1.0e-4, tol_idx, 
                              std::string ("Target tolerance for linear solver"));

        prop->add_property_f (1, final_res_idx, std::string ("Solution residual"));

        prop->add_property_i (20, max_iters_idx, 
                              std::string ("Maximum allowed number of iterations"));

        prop->add_property_i (0, iters_idx, 
                              std::string ("Total number of used solver iterations"));

        prop->add_property_b (false, success_idx, 
                              std::string ("True if solver successfully convergent"));
      }

    
    int cgs_solver::solve (sp_matrix_t matrix, spv_double sp_rhs, spv_double sp_sol)
    {
#ifdef _DEBUG
      BOSOUT (section::solvers, level::debug) << "CGS\n" << bs_end;
#endif
      BS_ERROR (sp_rhs->size (), "cgs_solve");
      BS_ERROR (sp_sol->size (), "cgs_solve");
      BS_ERROR (prop, "cgs_solve");

      t_double rho_1, rho_2 = 1, alpha = 1, beta, sigma;
      int iter;
      const double epsmac = 1e-24;
      t_double r_norm, b_norm, den_norm;
      //fp_type *x = solution;
      t_double *rhs = &(*sp_rhs)[0];
      t_double *sol = &(*sp_sol)[0];

      const t_double one = 1.0;
      //OMP_TIME_MEASURE_START (cgs_solve_timer);

      t_double tol = prop->get_f (tol_idx);
      tol *= tol;
      //resid = prop->get_residuals ();
      //convergence_rate = prop->get_convergence_rate ();

      int max_iter  = prop->get_i (max_iters_idx);
      t_long n    = matrix->get_n_rows () * matrix->get_n_block_size ();
      BS_ASSERT (n == (t_long)sp_sol->size ());
      
      t_double *p               = &(*sp_p)[0];
      t_double *phat            = &(*sp_phat)[0];
      t_double *v               = &(*sp_v)[0];
      t_double *tmp             = &(*sp_tmp)[0];
      t_double *q               = &(*sp_q)[0];
      t_double *u               = &(*sp_u)[0];
      t_double *d               = &(*sp_d)[0];
      t_double *dhat            = &(*sp_dhat)[0];
      t_double *r               = &(*sp_r)[0];
      t_double *rtilde          = &(*sp_rtilde)[0];
      //t_double *r_old           = &(*sp_r_old)[0];

      prop->set_b (success_idx, false);

      // solution = {0}
      //assign (solution, n, 0);
      memset (sol, 0, sizeof (t_double) * n);
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
      memset (r, 0, sizeof (t_double) * n);
      memset (tmp, 0, sizeof (t_double) * n);
      memset (p, 0, sizeof (t_double) * n);
      memset (v, 0, sizeof (t_double) * n);
      memset (q, 0, sizeof (t_double) * n);

      
       // TODO:paste
      //tmp.assign (n, 0);
      //p.assign (n, 0);
      //v.assign (n, 0);
      //q.assign (n, 0);

      // p0 = u0 = r0;
      //u.assign (r.begin (), r.end ());
      memcpy (u, r, sizeof (t_double) * n);
      // TODO:end

      // r = Ax0 - b
      matrix->calc_lin_comb (-1.0, 1.0, sp_sol, sp_rhs, sp_r);
      //rtilde.assign (r.begin (), r.end ());
      memcpy (rtilde, r, sizeof (t_double) * n);

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
               memcpy (phat, p, sizeof (t_double) * n);
               //phat.assign (p.begin (), p.end ());
             }

          // v = A * phat = A * p, if no precondition;
          //v.assign (n, 0);
          memset (v, 0, sizeof (t_double) * n);

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
              memcpy (dhat, d, sizeof (t_double) * n);
            }

          //tmp.assign (n, 0);
          memset (tmp, 0, sizeof (t_double) * n);
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

    
    int cgs_solver::solve_prec (sp_matrix_t matrix, spv_double rhs, spv_double sol)
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
     int
    cgs_solver::setup (sp_matrix_t matrix)
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
    BLUE_SKY_TYPE_STD_CREATE (cgs_solver);
    BLUE_SKY_TYPE_STD_COPY (cgs_solver);

    BLUE_SKY_TYPE_IMPL (cgs_solver, lsolver_iface, "cgs_solver", "CGS linear solver", "CGS linear solver");
  }
