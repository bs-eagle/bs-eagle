/**
 * \file bicgstab_solver.cpp
 * \brief CGS linear solver impl
 * \author Salimgareeva E.M. (algorithm from "A fast Lanczos-type solver for nonsymmetric linear systems", P.Sonneveld)
 * \date 11.02.2009
 * */

#include "bicgstab_solver.h"
#include "mv_functions.h"


#include BS_FORCE_PLUGIN_IMPORT ()
#include "bos_report.h"
#include "matrix_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {

    //////////////////////////////////////////////////////////////////////////
    //  bicgstab_solver

    //! constructor
    
    bicgstab_solver::bicgstab_solver (bs_type_ctor_param /*param*/)
      : lsolver_iface ()
    {
      prop = BS_KERNEL.create_object ("prop");
      if (!prop)
        {
          bs_throw_exception ("Type (prop) not registered");
        }
      init_prop ();
      sp_p = BS_KERNEL.create_object (v_double::bs_type ());
      sp_phat = BS_KERNEL.create_object (v_double::bs_type ());
      sp_s = BS_KERNEL.create_object (v_double::bs_type ());
      sp_shat = BS_KERNEL.create_object (v_double::bs_type ());
      sp_t = BS_KERNEL.create_object (v_double::bs_type ());
      sp_v = BS_KERNEL.create_object (v_double::bs_type ());
      sp_r = BS_KERNEL.create_object (v_double::bs_type ());
      sp_rtilde = BS_KERNEL.create_object (v_double::bs_type ());

    }

    //! copy constructor
    
    bicgstab_solver::bicgstab_solver(const bicgstab_solver &solver)
      : bs_refcounter (), lsolver_iface ()
    {
      if (&solver != this)
        *this = solver;
    }

    //! destructor
    
    bicgstab_solver::~bicgstab_solver ()
    {}

    //! set solver's properties
    
    void bicgstab_solver::set_prop(sp_prop_t prop_)
    {
      prop = prop_;

      init_prop ();
    }

     void
    bicgstab_solver::init_prop ()
      {
        prop->add_property_f (1.0e-4, tol_idx, 
                              L"Target tolerance for linear solver");

        prop->add_property_f (1, final_res_idx, L"Solution residual");

        prop->add_property_i (20, max_iters_idx, 
                              L"Maximum allowed number of iterations");

        prop->add_property_i (0, iters_idx, 
                              L"Total number of used solver iterations");

        prop->add_property_b (false, success_idx, 
                              L"True if solver successfully convergent");
        
      }
    
    int bicgstab_solver::solve (sp_matrix_t matrix, spv_double sp_rhs, spv_double sp_sol)
    {
      BS_ERROR (matrix, "bicgstab_solve");
      BS_ERROR (sp_rhs->size (), "bicgstab_solve");
      BS_ERROR (sp_sol->size (), "bicgstab_solve");
      BS_ERROR (prop, "bicgstab_solve");

      t_double rho_1, rho_2 = 1, alpha = 1, beta, omega = 1;
      int iter;
      const double epsmac = 1e-24;
      t_double r_norm, b_norm, den_norm, s_norm;
      t_double *rhs = &(*sp_rhs)[0];
      t_double *sol = &(*sp_sol)[0];
      //t_double *x = solution;

      //OMP_TIME_MEASURE_START (bicgstab_solve_timer);

      t_double tol = prop->get_f (tol_idx);
      tol *= tol;
      //resid = prop->get_residuals ();
      //convergence_rate = prop->get_convergence_rate ();

      int max_iter  = prop->get_i (max_iters_idx);
      t_long n         = matrix->get_n_rows () * matrix->get_n_block_size ();

      sp_p->resize (n);
      sp_phat->resize (n);
      sp_s->resize (n);
      sp_shat->resize (n);
      sp_t->resize (n); 
      sp_v->resize (n); 
      sp_r->resize (n); 
      sp_rtilde->resize (n);

      t_double *p              = &(*sp_p)[0];
      t_double *phat           = &(*sp_phat)[0];
      t_double *s              = &(*sp_s)[0];
      t_double *shat           = &(*sp_shat)[0];
      t_double *t              = &(*sp_t)[0];
      t_double *v              = &(*sp_v)[0];
      t_double *r              = &(*sp_r)[0];
      t_double *rtilde         = &(*sp_rtilde)[0];

      prop->set_b (success_idx, false);

      memset (sol, 0, sizeof (t_double) * n);
      memset (r, 0, sizeof (t_double) * n);

      matrix->calc_lin_comb (-1.0, 1.0, sp_sol, sp_rhs, sp_r);
      r_norm = mv_vector_inner_product_n (r, r, n);
      if (r_norm <= tol) // initial guess quite good
        return 0;

      rho_1 = r_norm;
      b_norm = sqrt (mv_vector_inner_product_n (rhs, rhs, n));

      memcpy (p, r, sizeof (t_double) * n);
      memcpy (rtilde, r, sizeof (t_double) * n);
      memset (v, 0, sizeof (t_double) * n);

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
          //printf ("BiCGStab iteration: %d, resid = %le\n", iter, r_norm);
          //fflush (stdout);

          if (iter)
            {
              rho_1 = mv_vector_inner_product_n (r, rtilde, n); //in first iter equals to r_norm
              if (rho_1 == 0) // failure
                {
                  if (den_norm > epsmac)
                    prop->set_f(final_res_idx, r_norm / den_norm);
                  else
                    prop->set_f (final_res_idx, r_norm);
                  //printf ("BiCGStab failure: rho_1 = 0! resid = %le\n", prop->final_resid);
                  bs_throw_exception ("BICGSTAB: Failure - rho_1 == 0");
                }
              beta = (rho_1 / rho_2) * (alpha / omega);
              // p = r + beta * (p - omega * v);
              //AXPY_AYPX (p, -omega, v, beta, r, k, n);
              axpy_aypx_n (p, -omega, v, beta, r, n);
            }

          // phat = M^(-1) * p;
          if (prec)
            {
              if (prec->solve_prec (matrix, sp_p, sp_phat))
                {
                  bs_throw_exception ("BICGSTAB: Preconditioner failed");
                }
            }
          else // no precondition (preconditioner=identity_matrix)
            {
              memcpy (phat, p, sizeof (t_double) * n);
              //phat.assign (p.begin (), p.end ());
            }

          // v = A * phat;
          memset (v, 0, sizeof (t_double) * n);
          //v.assign (n, 0);
          matrix->matrix_vector_product (sp_phat, sp_v);

          alpha = mv_vector_inner_product_n (rtilde, v, n);
          if (alpha > epsmac || alpha < -epsmac)
            alpha = rho_1 / alpha;
          else // failure
            {
              if (den_norm > epsmac)
                prop->set_f(final_res_idx, r_norm / den_norm);
              else
                prop->set_f (final_res_idx, r_norm);

              //printf ("BiCGStab failure: (rtilde, v) = 0! resid = %le\n", prop->final_resid);
              bs_throw_exception ("BICGSTAB: Failure - (rtilde, v) == 0");
            }

          // s = r - alpha * v;

          memcpy (s, r, sizeof (t_double) * n);
          //s.assign (r.begin (), r.end ());
          axpy_n (s, v, -alpha, n);
          //AXPY (s, -alpha, v, k, n);

          //x = x + alpha * phat;
          //AXPY (x, alpha, phat, k, n);
          //axpy (x, phat, alpha);
          axpy_n (sol, phat, alpha, n);

          s_norm = mv_vector_inner_product_n (s, s, n);
          if (s_norm < tol)
            {
              //check convergence
              //matrix->calc_lin_comb (-1.0, 1.0, x, rhs, t);// t is buffer
              matrix->calc_lin_comb (-1.0, 1.0, sp_sol, sp_rhs, sp_t);// t is buffer
              r_norm = mv_vector_inner_product_n (t, t, n);
              if (r_norm <= tol)
                break;
            }

          // shat = M^(-1) * s;
          if (prec)
            {
              if (prec->solve_prec (matrix, sp_s, sp_shat))
                {
                  bs_throw_exception ("BICGSTAB: Preconditioner failed");
                }
            }
          else // no precondition (preconditioner=identity_matrix)
            {
              memcpy (shat, s, sizeof (t_double) * n);
              //shat.assign (s.begin (), s.end ());
            }

          // t = A * shat;
          memset (t, 0, sizeof (t_double) * n);
          //t.assign (n, 0);
          matrix->matrix_vector_product (sp_shat, sp_t);

          // omega = (t,s) / (t,t);
          omega = mv_vector_inner_product_n (t, t, n);
          if (omega > epsmac)
            {
              omega = mv_vector_inner_product_n (t, s, n) / omega;
            }
          else // failure
            {
              if (den_norm > epsmac)
                prop->set_f(final_res_idx, r_norm / den_norm);
              else
                prop->set_f (final_res_idx, r_norm);
              //printf ("BiCGStab failure: (t, t) = 0! resid = %le\n", prop->final_resid);
              bs_throw_exception ("BICGSTAB: Failure - (t, t) == 0");
            }

          if (omega < epsmac) // failure
            {
              if (den_norm > epsmac)
                prop->set_f(final_res_idx, r_norm / den_norm);
              else
                prop->set_f (final_res_idx, r_norm);
              //printf ("BiCGStab failure: omega = 0! resid = %le\n", prop->final_resid);
              bs_throw_exception ("BICGSTAB: Failure - omega == 0");
            }

          //x = x + omega * shat;
          //AXPY (x, omega, shat, k, n);
          //axpy (x, shat, omega);
          axpy_n (sol, shat, omega, n);

          //r = s - omega * t;
          //memcpy (r, s, n * sizeof (t_double));
          //AXPY (r, -omega, t, k, n);
          memcpy (r, s, n * sizeof (t_double));
          //r.assign (s.begin (), s.end ());
          axpy_n (r, t, -omega, n);
          /*
          //additional check convergence
          mv_calc_lin_comb (matrix, -1.0, 1.0, x, rhs, s); // s is buffer
          r_norm = mv_vector_inner_product (s, s, n);
          */
          r_norm = mv_vector_inner_product_n (r, r, n);
          if (r_norm <= tol)
            break;

          rho_2 = rho_1;
        } // end of main loop

      prop->set_i (iters_idx, iter + 1);
      prop->set_b (success_idx, true);

      //printf ("BiCGStab after iteration: %d, resid = %le\n", iter, r_norm);
      /*
      //additional checking convergence
      mv_calc_lin_comb (matrix, -1.0, 1.0, solution, rhs, r);
      r_norm = mv_vector_inner_product (r, r, n);
      */
      if (den_norm > epsmac)
        prop->set_f(final_res_idx, r_norm / den_norm);
      else
        prop->set_f (final_res_idx, r_norm);

      BOSOUT (section::solvers, level::low) << "r_norm = " << r_norm << " r_norm / den_norm = " << r_norm / den_norm << " iter = " << iter << bs_end;

      // printf ("BiCGStab OK! iters = %d, resid = %le\n", prop->iters, prop->final_resid);
      //OMP_TIME_MEASURE_END (bicgstab_solve_timer);

      return 0;
    }

    
    int bicgstab_solver::solve_prec (sp_matrix_t matrix, spv_double rhs, spv_double solution)
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
    bicgstab_solver::setup (sp_matrix_t matrix)
    {
      if (!matrix)
        {
          bs_throw_exception ("BICGStab: Passed matrix is null");
        }

      BS_ASSERT (prop);
      if (prec)
        {
          return prec->setup (matrix);
        }

      return 0;
    }

    //////////////////////////////////////////////////////////////////////////
    BLUE_SKY_TYPE_STD_CREATE (bicgstab_solver);
    BLUE_SKY_TYPE_STD_COPY (bicgstab_solver);

    BLUE_SKY_TYPE_IMPL (bicgstab_solver, lsolver_iface, "bicgstab_solver", "BiCGStab linear solver", "BiCGStab linear solver");
  }
