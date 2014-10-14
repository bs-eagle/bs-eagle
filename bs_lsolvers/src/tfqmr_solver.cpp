/**
 * \file tfqmr_solver.cpp
 * \brief TFQMR linear solver impl
 * \author Salimgareeva E.M. (algorithm from "A fast Lanczos-type solver for nonsymmetric linear systems", P.Sonneveld)
 * \date 11.02.2009
 * */

#include "tfqmr_solver.h"
#include "mv_functions.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "bos_report.h"
#include "matrix_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {

    //////////////////////////////////////////////////////////////////////////
    //  tfqmr_solver

    //! constructor
    
    tfqmr_solver::tfqmr_solver (bs_type_ctor_param /*param*/)
      : lsolver_iface ()
    {
      prop = BS_KERNEL.create_object ("prop");
      if (!prop)
        {
          bs_throw_exception ("Type (prop) not registered");
        }
      init_prop ();
      sp_p      = BS_KERNEL.create_object (v_double::bs_type ());
      sp_v      = BS_KERNEL.create_object (v_double::bs_type ());
      sp_w      = BS_KERNEL.create_object (v_double::bs_type ());
      sp_u      = BS_KERNEL.create_object (v_double::bs_type ());
      sp_q      = BS_KERNEL.create_object (v_double::bs_type ());
      sp_d      = BS_KERNEL.create_object (v_double::bs_type ());
      sp_res    = BS_KERNEL.create_object (v_double::bs_type ());
      sp_r      = BS_KERNEL.create_object (v_double::bs_type ());
      sp_rtilde = BS_KERNEL.create_object (v_double::bs_type ());
      sp_tmp    = BS_KERNEL.create_object (v_double::bs_type ());
      sp_rhat   = BS_KERNEL.create_object (v_double::bs_type ());
      sp_y      = BS_KERNEL.create_object (v_double::bs_type ());
    }

     void
    tfqmr_solver::init_prop ()
      {
        prop->add_property_f (1.0e-4, tol_idx, L"Target tolerance for linear solver");
        prop->add_property_f (1, final_res_idx, L"Solution residual");
        prop->add_property_i (20, max_iters_idx, L"Maximum allowed number of iterations");
        prop->add_property_i (0, iters_idx, L"Total number of used solver iterations");
        prop->add_property_b (false, success_idx, L"True if solver successfully convergent");
      }
    //! copy constructor
    
    tfqmr_solver::tfqmr_solver(const tfqmr_solver &solver)
      : bs_refcounter (), lsolver_iface ()
    {
      if (&solver != this)
        *this = solver;
    }

    //! destructor
    
    tfqmr_solver::~tfqmr_solver ()
    {}

    //! set solver's properties
    
    void tfqmr_solver::set_prop (sp_prop_t prop_)
    {
      prop = prop_;

      init_prop ();
    }

    
    int tfqmr_solver::solve (sp_matrix_t matrix, spv_double sp_rhs, spv_double sp_sol)
    {
      BOSOUT (section::solvers, level::debug) << "TFQMR\n" << bs_end;

      BS_ERROR (matrix, "tfqmr_solve");
      BS_ERROR (sp_rhs->size (), "tfqmr_solve");
      BS_ERROR (sp_sol->size (), "tfqmr_solve");
      BS_ERROR (prop, "tfqmr_solve");


      t_double rho_1, rho_2 = 1, alpha = 1, beta, sigma;
      int iter;
      const double epsmac = 1e-24;
      t_double r_norm, b_norm, den_norm, w_norm, eta, nu, tau, c;
      //fp_type *x = solution;

      //OMP_TIME_MEASURE_START (tfqmr_solve_timer);

      t_double *rhs = &(*sp_rhs)[0];
      t_double *sol = &(*sp_sol)[0];
      t_double tol = prop->get_f (tol_idx);
      tol *= tol;
      //resid = prop->get_residuals ();
      //convergence_rate = prop->get_convergence_rate ();

      int max_iter  = prop->get_i (max_iters_idx);
      t_long n         = matrix->get_n_rows () * matrix->get_n_block_size ();

      sp_p->resize (n);
      sp_v->resize (n);
      sp_w->resize (n);
      sp_u->resize (n);
      sp_q->resize (n);
      sp_d->resize (n);
      sp_res->resize (n);
      sp_r->resize (n);
      sp_rtilde->resize (n);
      sp_tmp->resize (n);
      sp_rhat->resize (n);
      sp_y->resize (n);
      //x_cgs = y + n;

      prop->set_b (success_idx, false);

      // solution = {0}
      sp_sol->assign (0);
      // r = {0}
      sp_r->assign (0);
       // TODO:paste
      sp_tmp->assign (0);
      sp_p->assign (0);
      sp_v->assign (0);
      sp_q->assign (0);
      sp_d->assign (0);

      t_double *r              = &(*sp_r)[0];
      t_double *p              = &(*sp_p)[0];
      t_double *v              = &(*sp_v)[0];
      t_double *w              = &(*sp_w)[0];
      t_double *u              = &(*sp_u)[0];
      t_double *q              = &(*sp_q)[0];
      t_double *d              = &(*sp_d)[0];
      t_double *res            = &(*sp_res)[0];
      t_double *rtilde         = &(*sp_rtilde)[0];
      t_double *tmp            = &(*sp_tmp)[0];
      t_double *rhat           = &(*sp_rhat)[0];
      t_double *y              = &(*sp_y)[0];

      // r = Ax0 - b
      matrix->calc_lin_comb (-1.0, 1.0, sp_sol, sp_rhs, sp_r);
      memcpy (rtilde, r, sizeof (t_double) * n);
      //rtilde.assign (r.begin (), r.end ());

      // p0 = u0 = r0;
      //memcpy (p, r, n * sizeof (double));
      memcpy (u, r, sizeof (t_double) * n);
      memcpy (p, r, sizeof (t_double) * n);
      memcpy (w, r, sizeof (t_double) * n);
      //u.assign (r.begin (), r.end ());
      //p.assign (r.begin (), r.end ());
      //w.assign (r.begin (), r.end ());

      // tmp = M^(-1) * u;
      if (prec)
        {
          if (prec->solve_prec (matrix, sp_u, sp_tmp))
            {
              bs_throw_exception ("TFQMR: Preconditioner failed");
            }
              memcpy (u, tmp, sizeof (t_double) * n);
              memcpy (p, tmp, sizeof (t_double) * n);
	      //u.assign (tmp.begin (), tmp.end ());
	      //p.assign (u.begin (), u.end ());
        }

      matrix->matrix_vector_product (sp_p, sp_v);

      //tools::save_seq_vector (tools::string_formater ("1_well_bhp.%s.txt", it->first).str).save (it->second);

      r_norm = mv_vector_inner_product_n (r, r, n);


      if (r_norm <= tol) // initial guess quite good
        return 0;

      tau = sqrt (r_norm);
      rho_1 = r_norm;
      rho_2 = r_norm;
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

      // set up initial norm and convergense factor
      //prop->set_relative_factor (den_norm);

      int m, count;
      // main loop
      for (iter = 0; iter < max_iter; ++iter)
        {
          //printf ("TFQMR iteration: %d, resid = %le\n", iter, r_norm);
          //fflush (stdout);
          // TODO: paste
          if (iter)
            {
             //rho_1 = mv_vector_inner_product (r, rtilde, n);//in first iter equals to r_norm
              if (rho_1 == 0) // failure
                {
                  if (den_norm > epsmac)
                    prop->set_f (final_res_idx, r_norm / den_norm);
                  else
                    prop->set_f (final_res_idx, r_norm);
                  
                  bs_throw_exception ("TFQMR: Failure - rho_1 == 0");
                }
               sum_vector_n (u, (t_double)1., res, beta, p, n); //p[n] = u[n]+beta*res

               //v.assign (n, 0);
               memset (v, 0, sizeof (t_double) * n);
               matrix->matrix_vector_product (sp_p, sp_v); //v[n]=Ap[n]
             }

           sigma = mv_vector_inner_product_n (rtilde, v, n); //sigma=(rtilde,v[n-1])

           alpha = rho_1/sigma;

           // tmp = M^(-1)*v
           if (prec)
             {
               if (prec->solve_prec (matrix, sp_v, sp_tmp))
                {
                  bs_throw_exception ("TFQMR: Preconditioner failed");
                }
                   memcpy (v, tmp, sizeof (t_double) * n);
	           //v.assign (tmp.begin (), tmp.end ());
             }

           sum_vector_n (u, (t_double)1., v, -alpha, q, n); //q[n] = u[n-1]-alpha*v[n-1]
           sum_vector_n (u, (t_double)1., q, (t_double)1., res, n); //res = u[n-1]+q[n]

           //tmp.assign (n, 0);
           memset (tmp, 0, sizeof (t_double) * n);
           matrix->matrix_vector_product (sp_res, sp_tmp);// tmp=A*res
           sum_vector_n (r, (t_double)1., tmp, -alpha, r, n);// r=r-alpha*res

           //r_norm_old = r_norm;
           r_norm = mv_vector_inner_product_n (r, r, n);

           for (m = 1; m <= 2 ; m++)
             {
               if (m == 1) // m is odd
                 {
                   memcpy (y, u, sizeof (t_double) * n);
                   //y.assign (u.begin (), u.end ());
                   w_norm = sqrt(r_norm * r_norm);
                 }
               else // m is even
                 {
                   memcpy (y, q, sizeof (t_double) * n);
                   //y.assign (q.begin (), q.end ());
                   w_norm = sqrt(r_norm);
                 }

               sum_vector_n (y, (t_double)1., d, eta*nu*nu/alpha, d, n); //d[m] = y[m] + (eta[m-1]*nu[m-1]^2/alpha[n-1])*d[m-1]
               nu = w_norm/tau; //nu[m]=||w[m+1]||/tau[m-1]
               c = 1./sqrt (1. + nu*nu);
               tau = tau*c*nu; //tau[m]=tau[m-1]nu[m]c[m]
               eta = c*c*alpha; //eta[m]=c[m]^2*alpha[n-1]
               //SUM_VECTOR(x,d,1,alpha,x_cgs,k,n); //x_cgs[n] = x[2n-1]+alpha[n-1]*d[2n]
               sum_vector_n (sol, (t_double)1., d, eta, sol, n); //x[m] = x[m-1]+eta[m]*d[m]
               if (r_norm <= tol)
                 {
                   count = 1;
                   break;
                 }
             }

           if (r_norm <= tol)
             break;

           rho_1 = mv_vector_inner_product_n (r, rtilde, n);//in first iter equals to r_norm
           beta = rho_1 / rho_2;

           // rhat = M^(-1) * r;
           if (prec)
             {
               if (prec->solve_prec (matrix, sp_r, sp_rhat))
                {
                  bs_throw_exception ("TFQMR: Preconditioner failed");
                }
             }
           else // no precondition (preconditioner=identity_matrix)
             {
               memcpy (rhat, r, sizeof (t_double) * n);
               //rhat.assign (r.begin (), r.end ());
             }

          sum_vector_n (rhat, (t_double)1., q, beta, u, n); //u[n] = r[n]+beta*q[n]
          sum_vector_n (q, (t_double)1., p, beta, res, n); //res = q[n]+beta*p[n-1]

          rho_2 = rho_1;
     }

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

      //printf ("TFQMR OK! iters = %d, resid = %le\n", lprop->iters, lprop->final_resid);
      //OMP_TIME_MEASURE_END (tfqmr_solve_timer);

      return 0;
    }

    
    int tfqmr_solver::solve_prec (sp_matrix_t matrix, spv_double rhs, spv_double solution)
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
    tfqmr_solver::setup (sp_matrix_t matrix)
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
    BLUE_SKY_TYPE_STD_CREATE (tfqmr_solver);
    BLUE_SKY_TYPE_STD_COPY (tfqmr_solver);

    BLUE_SKY_TYPE_IMPL (tfqmr_solver, lsolver_iface, "tfqmr_solver", "TFQMR linear solver", "TFQMR linear solver");
  }
