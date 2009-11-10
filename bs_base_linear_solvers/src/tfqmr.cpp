/**
 * \file tfqmr.cpp
 * \brief TFQMR linear solver impl
 * \author Salimgareeva E.M. (algorithm from "A fast Lanczos-type solver for nonsymmetric linear systems", P.Sonneveld)
 * \date 11.02.2009
 * */
#include "bs_base_linear_solvers_stdafx.h"
#include "linear_solvers.h"
#include "tfqmr.h"
#include "lin_solv_macro.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "save_seq_vector.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {

    //////////////////////////////////////////////////////////////////////////
    //  tfqmr_solver

    //! constructor
    template <class strat_t>
    tfqmr_solver<strat_t>::tfqmr_solver (bs_type_ctor_param param)
      : linear_solver_base<strat_t> (param)
    {}

    //! copy constructor
    template <class strat_t>
    tfqmr_solver<strat_t>::tfqmr_solver(const tfqmr_solver &solver)
      : bs_refcounter (), linear_solver_base<strat_t> (solver)
    {
      if (&solver != this)
        *this = solver;
    }

    //! destructor
    template <class strat_t>
    tfqmr_solver<strat_t>::~tfqmr_solver ()
    {}


    template <class strat_t>
    int tfqmr_solver<strat_t>::solve(matrix_t *matrix, rhs_item_array_t &rhs, item_array_t &solution)
    {
      return templ_solve (matrix, rhs, solution);
    }

    template <class strat_t>
    int tfqmr_solver<strat_t>::solve_prec(matrix_t *matrix, item_array_t &rhs, item_array_t &solution)
    {
      return templ_solve (matrix, rhs, solution);
    }



    /*!
    * \brief TFQMR linear solver
    *
    * \param[in] matrix -- pointer to the matrix
    * \param[in] rhs -- right hand side
    * \param[out] solution -- solution
    *
    * \return 0 if success
    */
    template <class strat_t> template <class rhs_t>
    int tfqmr_solver<strat_t>::templ_solve(matrix_t *matrix1, rhs_t &rhs, item_array_t &solution)
    {
      BOSOUT (section::solvers, level::debug) << "TFQMR\n" << bs_end;
      typedef item_t fp_type;

      BS_ERROR (matrix1, "tfqmr_solve");
      setup_preconditioner<matrix_t> *setup_prec = dynamic_cast<setup_preconditioner<matrix_t> *> (static_cast <matrix_t *> (matrix1));

      sp_bcsr_matrix_t matrix;
      if (setup_prec)
      {
        setup_prec->prepare_matrix ();
        matrix = setup_prec->get_prepared_matrix ();
        BS_ASSERT (matrix);
      }
      else
      {
        matrix = sp_bcsr_matrix_t (matrix1, bs_dynamic_cast ());
      }

      BS_ERROR (rhs.size (), "tfqmr_solve");
      BS_ERROR (solution.size (), "tfqmr_solve");
      BS_ERROR (base_t::prop, "tfqmr_solve");

      const smart_ptr<linear_solver_prop> &lprop(this->prop);

      fp_type rho_1, rho_2 = 1, alpha = 1, beta, sigma;
      int iter;
      const double epsmac = 1e-24;
      fp_type r_norm = 0;
      fp_type b_norm = 0;
      fp_type den_norm = 0;
      fp_type r_norm_old = 0;
      fp_type w_norm = 0;
      fp_type eta = 0;
      fp_type nu = 0;
      fp_type tau = 0;
      fp_type c = 0;
      //fp_type *x = solution;

      //OMP_TIME_MEASURE_START (tfqmr_solve_timer);

      item_t tol = this->prop->get_tolerance ();
      tol *= tol;
      //resid = prop->get_residuals ();
      //convergence_rate = prop->get_convergence_rate ();

      index_t max_iter  = this->prop->get_max_iters ();
      index_t n         = matrix->n_rows * matrix->n_block_size;

      item_array_t p (n);
      item_array_t v (n);
      item_array_t w (n);
      item_array_t u (n);
      item_array_t q (n);
      item_array_t d (n);
      item_array_t res (n);
      item_array_t r (n);
      item_array_t rtilde (n);
      item_array_t tmp (n);
      item_array_t rhat (n);
      item_array_t y (n);
      //x_cgs = y + n;

      lprop->set_success (0);

      // solution = {0}
      solution.assign (n, 0);
      // r = {0}
      r.assign (n, 0);
       // TODO:paste
      tmp.assign (n, 0);
      p.assign (n, 0);
      v.assign (n, 0);
      q.assign (n, 0);
      d.assign (n, 0);

      // r = Ax0 - b
      matrix->calc_lin_comb (-1.0, 1.0, solution, rhs, r);
      rtilde.assign (r.begin (), r.end ());

      // p0 = u0 = r0;
      //memcpy (p, r, n * sizeof (double));
      u.assign (r.begin (), r.end ());
      p.assign (r.begin (), r.end ());
      w.assign (r.begin (), r.end ());

      // tmp = M^(-1) * u;
      if (this->prec)
        {
          if (base_t::prec->solve_prec (matrix, u, tmp))
            {
              bs_throw_exception ("TFQMR: Preconditioner failed");
            }
	      u.assign (tmp.begin (), tmp.end ());
	      p.assign (u.begin (), u.end ());
        }

      matrix->matrix_vector_product (p, v);

      //tools::save_seq_vector (tools::string_formater ("1_well_bhp.%s.txt", it->first).str).save (it->second);

      r_norm = bos_helper::mv_vector_inner_product (r, r);


      if (r_norm <= tol) // initial guess quite good
        return 0;

      tau = sqrt (r_norm);
      rho_1 = r_norm;
      rho_2 = r_norm;
      b_norm = sqrt (bos_helper::mv_vector_inner_product (rhs, rhs));


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
      lprop->set_relative_factor (den_norm);

      index_t m, count;
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
                    lprop->set_final_resid (r_norm / den_norm);
                  else
                    lprop->set_final_resid (r_norm);
                  
                  bs_throw_exception ("TFQMR: Failure - rho_1 == 0");
                }
               sum_vector (u, 1., res, beta, p); //p[n] = u[n]+beta*res

               v.assign (n, 0);
               matrix->matrix_vector_product (p, v); //v[n]=Ap[n]
             }

           sigma = bos_helper::mv_vector_inner_product (rtilde, v); //sigma=(rtilde,v[n-1])

           alpha = rho_1/sigma;

           // tmp = M^(-1)*v
           if (this->prec)
             {
               if (base_t::prec->solve_prec (matrix, v, tmp))
	              {
                  bs_throw_exception ("TFQMR: Preconditioner failed");
                }
	           v.assign (tmp.begin (), tmp.end ());
             }

           sum_vector (u, 1., v, -alpha, q); //q[n] = u[n-1]-alpha*v[n-1]
           sum_vector (u, 1., q, 1., res); //res = u[n-1]+q[n]

           tmp.assign (n, 0);
           matrix->matrix_vector_product (res, tmp);// tmp=A*res
           sum_vector (r, 1., tmp, -alpha, r);// r=r-alpha*res

           //r_norm_old = r_norm;
           r_norm = bos_helper::mv_vector_inner_product (r, r);

           for (m = 1; m <= 2 ; m++)
             {
               if (m == 1) // m is odd
                 {
                   y.assign (u.begin (), u.end ());
                   w_norm = sqrt(r_norm * r_norm_old);
                 }
               else // m is even
                 {
                   y.assign (q.begin (), q.end ());
                   w_norm = sqrt(r_norm);
                 }

               sum_vector (y, 1., d, eta*nu*nu/alpha, d); //d[m] = y[m] + (eta[m-1]*nu[m-1]^2/alpha[n-1])*d[m-1]
               nu = w_norm/tau; //nu[m]=||w[m+1]||/tau[m-1]
               c = 1./sqrt (1. + nu*nu);
               tau = tau*c*nu; //tau[m]=tau[m-1]nu[m]c[m]
               eta = c*c*alpha; //eta[m]=c[m]^2*alpha[n-1]
               //SUM_VECTOR(x,d,1,alpha,x_cgs,k,n); //x_cgs[n] = x[2n-1]+alpha[n-1]*d[2n]
               sum_vector (solution, 1., d, eta, solution); //x[m] = x[m-1]+eta[m]*d[m]
               if (r_norm <= tol)
                 {
                   count = 1;
                   break;
                 }
             }

           if (r_norm <= tol)
             break;

           rho_1 = bos_helper::mv_vector_inner_product (r, rtilde, n);//in first iter equals to r_norm
           beta = rho_1 / rho_2;

           // rhat = M^(-1) * r;
           if (this->prec)
             {
               if (base_t::prec->solve_prec (matrix, r, rhat))
	              {
                  bs_throw_exception ("TFQMR: Preconditioner failed");
                }
             }
           else // no precondition (preconditioner=identity_matrix)
             {
               rhat.assign (r.begin (), r.end ());
             }

          sum_vector (rhat, 1., q, beta, u); //u[n] = r[n]+beta*q[n]
          sum_vector (q, 1., p, beta, res); //res = q[n]+beta*p[n-1]

          rho_2 = rho_1;
     }

     //TODO: end
     lprop->set_iters (iter + 1);
     lprop->set_success (1);

      /*
      //additional checking convergence
      mv_calc_lin_comb (matrix, -1.0, 1.0, solution, rhs, r);
      r_norm = mv_vector_inner_product (r, r, n);
      */
      if (den_norm > epsmac)
        lprop->set_final_resid (r_norm / den_norm);
      else
        lprop->set_final_resid (r_norm);

      //printf ("TFQMR OK! iters = %d, resid = %le\n", lprop->iters, lprop->final_resid);
      //OMP_TIME_MEASURE_END (tfqmr_solve_timer);

      return 0;
    }

    /**
    * @brief setup for TFQMR
    *
    * @param matrix -- input matrix
    *
    * @return 0 if success
    */
    template <class strat_t> int
    tfqmr_solver<strat_t>::setup (matrix_t *matrix)
    {
      if (!matrix)
        {
          bs_throw_exception ("TFQMR: Passed matrix is null");
        }

      BS_ASSERT (base_t::prop);
      if (base_t::prec)
        {
          return base_t::prec->setup (matrix);
        }

      return 0;
    }

    //////////////////////////////////////////////////////////////////////////
    BLUE_SKY_TYPE_STD_CREATE_T_DEF(tfqmr_solver, (class));
    BLUE_SKY_TYPE_STD_COPY_T_DEF(tfqmr_solver, (class));

    BLUE_SKY_TYPE_IMPL_T_EXT(1, (tfqmr_solver<base_strategy_fi>) , 1, (linear_solver_base<base_strategy_fi>), "tfqmr_solver_base_fi", "tfqmr linear solver", "TFQMR linear solver", false);
    BLUE_SKY_TYPE_IMPL_T_EXT(1, (tfqmr_solver<base_strategy_di>) , 1, (linear_solver_base<base_strategy_di>), "tfqmr_solver_base_di", "tfqmr linear solver", "TFQMR linear solver", false);
    BLUE_SKY_TYPE_IMPL_T_EXT(1, (tfqmr_solver<base_strategy_mixi>) , 1, (linear_solver_base<base_strategy_mixi>), "tfqmr_solver_base_mixi", "TFQMR linear solver", "TFQMR linear solver", false);


  }
