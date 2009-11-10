/**
 * \file cgs.cpp
 * \brief CGS linear solver impl
 * \author Salimgareeva E.M. (algorithm from "A fast Lanczos-type solver for nonsymmetric linear systems", P.Sonneveld)
 * \date 11.02.2009
 * */
#include "bs_base_linear_solvers_stdafx.h"
#include "linear_solvers.h"
#include "cgs.h"
#include "lin_solv_macro.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "save_seq_vector.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {

    //////////////////////////////////////////////////////////////////////////
    //  cgs_solver

    //! constructor
    template <class strat_t>
    cgs_solver<strat_t>::cgs_solver (bs_type_ctor_param param)
      : linear_solver_base<strat_t> (param)
    {}

    //! copy constructor
    template <class strat_t>
    cgs_solver<strat_t>::cgs_solver(const cgs_solver &solver)
      : bs_refcounter (), linear_solver_base<strat_t> (solver)
    {
      if (&solver != this)
        *this = solver;
    }

    //! destructor
    template <class strat_t>
    cgs_solver<strat_t>::~cgs_solver ()
    {}


    template <class strat_t>
    int cgs_solver<strat_t>::solve(matrix_t *matrix, rhs_item_array_t &rhs, item_array_t &solution)
    {
      return templ_solve (matrix, rhs, solution);
    }

    template <class strat_t>
    int cgs_solver<strat_t>::solve_prec(matrix_t *matrix, item_array_t &rhs, item_array_t &solution)
    {
      return templ_solve (matrix, rhs, solution);
    }


    /*!
    * \brief CGS linear solver
    *
    * \param[in] matrix -- pointer to the matrix
    * \param[in] rhs -- right hand side
    * \param[out] solution -- solution
    *
    * \return 0 if success
    */

    //template <class strat_t>
    //int cgs_solver<strat_t>::solve(matrix_t *matrix, rhs_item_array_t &rhs, item_array_t &solution)
    template <class strat_t> template <class rhs_t> int
    cgs_solver<strat_t>::templ_solve (matrix_t *matrix, rhs_t &rhs, item_array_t &solution)
    {
#ifdef _DEBUG
      BOSOUT (section::solvers, level::debug) << "CGS\n" << bs_end;
#endif
      typedef item_t fp_type;

      BS_ERROR (rhs.size (), "cgs_solve");
      BS_ERROR (solution.size (), "cgs_solve");
      BS_ERROR (base_t::prop, "cgs_solve");

      const smart_ptr<linear_solver_prop> &lprop(this->prop);

      fp_type rho_1, rho_2 = 1, alpha = 1, beta, sigma;
      int iter;
      const double epsmac = 1e-24;
      fp_type r_norm, b_norm, den_norm;
      //fp_type *x = solution;

      //OMP_TIME_MEASURE_START (cgs_solve_timer);

      item_t tol = this->prop->get_tolerance ();
      tol *= tol;
      //resid = prop->get_residuals ();
      //convergence_rate = prop->get_convergence_rate ();

      index_t max_iter  = this->prop->get_max_iters ();
      index_t n         = matrix->n_rows * matrix->n_block_size;

      item_array_t p (n);
      item_array_t phat (n);
      item_array_t v (n);
      item_array_t tmp (n);
      item_array_t q (n);
      item_array_t u (n);
      item_array_t d (n);
      item_array_t dhat (n);
      item_array_t r (n);
      item_array_t rtilde (n);
      item_array_t r_old (n);

      lprop->set_success (0);

      // solution = {0}
      assign (solution, n, 0);

      // r = {0}
      assign (r, n, 0);
       // TODO:paste
      assign (tmp, n, 0);
      assign (p, n, 0);
      assign (v, n, 0);
      assign (q, n, 0);

      // p0 = u0 = r0;
      //memcpy (p, r, n * sizeof (double));
      u.assign (r.begin (), r.end ());
      // TODO:end

      // r = Ax0 - b
      matrix->calc_lin_comb (-1.0, 1.0, solution, rhs, r);
      rtilde.assign (r.begin (), r.end ());

      //tools::save_seq_vector (tools::string_formater ("1_well_bhp.%s.txt", it->first).str).save (it->second);

      r_norm = bos_helper::mv_vector_inner_product (r, r);

      if (r_norm <= tol) // initial guess quite good
        return 0;

      rho_1 = r_norm;
      b_norm = sqrt (bos_helper::mv_vector_inner_product (rhs, rhs));

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
      lprop->set_relative_factor (den_norm);

      // main loop
      for (iter = 0; iter < max_iter; ++iter)
        {
          //printf ("CGS iteration: %d, resid = %le\n", iter, r_norm);
          //fflush (stdout);
          // TODO: paste
          if (iter)
            {
              //rho_1 = (r,rtilde)
              rho_1 = bos_helper::mv_vector_inner_product (r, rtilde); //in first iter equals to r_norm
              if (rho_1 == 0) // failure
                {
                  if (den_norm > epsmac)
                    lprop->set_final_resid (r_norm / den_norm);
                  else
                    lprop->set_final_resid (r_norm);
                  bs_throw_exception ("CGS: Failure - rho_1 == 0");
                }
            }

           beta = rho_1/rho_2; // beta = rho_n/rho_n-1
           rho_2 = rho_1;

           // u = r + beta*q
           sum_vector(r, 1., q, beta, u);
           // tmp = q+beta*p_old
           sum_vector(q, 1., p, beta, tmp);
           // p_new = u + beta*tmp
           sum_vector(u, 1., tmp, beta, p);

           //temp_p.assign (p.begin (), p.end ());
           if (this->prec)
             {
               if (base_t::prec->solve_prec (matrix, p, phat))
                 {
                   bs_throw_exception ("CGS: Preconditioner failed");
                 }
             }
           else // no precondition (preconditioner=identity_matrix)
             {
               phat.assign (p.begin (), p.end ());
             }

          // v = A * phat = A * p, if no precondition;
          assign (v, n, 0);

          matrix->matrix_vector_product (phat, v);
          // sigma = (v,rtilde)
          sigma = bos_helper::mv_vector_inner_product (rtilde, v);

          if (sigma > epsmac || sigma < -epsmac)
          // alpha = rho_1/sigma
            alpha = rho_1/sigma;
          else // failure
            {
              if (den_norm > epsmac)
                lprop->set_final_resid (r_norm / den_norm);
              else
                lprop->set_final_resid (r_norm);
              
              bs_throw_exception ("CGS: Failure - sigma == 0");
            }

          // q = u - alpha*v
          sum_vector(u, 1., v, -alpha, q);
          // d = u + q
          sum_vector(u, 1., q, 1., d);

          // dhat = M^(-1) * d;
          //temp_d.assign (d.begin (), d.end ());
          if (this->prec)
            {
              if(base_t::prec->solve_prec (matrix, d, dhat))
	              {
	                bs_throw_exception ("CGS: Preconditioner failed");
	              }
            }
          else // no precondition (preconditioner=identity_matrix)
            {
              dhat.assign (d.begin (), d.end ());
            }

         assign (tmp, n, 0);
         // tmp = A*d
         matrix->matrix_vector_product (dhat, tmp);

         // r = r - alpha*tmp
         sum_vector(r, 1., tmp, -alpha, r);
         // x = x + alpha*dhat
         sum_vector(solution, 1., dhat, alpha, solution);

         r_norm = bos_helper::mv_vector_inner_product (r, r, n);


         if (r_norm <= tol) // && check_resid_for_matbalance (n_rows, nb, r, matb_tol))
          break;
     }

     //tools::save_seq_vector ("solution.txt").save(solution);

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

      //printf ("CGS OK! iters = %d, resid = %le\n", lprop->iters, lprop->final_resid);
      //OMP_TIME_MEASURE_END (bicgstab_solve_timer);

      return 0;
    }

    /**
    * @brief setup for CGS
    *
    * @param matrix -- input matrix
    *
    * @return 0 if success
    */
    template <class strat_t> int
    cgs_solver<strat_t>::setup (matrix_t *matrix)
    {
      if (!matrix)
        {
          bs_throw_exception ("CGS: Passed matrix is null");
        }

      BS_ASSERT (base_t::prop);
      if (base_t::prec)
        {
          return base_t::prec->setup (matrix);
        }

      return 0;
    }

    //////////////////////////////////////////////////////////////////////////
    BLUE_SKY_TYPE_STD_CREATE_T_DEF(cgs_solver, (class));
    BLUE_SKY_TYPE_STD_COPY_T_DEF(cgs_solver, (class));

    BLUE_SKY_TYPE_IMPL_T_EXT(1, (cgs_solver<base_strategy_fi>) , 1, (linear_solver_base<base_strategy_fi>), "cgs_solver_base_fi", "CGS linear solver", "CGS linear solver", false);
    BLUE_SKY_TYPE_IMPL_T_EXT(1, (cgs_solver<base_strategy_di>) , 1, (linear_solver_base<base_strategy_di>), "cgs_solver_base_di", "CGS linear solver", "CGS linear solver", false);
    BLUE_SKY_TYPE_IMPL_T_EXT(1, (cgs_solver<base_strategy_mixi>) , 1, (linear_solver_base<base_strategy_mixi>), "cgs_solver_base_mixi", "CGS linear solver", "CGS linear solver", false);
  }
