/**
 * \file gmres_solver.cpp
 * \brief CGS linear solver impl
 * \author Salimgareeva E.M. (algorithm from "A fast Lanczos-type solver for nonsymmetric linear systems", P.Sonneveld)
 * \date 11.02.2009
 * */

#include "bs_kernel.h"
#include "bs_assert.h"

#include "gmres_solver.h"
#include "mv_functions.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "matrix_iface.h"
#include "bos_report.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {

    //////////////////////////////////////////////////////////////////////////
    //  gmres_solver

    //! constructor
    gmres_solver::gmres_solver (bs_type_ctor_param /*param*/)
      : lsolver_iface ()
    {
      prop = BS_KERNEL.create_object ("prop");
      if (!prop)
        {
          bs_throw_exception ("Type (prop) not registered");
        }
      init_prop ();

      sp_w = BS_KERNEL.create_object (v_double::bs_type ());
      sp_r = BS_KERNEL.create_object (v_double::bs_type ());
      sp_s = BS_KERNEL.create_object (v_double::bs_type ());
      sp_c = BS_KERNEL.create_object (v_double::bs_type ());
      sp_rs = BS_KERNEL.create_object (v_double::bs_type ());
      sp_hh = BS_KERNEL.create_object (v_double::bs_type ());
    }

    //! copy constructor
    gmres_solver::gmres_solver(const gmres_solver &solver)
      : bs_refcounter (), lsolver_iface ()
    {
      if (&solver != this)
        *this = solver;
    }

    //! destructor
    gmres_solver::~gmres_solver ()
    {}

    //! set solver's properties
    void gmres_solver::set_prop (sp_prop_t prop_)
    {
      prop = prop_;

      init_prop ();
    }

    //! set properties by default
    void
    gmres_solver::init_prop ()
      {
        prop->add_property_f (1.0e-4, tol_idx,
                              std::string ("Target tolerance for linear solver"));
        prop->add_property_f (1, final_res_idx,
                              std::string ("Solution residual"));
        prop->add_property_i (200, max_iters_idx,
                              std::string ("Maximum allowed number of iterations"));
        prop->add_property_i (0,iters_idx,
                              std::string ("Total number of used solver iterations"));
        prop->add_property_i (30, m_idx,
                              std::string ("Number of vectors used for ortogonalization"));
        prop->add_property_b (false, success_idx,
                              std::string ("True if solver successfully convergent"));
        prop->add_property_i (0, ortonorm_vlen,
                              std::string ("ORTONORM_VLEN (for JACOBIAN)"));
      }

    /**
    * @brief solve for GMRES
    * @param matrix -- input smart pointer to matrix
    * @param sp_rhs -- input smart pointer to right hand side vector
    * @param sp_sol -- input smart pointer to solution vector
    * @return 0 if success
    */
    int gmres_solver::solve (sp_matrix_t matrix, spv_double sp_rhs, spv_double sp_sol)
    {
      BS_ASSERT (matrix);
      BS_ASSERT (sp_rhs->size ());
      BS_ASSERT (sp_sol->size ());
      BS_ASSERT (sp_rhs->size () == sp_sol->size ()) (sp_rhs->size ()) (sp_sol->size ());
      BS_ASSERT (prop);

      //t_double *rs, *hh, *c, *s;
      t_long i, j, k, cnt;
      t_double gamma, t, r_norm, b_norm, den_norm;
      t_double *cur_h;
      const double epsmac = 1.e-16;
      const int m = prop->get_i (m_idx);

      t_double *rhs = &(*sp_rhs)[0];
      t_double *sol = &(*sp_sol)[0];
      // do barrier at mpi or nothing do at seq (base)
      // in future we can hold instance of barrier_t as a data member
      //barrier_t ().barrier (barrier_t::comm_world_v ());

      t_long n = matrix->get_n_rows () * matrix->get_n_block_size ();
      t_double tol = prop->get_f (tol_idx);

      // check workspace
      //wksp.assign (n * (m + 3) + (m + 2) * (m + 1) + 2 * m, 0);

      //vec_p.assign (m + 1, fp_vector_type_t ());

      //vec_p.reserve (m + 1);
      cnt = m + 1;
      int old_n = vec_p.size ();
      vec_p.resize (cnt);
      for (i = old_n; i < cnt; ++i)
        {
          vec_p[i] = BS_KERNEL.create_object (v_double::bs_type ());
        }
      for (i = 0, cnt = m + 1; i < cnt; ++i)
        {
          matrix->init_vector (vec_p[i]);
        }

      matrix->init_vector (sp_w);
      matrix->init_vector (sp_r);

      sp_s->resize (m);
      sp_s->assign (0);
      sp_c->resize (m);
      sp_c->assign (0);
      sp_rs->resize (m + 1);
      sp_rs->assign (0);
      sp_hh->resize ((m + 1) * (m + 1));
      sp_hh->assign (0);

      t_double *vec_s          = &(*sp_s)[0];
      t_double *vec_w          = &(*sp_w)[0];
      t_double *vec_r          = &(*sp_r)[0];
      t_double *vec_c          = &(*sp_c)[0];
      t_double *vec_rs          = &(*sp_rs)[0];
      t_double *vec_hh          = &(*sp_hh)[0];

      // initialize work arrays
      //s = &wksp[0];
      //c = s + m;
      //rs = c + m;
      //hh = rs + m + 1;

      prop->set_b (success_idx, false);

      sp_sol->assign (0);

      //calculate_init_r (n, matrix, solution, rhs, p);
      matrix->calc_lin_comb (-1.0, 1.0, sp_sol, sp_rhs, vec_p[0]);

      b_norm = sqrt (mv_vector_inner_product_n (rhs, rhs, n));
      //BOSOUT (section::solvers, level::low) << "b_norm = " << b_norm << bs_end;


      /* Since it is does not diminish performance, attempt to return an error flag
      and notify users when they supply bad input. */
      t_double *vec_po = &(*vec_p[0])[0];
      r_norm = sqrt (mv_vector_inner_product_n (vec_po, vec_po, n));

      //BOSOUT (section::solvers, level::low) << "r_norm = " << r_norm << bs_end;

      //printf ("Initial_residual %le\n", r_norm / b_norm);

      if (b_norm > epsmac)
        {
          /* convergence criterion |r_i|/|b| <= accuracy if |b| > 0 */
          tol *= b_norm;
          den_norm = b_norm;
        }
      else
        {
          /* convergence criterion |r_i|/|r0| <= accuracy if |b| = 0 */
          tol *= r_norm;
          den_norm = r_norm;
        }
      // set up initial norm and convergense factor
      //lprop->set_relative_factor (den_norm);

      int iter = 0;
      int max_iter = prop->get_i (max_iters_idx);
      for (iter = 0; iter < max_iter;)
        {
          /* initialize first term of hessenberg system */
          vec_rs[0] = r_norm;

          if (r_norm < epsmac || r_norm <= tol)
            break;

          t = 1.0 / r_norm;

          scale_vector_n (vec_po, t, n);

          for (i = 1; i < m && iter < max_iter; ++i, ++iter)
            {
              //barrier_t ().barrier (barrier_t::comm_world_v ());

              //memset (r, 0, sizeof (double) * n);
              if (prec)
                {
                  if (prec->solve_prec (matrix, vec_p[i - 1], sp_r))
                    {
                      bs_throw_exception ("GMRES: Preconditioner failed");
                    }
                }
              else // no precondition (preconditioner=identity_matrix)
                {
                  //vec_r.assign (vec_p[i - 1].begin (), vec_p[i - 1].end ());
                  memcpy (vec_r, &(*vec_p[i - 1])[0], n * sizeof (t_double));
                }

              //barrier_t ().barrier (barrier_t::comm_world_v ());

              vec_p[i]->assign (0);
              matrix->matrix_vector_product (sp_r, vec_p[i]);

              /* modified Gram_Schmidt */
              cur_h = &vec_hh[(i - 1) * (m + 1)];
              t_double *vec_pi = &(*vec_p[i])[0];
              for (j = 0; j < i; ++j)
                {
                  t_double *vec_pj = &(*vec_p[j])[0];
                  cur_h[j] = mv_vector_inner_product_n (vec_pj, vec_pi, n);
                  t = -cur_h[j];

                  axpy_n (vec_pi, vec_pj, t, n);
                }
              t = sqrt (mv_vector_inner_product_n (vec_pi, vec_pi, n));
              cur_h[i] = t;
              if (t > epsmac)
                {
                  t = 1.0 / t;
                  scale_vector_n (vec_pi, t, n);
                }
              /* done with modified Gram_schmidt and Arnoldi step.
              update factorization of hh */
              for (j = 1; j < i; j++)
                {
                  t = cur_h[j - 1];
                  cur_h[j - 1] = vec_c[j - 1] * t + vec_s[j - 1] * cur_h[j];
                  cur_h[j] = -vec_s[j - 1] * t + vec_c[j - 1] * cur_h[j];
                }

              gamma = sqrt (cur_h[i - 1] * cur_h[i - 1] + cur_h[i] * cur_h[i]);

              if (gamma <epsmac)
                gamma = epsmac;

              vec_c[i - 1]  = cur_h[i - 1] / gamma;
              vec_s[i - 1]  = cur_h[i]     / gamma;

              vec_rs[i]     = -vec_s[i - 1] * vec_rs[i - 1];
              vec_rs[i - 1] =  vec_c[i - 1] * vec_rs[i - 1];

              /* determine residual norm */
              cur_h[i - 1] = vec_c[i - 1] * cur_h[i - 1] + vec_s[i - 1] * cur_h[i];

              r_norm = fabs (vec_rs[i]);
              if (r_norm <= tol)
                break;
            }

          if (i == m || iter == max_iter)
            {
              i = i - 1;
            }

          /* now compute solution, first solve upper triangular system */
          vec_rs[i - 1] = vec_rs[i - 1] / vec_hh[i - 1 + (i - 1) * (m + 1)];
          for (k = i - 2; k >= 0; --k)
            {
              t = vec_rs[k];
              for (j = k + 1; j < i; j++)
                {
                  t -= vec_hh[k + j * (m + 1)] * vec_rs[j];
                }
              vec_rs[k] = t / vec_hh[k + k * (m + 1)];
            }

          /* form linear combination of p's to get solution */
          //vec_w.assign (vec_p[0].begin (), vec_p[0].end ());
          memcpy (&(vec_w[0]), vec_po, sizeof (t_double) * n);
          t = vec_rs[0];

          scale_vector_n (vec_w, t, n);

          for (j = 1; j < i; j++)
            {
              t = vec_rs[j];
              t_double *vec_pj = &(*vec_p[j])[0];
              axpy_n (vec_w, vec_pj, t, n);
            }

          //memset (r, 0, sizeof (double) * n);

          //barrier_t ().barrier (barrier_t::comm_world_v ());

          if (prec)
            {
              if (prec->solve_prec (matrix, sp_w, sp_r))
                {
                  bs_throw_exception ("GMRES: Preconditioner failed");
                }
            }
          else // no precondition (preconditioner=identity_matrix)
            {
              memcpy (vec_r, vec_w, n * sizeof (t_double));
              //vec_r.assign (vec_w.begin (), vec_w.end ());
            }

          //barrier_t ().barrier (barrier_t::comm_world_v ());

          axpy_n (sol, vec_r, (t_double)1.0, n);

          /* check for convergence, evaluate actual residual */
          if (r_norm <= tol)
            {
              //break;
#if 1
              //calculate_init_r (n, matrix, solution, rhs, r);
              matrix->calc_lin_comb (-1.0, 1.0, sp_sol, sp_rhs, sp_r);

              r_norm = sqrt (mv_vector_inner_product_n (vec_r, vec_r, n));

              if (r_norm <= tol)
                break;
              else
                {
                  //vec_p[0].assign (vec_r.begin (), vec_r.end ());
                  memcpy (vec_po, vec_r, sizeof (t_double) * n);
                  i = 0;
                  ++iter;
                }
#endif //0
            }

          /* compute residual vector and continue loop */
          for (j = i; j > 0; --j)
            {
              vec_rs[j - 1] = -vec_s[j - 1] * vec_rs[j];
              vec_rs[j] = vec_c[j - 1] * vec_rs[j];
            }

          if (i)
            {
              t = vec_rs[0];
              scale_vector_n (vec_po, t, n);

            }
          for (j = 1; j < i + 1; ++j)
            {
              t = vec_rs[j];
              t_double *vec_pj = &(*vec_p[j])[0];
              axpy_n (vec_po, vec_pj, t, n);
            }
        }

      prop->set_i (iters_idx, iter + 1);
      if (iter < max_iter)
        prop->set_b (success_idx, true);

      if (den_norm > 1.0e-12)
        prop->set_f (final_res_idx, r_norm / den_norm);
      else
        prop->set_f (final_res_idx, r_norm);

#ifdef _DEBUG
      BOSOUT (section::solvers, level::low) << "r_norm = " << r_norm << " r_norm / den_norm = " << r_norm / den_norm << " iter = " << (iter + 1) << bs_end;
#endif

      //barrier_t ().barrier (barrier_t::comm_world_v ());
      return 0;
    }


    int gmres_solver::solve_prec (sp_matrix_t matrix, spv_double rhs, spv_double sol)
    {
      return solve (matrix, rhs, sol);
    }

    /**
    * @brief setup for GMRES
    *
    * @param matrix -- input matrix
    *
    * @return 0 if success
    */
    int
    gmres_solver::setup (sp_matrix_t matrix)
    {
      if (!matrix)
        {
          bs_throw_exception ("GMRES: Passed matrix is null");
        }

      BS_ASSERT (prop);
      if (prec)
        {
          return prec->setup (matrix);
        }

      return 0;
    }

    //////////////////////////////////////////////////////////////////////////
    BLUE_SKY_TYPE_STD_CREATE (gmres_solver);
    BLUE_SKY_TYPE_STD_COPY (gmres_solver);

    BLUE_SKY_TYPE_IMPL (gmres_solver, lsolver_iface, "gmres_solver", "GMRES linear solver", "GMRES linear solver");
  }
