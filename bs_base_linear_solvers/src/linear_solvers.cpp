/*!
* \file linear_solvers.cpp
* \brief implementation of linear solvers
* \author Borschuk Oleg
* \date 2006-07-26
*/
#include "bs_base_linear_solvers_stdafx.h"

#include "linear_solvers.h"
#include "lin_solv_macro.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "dummy_base.h"
#include "bos_report.h"
#include "vector_assign.h"

#ifdef _MPI
#include "mpi_csr_matrix.h"
#include "mpi_vector.h"
#endif  // #ifdef _MPI
#include BS_STOP_PLUGIN_IMPORT ()



namespace blue_sky
{
namespace bos_helper
{
namespace helper
{
namespace tmp_ns
{
      template <typename T> struct is_int
        {
          enum { value = 0};
        };
      template <>						struct is_int <int>
        {
          enum { value = 1};
        };
      template <>						struct is_int <unsigned int>
        {
          enum { value = 1};
        };
      template <>						struct is_int <long>
        {
          enum { value = 1};
        };
      template <>						struct is_int <unsigned long>
        {
          enum { value = 1};
        };

template <class vector_v1_t, class vector_v2_t> inline typename vector_v1_t::value_type
    mv_vector_inner_product (const vector_v1_t &v1, const vector_v2_t &v2, int /* obsolete */ = 0)
    {
      BOOST_STATIC_ASSERT (helper::is_int <typename vector_v1_t::value_type>::value == 0);
      BOOST_STATIC_ASSERT (helper::is_int <typename vector_v2_t::value_type>::value == 0);

      typename vector_v1_t::value_type sum = 0;
      size_t i = 0;
      size_t n = v1.size ();

      BS_ASSERT (v1.size () == v2.size ());

#ifdef MV_VECTOR_INNER_PRODUCT_PARALLEL
#pragma omp parallel for reduction (+: sum)
#endif //MV_VECTOR_INNER_PRODUCT_PARALLEL
      for (i = 0; i < n; ++i)
        {
          sum += v1[i] * v2[i];
        }

      return sum;
    }
}
}
}
}

namespace blue_sky
  {
 //namespace tmp_ns
 //{

  //////////////////////////////////////////////////////////////////////////
  //  linear_solver_prop

  //! constructor
  linear_solver_prop::linear_solver_prop (bs_type_ctor_param /*param*/)
  {
    this->params_names.resize(LS_TOTAL);
    EREG(linear_solver_prop,FP_RELATIVE_FACTOR,"relative factor");
    EREG(linear_solver_prop,FP_TOLERANCE,"tolerance");
    EREG(linear_solver_prop,FP_MATBAL_TOLERANCE,"matbal tolerance");
    EREG(linear_solver_prop,FP_FINAL_RESID,"final resid");

    EREG(linear_solver_prop,I_ITERS,"iterations count");
    EREG(linear_solver_prop,I_MAX_ITERS,"max iterations count");
    EREG(linear_solver_prop,I_SUCCESS,"success");

    resize (LS_TOTAL);

    set_max_iters (30);
    set_tolerance (1.0e-5);
    set_relative_factor (1.0);
    set_success (0);
    set_iters (0);
    set_final_resid (0);


  }

  //! copy constructor
  linear_solver_prop::linear_solver_prop (const linear_solver_prop &prop)
        : bs_refcounter (), named_pbase ()
  {
    if (&prop != this)
      *this = prop;
  }

  //! destructor
  linear_solver_prop::~linear_solver_prop ()
  {
  }


  /**
  * @brief set maximum number of iterations
  *
  * @param n_iters -- number of iterations
  *
  * @return 0 if success
  */
  int
  linear_solver_prop::set_max_iters (int n_iters)
  {
    int r_code = 0;

    if (n_iters < 0)
      n_iters = 1;

    //FI_DOUBLE_ARRAY_REALLOCATOR (resid, n_iters + 2, r_code);
    //FI_DOUBLE_ARRAY_REALLOCATOR (convergence_rate, n_iters + 2, r_code);

    set_param(I_MAX_ITERS,n_iters);
    return r_code;
  }

  /**
  * @brief get parameter name
  *
  * @param idx - parameter name
  *
  * @return parameter name
  */
  const std::string & linear_solver_prop::get_params_name (idx_type idx)
  {
    return params_names[idx].name;
  }

  //////////////////////////////////////////////////////////////////////////
  //  linear_solver_base

  //! constructor
  template <class strategy_t>
  linear_solver_base<strategy_t>::linear_solver_base (bs_type_ctor_param /*param*/, bs_node::sp_node node)
  :bs_node (node)
  {
    init();
  }

  template <class strategy_t>
  linear_solver_base<strategy_t>::linear_solver_base (bs_type_ctor_param)
  :  bs_node(bs_node::create_node (new solver_trait ()))
  {
    init();
  }

  template <class strategy_t>
  void linear_solver_base<strategy_t>::init()
  {
    prop = BS_KERNEL.create_object (linear_solver_prop::bs_type ());
    BS_ASSERT (prop);

    //sp_link prec_link = bs_link::create (bs_node::create_node (), "prec");
    //sp_link prop_link = bs_link::create (bs_node::create_node (), "prop");

    //bs_node::insert (prec_link, false);
    //bs_node::insert (prop_link, false);
  }

  template <class strategy_t>
  linear_solver_base<strategy_t>::linear_solver_base (const linear_solver_base &solver)
        : bs_refcounter (), bs_node ()
  {
    if (&solver != this)
      *this = solver;
  }

  //! destructor
  template <class strategy_t>
  linear_solver_base<strategy_t>::~linear_solver_base ()
  {
    //prec = 0;
    //prop = 0;
  }

  //////////////////////////////////////////////////////////////////////////
  //  gmres2_solver

  //! constructor
  template <class strat_t>
  gmres_solver2<strat_t>::gmres_solver2 (bs_type_ctor_param param)
  : linear_solver_base<strat_t> (param)
  {
    m = 30;
  }

  //! copy constructor
  template <class strat_t>
  gmres_solver2<strat_t>::gmres_solver2(const gmres_solver2 &solver)
        : bs_refcounter (), linear_solver_base<strat_t> (solver)
  {
    if (&solver != this)
      *this = solver;
  }

  //! destructor
  template <class strat_t>
  gmres_solver2<strat_t>::~gmres_solver2 ()
  {
    m = 30;
  }


   template <class strat_t>
   int gmres_solver2<strat_t>::solve(matrix_t *matrix, rhs_item_array_t &rhs, item_array_t &solution)
   {
     return templ_solve (matrix, rhs, solution);
   }

   template <class strat_t>
   int gmres_solver2<strat_t>::solve_prec(matrix_t *matrix, item_array_t &rhs, item_array_t &solution)
   {
     return templ_solve (matrix, rhs, solution);
   }


  /**
  * @brief restarted GMRES linear solver
  *
  * @param matrix -- matrix
  * @param rhs -- right hand side
  * @param solution -- solution
  *
  * @return 0 if success
  */
  template <class strat_t> template <class rhs_t> int
  gmres_solver2<strat_t>::templ_solve (matrix_t *matrix, rhs_t &rhs, item_array_t &solution)
  {
    //tools::save_seq_vector ("rhs_new").save (rhs);
    typedef typename strat_t::item_t fp_type;

    BS_ASSERT (matrix);
    BS_ASSERT (rhs.size ());
    BS_ASSERT (solution.size ());
    BS_ASSERT (rhs.size () == solution.size ()) (rhs.size ()) (solution.size ());
    BS_ASSERT (base_t::prop);

    const smart_ptr<linear_solver_prop > &lprop (base_t::prop);

    fp_type *rs, *hh, *c, *s;
    fp_type gamma, t, r_norm, b_norm, den_norm;
    fp_type *cur_h;
    const double epsmac = 1.e-16;

    // do barrier at mpi or nothing do at seq (base)
    // in future we can hold instance of barrier_t as a data member
    barrier_t ().barrier (barrier_t::comm_world_v ());

    index_t n = matrix->n_rows * matrix->n_block_size;
    item_t tol = lprop->get_tolerance ();

    // check workspace
    assign (base_t::wksp, n * (m + 3) + (m + 2) * (m + 1) + 2 * m, 0);

    vec_p.assign (m + 1, item_array_t ());
    for (int i = 0, cnt = (int)vec_p.size (); i < cnt; ++i)
      {
        vec_p[i] = item_array_t ();
        vec_p[i].template init_by_matrix<matrix_t> (matrix);
      }

    vec_w.template init_by_matrix<matrix_t> (matrix);
    vec_r.template init_by_matrix<matrix_t> (matrix);

    // initialize work arrays
    s = &base_t::wksp[0];
    c = s + m;
    rs = c + m;
    hh = rs + m + 1;

    lprop->set_success (0);

    assign (solution, n, 0);

    //calculate_init_r (n, matrix, solution, rhs, p);
    matrix->calc_lin_comb (-1.0, 1.0, solution, rhs, vec_p[0]);

    b_norm = sqrt (bos_helper::helper::tmp_ns::mv_vector_inner_product (rhs, rhs));
    //BOSOUT (section::solvers, level::low) << "b_norm = " << b_norm << bs_end;


    /* Since it is does not diminish performance, attempt to return an error flag
    and notify users when they supply bad input. */

    r_norm = sqrt (bos_helper::helper::tmp_ns::mv_vector_inner_product (vec_p[0], vec_p[0]));

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
    lprop->set_relative_factor (den_norm);

    index_t iter = 0;
    index_t max_iter = lprop->get_max_iters ();
    for (iter = 0; iter < max_iter;)
      {
        /* initialize first term of hessenberg system */
        rs[0] = r_norm;

        if (r_norm < epsmac || r_norm <= tol)
          break;

        t = 1.0 / r_norm;

        scale_vector (vec_p[0], t);

        index_t i = 0;
        for (i = 1; i < m && iter < max_iter; ++i, ++iter)
          {
            barrier_t ().barrier (barrier_t::comm_world_v ());

            //memset (r, 0, sizeof (double) * n);
            if (this->prec)
              {
                if (this->prec->solve_prec (matrix, vec_p[i - 1], vec_r))
                  {
                    bs_throw_exception ("GMRES: Preconditioner failed");
                  }
              }
            else // no precondition (preconditioner=identity_matrix)
              {
                vec_r.assign (vec_p[i - 1].begin (), vec_p[i - 1].end ());
              }

            barrier_t ().barrier (barrier_t::comm_world_v ());

            assign (vec_p[i], n, 0);
            matrix->matrix_vector_product (vec_r, vec_p[i]);

            /* modified Gram_Schmidt */
            cur_h = hh + (i - 1) * (m + 1);
            for (index_t j = 0; j < i; ++j)
              {
                cur_h[j] = bos_helper::helper::tmp_ns::mv_vector_inner_product (vec_p[j], vec_p[i], n);
                t = -cur_h[j];

                axpy (vec_p[i], vec_p[j], t);
              }
            t = sqrt (bos_helper::helper::tmp_ns::mv_vector_inner_product (vec_p[i], vec_p[i], n));
            cur_h[i] = t;
            if (t > epsmac)
              {
                t = 1.0 / t;

                scale_vector (vec_p[i], t);

              }
            /* done with modified Gram_schmidt and Arnoldi step.
            update factorization of hh */
            for (index_t j = 1; j < i; j++)
              {
                t = cur_h[j - 1];
                cur_h[j - 1] = c[j - 1] * t + s[j - 1] * cur_h[j];
                cur_h[j] = -s[j - 1] * t + c[j - 1] * cur_h[j];
              }

            gamma = sqrt (cur_h[i - 1] * cur_h[i - 1] + cur_h[i] * cur_h[i]);

            if (gamma <epsmac)
              gamma = epsmac;

            c[i - 1]  = cur_h[i - 1] / gamma;
            s[i - 1]  = cur_h[i]     / gamma;

            rs[i]     = -s[i - 1] * rs[i - 1];
            rs[i - 1] =  c[i - 1] * rs[i - 1];

            /* determine residual norm */
            cur_h[i - 1] = c[i - 1] * cur_h[i - 1] + s[i - 1] * cur_h[i];

            r_norm = fabs (rs[i]);
            if (r_norm <= tol)
              break;

          }

        if (i == m || iter == max_iter)
          {
            i = i - 1;
          }

        /* now compute solution, first solve upper triangular system */
        rs[i - 1] = rs[i - 1] / hh[i - 1 + (i - 1) * (m + 1)];
        for (index_t k = i - 2; k >= 0; --k)
          {
            t = rs[k];
            for (index_t j = k + 1; j < i; j++)
              {
                t -= hh[k + j * (m + 1)] * rs[j];
              }
            rs[k] = t / hh[k + k * (m + 1)];
          }

        /* form linear combination of p's to get solution */
        //vec_w.assign (vec_p[0].begin (), vec_p[0].end ());
        memcpy (&vec_w[0], &vec_p[0][0], sizeof (item_t) * n);
        t = rs[0];

        scale_vector (vec_w, t);

        for (index_t j = 1; j < i; j++)
          {
            t = rs[j];
            axpy (vec_w, vec_p[j], t);
          }

        //memset (r, 0, sizeof (double) * n);

        barrier_t ().barrier (barrier_t::comm_world_v ());

        if (this->prec)
          {
            if (base_t::prec->solve_prec (matrix, vec_w, vec_r))
              {
                bs_throw_exception ("GMRES: Preconditioner failed");
              }
          }
        else // no precondition (preconditioner=identity_matrix)
          {
            //memcpy (r, w, n * sizeof (fp_type));
            vec_r.assign (vec_w.begin (), vec_w.end ());
          }

        barrier_t ().barrier (barrier_t::comm_world_v ());

        axpy (solution, vec_r, 1.0);

        /* check for convergence, evaluate actual residual */
        if (r_norm <= tol)
          {
            //break;
#if 1
            //calculate_init_r (n, matrix, solution, rhs, r);
            matrix->calc_lin_comb (-1.0, 1.0, solution, rhs, vec_r);

            r_norm = sqrt (bos_helper::helper::tmp_ns::mv_vector_inner_product (vec_r, vec_r, n));

            if (r_norm <= tol)
              break;
            else
              {
                //vec_p[0].assign (vec_r.begin (), vec_r.end ());
                memcpy (&vec_p[0][0], &vec_r[0], sizeof (item_t) * n);
                i = 0;
                ++iter;
              }
#endif //0
          }

        /* compute residual vector and continue loop */

        for (index_t j = i; j > 0; --j)
          {
            rs[j - 1] = -s[j - 1] * rs[j];
            rs[j] = c[j - 1] * rs[j];
          }

        if (i)
          {
            t = rs[0];
            scale_vector (vec_p[0], t);

          }
        for (index_t j = 1; j < i + 1; ++j)
          {
            t = rs[j];
            axpy (vec_p[0], vec_p[j], t);
          }
      }

    lprop->set_iters (iter + 1);
    lprop->set_success (1);

    if (den_norm > 1.0e-12)
      lprop->set_final_resid (r_norm / den_norm);
    else
      lprop->set_final_resid (r_norm);

#ifdef _DEBUG
    BOSOUT (section::solvers, level::low) << "r_norm = " << r_norm << " r_norm / den_norm = " << r_norm / den_norm << " iter = " << (iter + 1) << bs_end;
#endif

    barrier_t ().barrier (barrier_t::comm_world_v ());

    return 0;
  }

  /**
  * @brief setup GMRES
  *
  * @param matrix -- input matrix
  *
  * @return 0 if success
  */
  template <class strat_t> int
  gmres_solver2<strat_t>::setup (matrix_t *matrix)
  {
    BS_ASSERT (matrix);
    if (!matrix)
      {
        bs_throw_exception ("GMRES: Passed matrix is null");
      }

    if (base_t::prec)
      {
        return base_t::prec->setup (matrix);
      }

    return 0;
  }


  //////////////////////////////////////////////////////////////////////////
  //  bicgstab_solver

  //! constructor
  template <class strat_t>
  bicgstab_solver<strat_t>::bicgstab_solver (bs_type_ctor_param param)
      : linear_solver_base<strat_t> (param)
  {}

  //! copy constructor
  template <class strat_t>
  bicgstab_solver<strat_t>::bicgstab_solver(const bicgstab_solver &solver)
      : bs_refcounter (), linear_solver_base<strat_t> (solver)
  {
    if (&solver != this)
      *this = solver;
  }

  //! destructor
  template <class strat_t>
  bicgstab_solver<strat_t>::~bicgstab_solver ()
  {}



  template <class strat_t>
  int bicgstab_solver<strat_t>::solve(matrix_t *matrix, rhs_item_array_t &rhs, item_array_t &solution)
  {
    return templ_solve (matrix, rhs, solution);
  }

  template <class strat_t>
  int bicgstab_solver<strat_t>::solve_prec(matrix_t *matrix, item_array_t &rhs, item_array_t &solution)
  {
     return templ_solve (matrix, rhs, solution);
  }

  /*!
  * \brief BiCGStab linear solver
  *
  * \param[in] matrix -- pointer to the matrix
  * \param[in] rhs -- right hand side
  * \param[out] solution -- solution
  *
  * \return 0 if success
  */
  template <class strat_t> template <class rhs_t> int
  bicgstab_solver<strat_t>::templ_solve (matrix_t *matrix, rhs_t &rhs, item_array_t &solution)
  {
    typedef item_t fp_type;

    BS_ERROR (matrix, "bicgstab_solve");
    BS_ERROR (rhs.size (), "bicgstab_solve");
    BS_ERROR (solution.size (), "bicgstab_solve");
    BS_ERROR (base_t::prop, "bicgstab_solve");

    const smart_ptr<linear_solver_prop> &lprop(this->prop);

    fp_type rho_1, rho_2 = 1, alpha = 1, beta, omega = 1;
    int iter;
    const double epsmac = 1e-24;
    fp_type r_norm, b_norm, den_norm, s_norm;
    //fp_type *x = solution;

    //OMP_TIME_MEASURE_START (bicgstab_solve_timer);

    item_t tol = this->prop->get_tolerance ();
    tol *= tol;
    //resid = prop->get_residuals ();
    //convergence_rate = prop->get_convergence_rate ();

    index_t max_iter  = this->prop->get_max_iters ();
    index_t n         = matrix->n_rows * matrix->n_block_size;

    item_array_t p (n);
    item_array_t phat (n);
    item_array_t s (n);
    item_array_t shat (n);
    item_array_t t (n), v (n), r (n), rtilde (n);

    lprop->set_success (0);

    assign (solution, n, 0);
    assign (r, n, 0);

    matrix->calc_lin_comb (-1.0, 1.0, solution, rhs, r);
    r_norm = bos_helper::mv_vector_inner_product (r, r);
    if (r_norm <= tol) // initial guess quite good
      return 0;

    rho_1 = r_norm;
    b_norm = sqrt (bos_helper::mv_vector_inner_product (rhs, rhs));

    p.assign      (r.begin (), r.end ());
    rtilde.assign (r.begin (), r.end ());
    assign (v, n, 0);

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
        //printf ("BiCGStab iteration: %d, resid = %le\n", iter, r_norm);
        //fflush (stdout);

        if (iter)
          {
            rho_1 = bos_helper::mv_vector_inner_product (r, rtilde); //in first iter equals to r_norm
            if (rho_1 == 0) // failure
              {
                if (den_norm > epsmac)
                  lprop->set_final_resid (r_norm / den_norm);
                else
                  lprop->set_final_resid (r_norm);

                bs_throw_exception (boost::format ("BICGSTAB: Failure - rho_1 == 0, resid = %le") % lprop->get_final_resid ());
              }
            beta = (rho_1 / rho_2) * (alpha / omega);
            // p = r + beta * (p - omega * v);
            //AXPY_AYPX (p, -omega, v, beta, r, k, n);
            axpy_aypx (p, -omega, v, beta, r);
          }

        // phat = M^(-1) * p;
        if (this->prec)
          {
            if (base_t::prec->solve_prec (matrix, p, phat))
              {
                bs_throw_exception ("BICGSTAB: Preconditioner failed");
              }
          }
        else // no precondition (preconditioner=identity_matrix)
          {
            phat.assign (p.begin (), p.end ());
          }

        // v = A * phat;
        assign (v, n, 0);
        matrix->matrix_vector_product (phat, v);

        alpha = bos_helper::mv_vector_inner_product (rtilde, v);
        if (alpha > epsmac || alpha < -epsmac)
          alpha = rho_1 / alpha;
        else // failure
          {
            if (den_norm > epsmac)
              lprop->set_final_resid (r_norm / den_norm);
            else
              lprop->set_final_resid (r_norm);

            bs_throw_exception (boost::format ("BICGSTAB: Failure - (rtilde, v) == 0, resid = %le") % lprop->get_final_resid ());
          }

        // s = r - alpha * v;
        s.assign (r.begin (), r.end ());
        axpy (s, v, -alpha);
        //AXPY (s, -alpha, v, k, n);

        //x = x + alpha * phat;
        //AXPY (x, alpha, phat, k, n);
        //axpy (x, phat, alpha);
        axpy (solution, phat, alpha);

        s_norm = bos_helper::mv_vector_inner_product (s, s);
        if (s_norm < tol)
          {
            //check convergence
            //matrix->calc_lin_comb (-1.0, 1.0, x, rhs, t);// t is buffer
            matrix->calc_lin_comb (-1.0, 1.0, solution, rhs, t);// t is buffer
            r_norm = bos_helper::mv_vector_inner_product (t, t);
            if (r_norm <= tol)
              break;
          }

        // shat = M^(-1) * s;
        if (this->prec)
          {
            if (base_t::prec->solve_prec (matrix, s, shat))
              {
                bs_throw_exception ("BICGSTAB: Preconditioner failed");
              }
          }
        else // no precondition (preconditioner=identity_matrix)
          {
            shat.assign (s.begin (), s.end ());
          }

        // t = A * shat;
        assign (t, n, 0);
        matrix->matrix_vector_product (shat, t);

        // omega = (t,s) / (t,t);
        omega = bos_helper::mv_vector_inner_product (t, t);
        if (omega > epsmac)
          {
            omega = bos_helper::mv_vector_inner_product (t, s) / omega;
          }
        else // failure
          {
            if (den_norm > epsmac)
              lprop->set_final_resid (r_norm / den_norm);
            else
              lprop->set_final_resid (r_norm);

            bs_throw_exception (boost::format ("BICGSTAB: Failure - (t, t) == 0, resid = %le") % lprop->get_final_resid ());
          }

        if (omega < epsmac) // failure
          {
            if (den_norm > epsmac)
              lprop->set_final_resid (r_norm / den_norm);
            else
              lprop->set_final_resid (r_norm);

            bs_throw_exception (boost::format ("BICGSTAB: Failure - omega == 0, resid = %le") % lprop->get_final_resid ());
          }

        //x = x + omega * shat;
        //AXPY (x, omega, shat, k, n);
        //axpy (x, shat, omega);
        axpy (solution, shat, omega);

        //r = s - omega * t;
        //memcpy (r, s, n * sizeof (fp_type));
        //AXPY (r, -omega, t, k, n);
        r.assign (s.begin (), s.end ());
        axpy (r, t, -omega);
        /*
        //additional check convergence
        mv_calc_lin_comb (matrix, -1.0, 1.0, x, rhs, s); // s is buffer
        r_norm = mv_vector_inner_product (s, s, n);
        */
        r_norm = bos_helper::mv_vector_inner_product (r, r);
        if (r_norm <= tol)
          break;

        rho_2 = rho_1;
      } // end of main loop

    lprop->set_iters (iter + 1);
    lprop->set_success (1);

    //printf ("BiCGStab after iteration: %d, resid = %le\n", iter, r_norm);
    /*
    //additional checking convergence
    mv_calc_lin_comb (matrix, -1.0, 1.0, solution, rhs, r);
    r_norm = mv_vector_inner_product (r, r, n);
    */
    if (den_norm > epsmac)
      lprop->set_final_resid (r_norm / den_norm);
    else
      lprop->set_final_resid (r_norm);

    BOSOUT (section::solvers, level::low) << "r_norm = " << r_norm << " r_norm / den_norm = " << r_norm / den_norm << " iter = " << iter << bs_end;

    // printf ("BiCGStab OK! iters = %d, resid = %le\n", prop->iters, prop->final_resid);
    //OMP_TIME_MEASURE_END (bicgstab_solve_timer);

    return 0;
  }

  /**
  * @brief setup for BiCGStab
  *
  * @param matrix -- input matrix
  *
  * @return 0 if success
  */
  template <class strat_t> int
  bicgstab_solver<strat_t>::setup (matrix_t *matrix)
  {
    BS_ASSERT (matrix);
    if (!matrix)
      {
        bs_throw_exception ("BICGSTAB: Passed matrix is null");
      }

    BS_ASSERT (base_t::prop);
    if (base_t::prec)
      {
        return base_t::prec->setup (matrix);
      }

    return 0;
  }

  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE(linear_solver_prop);
  BLUE_SKY_TYPE_STD_COPY(linear_solver_prop);
  BLUE_SKY_TYPE_IMPL(linear_solver_prop, objbase, "linear_solver_prop", "Property for linear solvers", "Property for linear solvers");

  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE_T_DEF(linear_solver_base, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF(linear_solver_base, (class));

  BLUE_SKY_TYPE_IMPL_T_EXT(1, (linear_solver_base<base_strategy_fi>) , 1, (objbase), "linear_solver_base_fi", "linear solver", "linear solver", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (linear_solver_base<base_strategy_di>) , 1, (objbase), "linear_solver_base_di", "linear solver", "linear solver", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (linear_solver_base<base_strategy_mixi>) , 1, (objbase), "linear_solver_base_mixi", "linear solver", "linear solver", false);

#ifdef _MPI
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (linear_solver_base<mpi_strategy_di>) , 1, (objbase), "linear_solver_mpi_di", "linear solver", "linear solver", false);
#endif
  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE_T_DEF(gmres_solver2, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF(gmres_solver2, (class));

  BLUE_SKY_TYPE_IMPL_T_EXT(1, (gmres_solver2<base_strategy_fi>) , 1, (linear_solver_base<base_strategy_fi>), "gmres_solver2_base_fi", "GMRES linear solver", "GMRES linear solver", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (gmres_solver2<base_strategy_di>) , 1, (linear_solver_base<base_strategy_di>), "gmres_solver2_base_di", "GMRES linear solver", "GMRES linear solver", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (gmres_solver2<base_strategy_mixi>) , 1, (linear_solver_base<base_strategy_mixi>), "gmres_solver2_base_mixi", "GMRES linear solver", "GMRES linear solver", false);

#ifdef _MPI
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (gmres_solver2<mpi_strategy_di>) , 1, (linear_solver_base<mpi_strategy_di>), "gmres_solver2_mpi_di", "GMRES linear solver", "GMRES linear solver", false);
#endif
  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE_T_DEF(bicgstab_solver, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF(bicgstab_solver, (class));

  BLUE_SKY_TYPE_IMPL_T_EXT(1, (bicgstab_solver<base_strategy_fi>) , 1, (linear_solver_base<base_strategy_fi>), "bicgstab_solver_base_fi", "BiCG linear solver", "BiCG linear solver", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (bicgstab_solver<base_strategy_di>) , 1, (linear_solver_base<base_strategy_di>), "bicgstab_solver_base_di", "BiCG linear solver", "BiCG linear solver", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (bicgstab_solver<base_strategy_mixi>) , 1, (linear_solver_base<base_strategy_mixi>), "bicgstab_solver_base_mixi", "BiCG linear solver", "BiCG linear solver", false);

  //////////////////////////////////////////////////////////////////////////
  // register types

  bool
  linear_solver_prop_register_type (const blue_sky::plugin_descriptor &pd)
  {
    bool res = true;
    res &= BS_KERNEL.register_type (pd, linear_solver_prop::bs_type ());

    return res;
  }
  bool
  linear_solvers_register_type (const blue_sky::plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, linear_solver_base<base_strategy_di>::bs_type ());
    BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, linear_solver_base<base_strategy_fi>::bs_type ());
    BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, linear_solver_base<base_strategy_mixi>::bs_type ());
    BS_ASSERT (res);


    res &= BS_KERNEL.register_type (pd, gmres_solver2<base_strategy_di>::bs_type ());
    BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, gmres_solver2<base_strategy_fi>::bs_type ());
    BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, gmres_solver2<base_strategy_mixi>::bs_type ());
    BS_ASSERT (res);



#ifdef _MPI
    res &= BS_KERNEL.register_type (pd, gmres_solver2<mpi_strategy_di>::bs_type ());
#endif

    res &= BS_KERNEL.register_type (pd, bicgstab_solver<base_strategy_di>::bs_type ());
    BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, bicgstab_solver<base_strategy_fi>::bs_type ());
    BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, bicgstab_solver<base_strategy_mixi>::bs_type ());
    BS_ASSERT (res);



    return res;
  }

// } // namespace tmp_ns
} // namespace blue_sky
