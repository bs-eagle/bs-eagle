/**
 * @file amg_solver.cpp
 * @brief implementation of AMG solver
 * @author
 * @version
 * @date 2011-03-30
 */

#include "amg_solver.h"
#include "mv_functions.h"

#define AMG_N_LEVELS_RESERVE 10

namespace blue_sky
{
    //! constructor
    amg_solver::amg_solver (bs_type_ctor_param)
                : amg_solver_iface ()
    {
      prop = BS_KERNEL.create_object ("prop");
      if (!prop)
        {
          bs_throw_exception ("Type (prop) not registered");
        }
      init_prop ();

      //
      a.reserve (AMG_N_LEVELS_RESERVE);
      p.reserve (AMG_N_LEVELS_RESERVE);
      s.reserve (AMG_N_LEVELS_RESERVE);
      cf.reserve (AMG_N_LEVELS_RESERVE);
      rhs.reserve (AMG_N_LEVELS_RESERVE);
      sol.reserve (AMG_N_LEVELS_RESERVE);
      smbuilder_vec.reserve (AMG_N_LEVELS_RESERVE);
      coarser_vec.reserve (AMG_N_LEVELS_RESERVE);
      pbuilder_vec.reserve (AMG_N_LEVELS_RESERVE);
      pre_smoother_vec.reserve (AMG_N_LEVELS_RESERVE);
      post_smoother_vec.reserve (AMG_N_LEVELS_RESERVE);

      // set default
      smbuilder_vec.resize (1);
      coarser_vec.resize (1);
      pbuilder_vec.resize (1);
      pre_smoother_vec.resize (1);
      post_smoother_vec.resize (1);
      smbuilder_vec[0] = BS_KERNEL.create_object ("simple_smbuilder");
      coarser_vec[0] = BS_KERNEL.create_object ("pmis2_coarse");
      pbuilder_vec[0] = BS_KERNEL.create_object ("standart2_pbuild");
      pre_smoother_vec[0] = BS_KERNEL.create_object ("gs_solver");
      post_smoother_vec[0] = BS_KERNEL.create_object ("gs_solver");
      //set smoother prop
      pre_smoother_vec[0]->get_prop ()->set_b ("inverse", false);
      post_smoother_vec[0]->get_prop ()->set_b ("inverse", false);

      lu_solver = BS_KERNEL.create_object ("blu_solver");
      lu_fact = BS_KERNEL.create_object ("dens_matrix");
      wksp = BS_KERNEL.create_object (v_double::bs_type ());
      r_matrix = BS_KERNEL.create_object ("bcsr_matrix");
      BS_ASSERT (lu_solver);
      BS_ASSERT (lu_fact);
      BS_ASSERT (wksp);
      BS_ASSERT (r_matrix);
    }

    //! copy constructor
    amg_solver::amg_solver(const amg_solver &solver)
      : bs_refcounter (), amg_solver_iface ()
    {
      if (&solver != this)
        *this = solver;
    }

    //! destructor
    amg_solver::~amg_solver ()
    {}

    //! set solver's properties
    void amg_solver::set_prop (sp_prop_t prop_)
    {
      prop = prop_;
      init_prop ();
    }

    void
    amg_solver::init_prop ()
      {
        prop->add_property_f (1.0e-4, tol_idx,
                              std::string ("Target tolerance for linear solver"));
        prop->add_property_f (1, final_res_idx,
                              std::string ("Solution residual"));
        prop->add_property_i (200, max_iters_idx,
                              std::string ("Maximum allowed number of iterations"));
        prop->add_property_i (0,iters_idx,
                              std::string ("Total number of used solver iterations"));
        prop->add_property_b (false, success_idx,
                              std::string ("True if solver successfully convergent"));
        // AMG input props
        prop->add_property_f (0.75, strength_threshold_idx,
                              std::string ("Threshold for strength matrix"));
        prop->add_property_f (0.01, max_row_sum_idx,
                              std::string ("Row sum threshold for strength matrix"));
        prop->add_property_i (100, n_last_level_points_idx,
                              std::string ("Minimal number of points in coarse grid"));
        prop->add_property_i (1, n_pre_smooth_iters_idx,
                              std::string ("Number of pre-smooth iterations"));
        prop->add_property_i (1, n_post_smooth_iters_idx,
                              std::string ("Number of post-smooth iterations"));
        // AMG output props
        prop->add_property_f (0, cop_idx,
                              std::string ("Operator compexity"));
        prop->add_property_i (0, n_levels_idx,
                              std::string ("Number of coarse levels"));
        prop->add_property_i (0, n_levels_max_idx,
                              std::string ("Maximal number of coarse levels"));
      }

    int amg_solver::solve (sp_matrix_t matrix_, spv_double sp_rhs, spv_double sp_sol)
    {
      const double epsmac = 1.e-16;
      BS_ASSERT (matrix_);
      BS_ASSERT (prop);

      BS_ASSERT (sp_rhs->size ());
      BS_ASSERT (sp_sol->size ());
      BS_ASSERT (sp_rhs->size () == sp_sol->size ()) (sp_rhs->size ()) (sp_sol->size ());

      sp_bcsr_t matrix (matrix_, bs_dynamic_cast ());
      if (!matrix)
        {
          bs_throw_exception ("AMG solve: Passed matrix is not BCSR");
        }

      //init amg props
      const t_long n_levels = prop->get_i (n_levels_idx);

      // rhs and solution on first level
      rhs[0] = sp_rhs;
      sol[0] = sp_sol;

      //initial guess
      sp_sol->assign (0);

      int level = 0;
      t_long n = matrix->get_n_rows ();
      t_double *rhs_ptr = &(*sp_rhs)[0];

      // max_iters>1 -- used as solver
      // max_iters=1 -- used as preconditioner
      int max_iter = prop->get_i (max_iters_idx);
      t_double tol = prop->get_f (tol_idx);
      prop->set_b (success_idx, false);

      t_double resid = tol + 1.0;
      if (max_iter > 1)
        {
          // calculate initial norm
          wksp->resize (n);
          t_double *wksp_ptr = &(*wksp)[0];
          matrix->calc_lin_comb (-1.0, 1.0, sp_sol, sp_rhs, wksp);
          t_double b_norm = sqrt (mv_vector_inner_product_n (rhs_ptr, rhs_ptr, n));
          t_double r_norm = sqrt (mv_vector_inner_product_n (wksp_ptr, wksp_ptr, n));
          if (b_norm > epsmac)
            { // convergence criterion |r_i|/|b| <= accuracy if |b| > 0
              tol *= b_norm;
            }
          else if (r_norm > epsmac)
            { // convergence criterion |r_i|/|r0| <= accuracy if |b| = 0
              tol *= r_norm;
            }
          resid = r_norm;
          //std::cout<<"n = "<<n<<" b_norm = "<<b_norm<<" r_norm = "<<r_norm
          //         <<" resid = "<<resid<<" tol = "<<tol<<"\n";
        }

      int iter;
      for (iter = 0; iter < max_iter && resid > tol; ++iter)
        {
          for (level = 0; level < n_levels; ++level)
            {
              sp_smooth_t pre_smoother = get_tool (pre_smoother_vec, level);
              BS_ASSERT (pre_smoother);

              t_long n = a[level]->get_n_rows ();
              wksp->resize (n);

              //initial guess on [level + 1] is sol=0
              sol[level + 1]->assign (0);

              // smooth all points
              pre_smoother->get_prop ()->set_i ("cf_type", 0);
              pre_smoother->smooth (a[level], NULL, prop->get_i (n_pre_smooth_iters_idx),
                                    rhs[level], sol[level]);
/*
              // smooth C-points
              pre_smoother->get_prop ()->set_i ("cf_type", 1);
              pre_smoother->smooth (a[level], cf[level], get_n_pre_smooth_iters (),
                                    rhs[level], sol[level]);
              // smooth F-points
              pre_smoother->get_prop ()->set_i ("cf_type", -1);
              pre_smoother->smooth (a[level], cf[level], get_n_pre_smooth_iters (),
                                    rhs[level], sol[level]);
*/
              // calculate r = b - Ax
              if (a[level]->calc_lin_comb (-1.0, 1.0, sol[level], rhs[level], wksp))
                return -1;

              // restriction: b^(k+1) = P^T * r
              rhs[level + 1]->assign (0);
              if (p[level]->matrix_vector_product_t (wksp, rhs[level + 1]))
                return -1;
            }

          lu_solver->solve (lu_fact, rhs[level], sol[level]);

          for (level = n_levels - 1; level >= 0; --level)
            {
              sp_smooth_t post_smoother = get_tool (post_smoother_vec, level);
              BS_ASSERT (post_smoother);

              //interpolation x = x + P * e^(k+1)
              if (p[level]->matrix_vector_product (sol[level + 1], sol[level]))
                return -6;

              // smooth all points
              post_smoother->get_prop ()->set_i ("cf_type", 0);
              post_smoother->smooth (a[level], NULL, prop->get_i (n_post_smooth_iters_idx),
                                     rhs[level], sol[level]);
/*
              // smooth F-points
              post_smoother->get_prop ()->set_i ("cf_type", -1);
              post_smoother->smooth (a[level], cf[level], get_n_post_smooth_iters (),
                                     rhs[level], sol[level]);
              // smooth C-points
              post_smoother->get_prop ()->set_i ("cf_type", 1);
              post_smoother->smooth (a[level], cf[level], get_n_post_smooth_iters (),
                                     rhs[level], sol[level]);
*/
            }

          if (max_iter > 1)
            {
              // calculate residual norm
              t_long n = matrix->get_n_rows ();
              wksp->resize (n);
              t_double *wksp_ptr = &(*wksp)[0];
              if (matrix->calc_lin_comb (-1.0, 1.0, sp_sol, sp_rhs, wksp))
                return -2;
              resid = sqrt (mv_vector_inner_product_n (wksp_ptr, wksp_ptr, n));
              //std::cout<<"AMG Iteration "<<iter + 1<<" resid = "<<resid<<"\n";
            }
        }//iter loop

      prop->set_i (iters_idx, iter);
      prop->set_b (success_idx, true);
      prop->set_f (final_res_idx, resid);

      return 0;
    }

    int amg_solver::solve_prec (sp_matrix_t matrix, spv_double rhs, spv_double sol)
    {
      return solve (matrix, rhs, sol);
    }

    /**
    * @brief setup for AMG
    * @param matrix -- input matrix
    * @return 0 if success
    */
    int amg_solver::setup (sp_matrix_t matrix_)
    {
      BS_ASSERT (prop);
      BS_ASSERT (lu_solver);

      //init amg props
      const t_double strength_threshold = prop->get_f (strength_threshold_idx);
      const t_double max_row_sum = prop->get_f (max_row_sum_idx);
      const t_long n_last_level_points = prop->get_i (n_last_level_points_idx);
      const bool update = false;//!

      if (!matrix_)
        {
          bs_throw_exception ("AMG setup: Passed matrix is null");
        }

      sp_bcsr_t matrix (matrix_, bs_dynamic_cast ());
      if (!matrix)
        {
          bs_throw_exception ("AMG setup: Passed matrix is not BCSR");
        }

      if (matrix->get_n_block_size () != 1)
        {
          bs_throw_exception ("AMG setup: Passed matrix block size != 1");
        }
      if (matrix->get_n_rows () != matrix->get_n_cols ())
        {
          bs_throw_exception ("AMG setup: Passed matrix is not square");
        }
      if (matrix->get_n_non_zeros () < 1)
        {
          bs_throw_exception ("AMG setup: Passed matrix has no nnz elements");
        }

      t_double nnz = matrix->get_n_non_zeros ();
      t_double cop = nnz;

      if (a.empty ()) //if setup called first time
        {
          // matrix on first level
          a.push_back (matrix);
          // first level solution and rhs vectors ptr=null
          // (will be replaced at solve())
          sol.push_back (NULL);
          rhs.push_back (NULL);
        }
      else
        {
          a[0] = matrix;
        }

      int level;
      for (level = 0;;++level)
        {
          t_long n = a[level]->get_n_rows ();
          //std::cout<<"AMG setup level = "<<level<<" n_rows = "<<n<<
          //" nnz = "<<a[level]->get_n_non_zeros ()<<"\n";

          if (n <= n_last_level_points)
            {
              break;
            }

          // init tools
          sp_smbuild_t smbuilder = get_tool (smbuilder_vec, level);
          sp_coarse_t  coarser   = get_tool (coarser_vec, level);
          sp_pbuild_t  pbuilder = get_tool (pbuilder_vec, level);
          BS_ASSERT (smbuilder);
          BS_ASSERT (coarser);
          BS_ASSERT (pbuilder);
          //std::cout<<"strength type: "<<s_builder->py_str ()<<"\n";
          //std::cout<<"coarse   type: "<<coarser->py_str ()<<"\n";
          //std::cout<<"interp   type: "<<p_builder->py_str ()<<"\n";

          // build strength matrix (fill s_markers)
          if (level >= get_n_levels_max ())
            {
              spv_long s_markers = BS_KERNEL.create_object (v_long::bs_type ());
              BS_ASSERT (s_markers);
              s.push_back (s_markers);
              spv_long cf_markers = BS_KERNEL.create_object (v_long::bs_type ());
              BS_ASSERT (cf_markers);
              cf.push_back (cf_markers);
            }

          t_long max_conn = smbuilder->build (a[level], strength_threshold,
                                              max_row_sum, s[level]);
          // calc s_nnz
          //t_long s_nnz = 0;
          //for (t_long ii = 0; ii < s_markers->size (); ++ii)
          //  {
          //    if ((*s_markers)[ii] > 0)
          //      s_nnz++;
          //  }
          //std::cout<<"a_nnz = "<<s_markers->size()<<" s_nnz = "<<s_nnz<<"\n";

          // coarse (fill cf_markers)
          cf[level]->resize (n);
          wksp->resize (n);//buffer for measure array
          t_long n_coarse_size = coarser->build (a[level], wksp, cf[level], s[level]);

          if (n_coarse_size == n || n_coarse_size < 1)
            {
              std::cout<<"coarse failed: n = "<<n<<" n_coarse = "<<n_coarse_size<<"\n";
              break;
            }

          // create objects
          if (level + 1 >= get_n_levels_max ())
            {
              // initialize prolongation matrix
              sp_bcsr_t p_matrix = BS_KERNEL.create_object ("bcsr_matrix");
              BS_ASSERT (p_matrix);
              p.push_back (p_matrix);
              // initialize next level matrix
              sp_bcsr_t a_matrix = BS_KERNEL.create_object ("bcsr_matrix");
              BS_ASSERT (a_matrix);
              a.push_back (a_matrix);
              // initialize next level solution and rhs vectors
              spv_double next_level_sol = BS_KERNEL.create_object (v_double::bs_type ());
              spv_double next_level_rhs = BS_KERNEL.create_object (v_double::bs_type ());
              BS_ASSERT (next_level_sol);
              BS_ASSERT (next_level_rhs);
              sol.push_back (next_level_sol);
              rhs.push_back (next_level_rhs);
            }

          // build prolongation (interpolation) matrix
          pbuilder->build (a[level], n_coarse_size, max_conn,
                            cf[level], s[level], p[level]);

          // initialize next level matrix
          r_matrix->build_transpose (p[level], 0, 0, 0);
          a[level + 1]->triple_matrix_product (r_matrix, a[level], p[level], update);

          sol[level + 1]->resize (n_coarse_size);
          rhs[level + 1]->resize (n_coarse_size);

          cop += a[level + 1]->get_n_non_zeros ();
        }

      set_n_levels (level);
      cop /= nnz;
      prop->set_f (cop_idx, cop);

      // init dense matrix from bcsr on last level
      // and build LU factorization for it
      lu_fact->init_by_matrix (a[level]);
      lu_solver->setup (lu_fact);
      return 0;
    }

    //////////////////////////////////////////////////////////////////////////
    BLUE_SKY_TYPE_STD_CREATE (amg_solver);
    BLUE_SKY_TYPE_STD_COPY (amg_solver);

    BLUE_SKY_TYPE_IMPL (amg_solver, amg_solver_iface, "amg_solver", "Algebraic Multigrid linear solver and preconditioner", "Algebraic Multigrid linear solver and preconditioner");

}  // blue_sky namespace
