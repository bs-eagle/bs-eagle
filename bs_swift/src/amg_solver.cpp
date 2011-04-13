/**
 * @file amg_solver.cpp
 * @brief implementation of AMG solver
 * @author
 * @version
 * @date 2011-03-30
 */

#include "amg_solver.h"

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
      smbuilder.reserve (AMG_N_LEVELS_RESERVE);
      coarser.reserve (AMG_N_LEVELS_RESERVE);
      pbuilder.reserve (AMG_N_LEVELS_RESERVE);

      // set default
      smbuilder.resize (1);
      coarser.resize (1);
      pbuilder.resize (1);
      pre_smoother.resize (1);
      post_smoother.resize (1);
      smbuilder[0] = BS_KERNEL.create_object ("simple_smbuilder");
      coarser[0] = BS_KERNEL.create_object ("cljp_coarse");
      pbuilder[0] = BS_KERNEL.create_object ("standart2_pbuild");
      pre_smoother[0] = BS_KERNEL.create_object ("gs_solver");
      post_smoother[0] = BS_KERNEL.create_object ("gs_solver");

      lu_solver = BS_KERNEL.create_object ("blu_solver");
      lu_fact = BS_KERNEL.create_object ("dens_matrix");
      wksp = BS_KERNEL.create_object (v_double::bs_type ());
      BS_ASSERT (lu_fact);
      BS_ASSERT (lu_solver);
      BS_ASSERT (wksp);
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
        prop->add_property_i (1, n_levels_idx,
                              std::string ("Number of coarse levels"));
      }

    int amg_solver::solve (sp_matrix_t matrix_, spv_double sp_rhs, spv_double sp_sol)
    {
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
      const t_long n_levels = get_n_levels ();
      std::cout<<"n_levels = "<<n_levels<<" a.size = "<<a.size ()<<" p.size = "<<p.size ()<<"\n";

      // rhs and solution on first level
      rhs.push_back (sp_rhs);
      sol.push_back (sp_sol);

      int level = 0;
      for (level = 0; level < n_levels; ++level)
        {
          sp_smooth_t pre_smoother = get_pre_smoother (level);
          t_long n = a[level]->get_n_rows ();
          sol[level]->resize (n);
          rhs[level]->resize (n);
          wksp->resize (n);
          std::cout<<"AMG solve level = "<<level<<" n_rows = "<<n<<"\n";

          //pre-smooth
          sp_prop_t smoother_prop;
          smoother_prop = pre_smoother->get_prop ();
          smoother_prop->set_b ("inverse", false);

          // smooth C-points
          smoother_prop->set_i ("cf_type", 1);
          pre_smoother->smooth (a[level], cf[level], get_n_pre_smooth_iters (),
                                rhs[level], sol[level]);
          // smooth F-points
          smoother_prop->set_i ("cf_type", -1);
          pre_smoother->smooth (a[level], cf[level], get_n_pre_smooth_iters (),
                                rhs[level], sol[level]);

          // calculate r = b - Ax
          if (a[level]->calc_lin_comb (-1.0, 1.0, sol[level], rhs[level], wksp))
            return -1;

          //set rhs[level + 1]=0
          std::fill (rhs[level + 1]->begin (), rhs[level + 1]->begin (), 0);
          // restriction: b^(k+1) = P^T * r
          if (p[level]->matrix_vector_product_t (wksp, rhs[level + 1]))
            return -1;
        }

      std::cout<<"LU solve...";
      lu_solver->solve (lu_fact, rhs[level], sol[level]);
      std::cout<<"OK\n";

      for (level = n_levels - 1; level >= 0; --level)
        {
          sp_smooth_t post_smoother = get_post_smoother (level);
          t_long n = a[level]->get_n_rows ();
          std::cout<<"AMG solve level = "<<level<<" n_rows = "<<n<<"\n";

          //set sol[level + 1]=0
          std::fill (sol[level + 1]->begin (), sol[level + 1]->begin (), 0);
          //interpolation x = x + P * e^(k+1)
          if (p[level]->calc_lin_comb (1.0, 1.0, sol[level + 1], sol[level], sol[level]))
            return -6;

          //post-smooth
          sp_prop_t smoother_prop;
          smoother_prop = post_smoother->get_prop ();
          smoother_prop->set_b ("inverse", false);

          // smooth F-points
          smoother_prop->set_i ("cf_type", -1);
          post_smoother->smooth (a[level], cf[level], get_n_post_smooth_iters (),
                                rhs[level], sol[level]);
          // smooth C-points
          smoother_prop->set_i ("cf_type", 1);
          post_smoother->smooth (a[level], cf[level], get_n_post_smooth_iters (),
                                rhs[level], sol[level]);
        }

      return 0;
    }

    int amg_solver::solve_prec (sp_matrix_t /*matrix*/, spv_double /*rhs*/, spv_double /*sol*/)
    {
      return 0;//prec->solve (matrix, rhs, sol);
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
      const t_double strength_threshold = get_strength_threshold ();
      const t_double max_row_sum = get_max_row_sum ();
      const t_long n_last_level_points = get_n_last_level_points ();
      const t_long max_connections = 0;//!
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

      t_double cop = matrix->get_n_non_zeros ();
      // matrix on first level
      a.push_back (matrix);

      int level;
      for (level = 0;;++level)
        {
          t_long n = matrix->get_n_rows ();
          std::cout<<"AMG setup level = "<<level<<" n_rows = "<<n<<"\n";

          if (n <= n_last_level_points)
            {
              break;
            }

          //init tools
          sp_smbuild_t s_builder = get_smbuilder (level);
          sp_coarse_t  coarser   = get_coarser   (level);
          sp_pbuild_t  p_builder = get_pbuilder  (level);
          BS_ASSERT (s_builder);
          BS_ASSERT (coarser);
          BS_ASSERT (p_builder);
          std::cout<<"strength type: "<<s_builder->py_str ()<<"\n";
          std::cout<<"coarse   type: "<<coarser->py_str ()<<"\n";
          std::cout<<"interp   type: "<<p_builder->py_str ()<<"\n";

          std::cout<<"strength...";
          // build strength matrix (fill s_markers)
          spv_long s_markers = BS_KERNEL.create_object (v_long::bs_type ());
          BS_ASSERT (s_markers);
          s.push_back (s_markers);

          s_builder->build (matrix, strength_threshold, max_row_sum, s_markers);
          std::cout<<"OK\ncoarse...";

          // coarse (fill cf_markers)
          spv_long cf_markers = BS_KERNEL.create_object (v_long::bs_type ());
          BS_ASSERT (cf_markers);
          cf.push_back (cf_markers);
          spv_double measure = BS_KERNEL.create_object (v_double::bs_type ());
          BS_ASSERT (measure);
          measure->resize (n);

          t_long n_coarse_size = coarser->build (matrix, measure, cf_markers, s_markers);
          if (n_coarse_size == n || n_coarse_size < 1)
            {
              std::cout<<"coarse failed: n = "<<n<<" n_coarse = "<<n_coarse_size<<"\n";
              break;
            }
          std::cout<<"OK\nbuild P...";
          // build prolongation (interpolation) matrix
          sp_bcsr_t p_matrix = BS_KERNEL.create_object ("bcsr_matrix");
          BS_ASSERT (p_matrix);
          p.push_back (p_matrix);

          p_builder->build (matrix, n_coarse_size, max_connections,
                            cf_markers, s_markers, p_matrix);
          std::cout<<"OK\ntriple...";
          // initialize next level matrix
          sp_bcsr_t a_matrix = BS_KERNEL.create_object ("bcsr_matrix");
          BS_ASSERT (a_matrix);
          a.push_back (a_matrix);

          sp_bcsr_t r_matrix = BS_KERNEL.create_object ("bcsr_matrix");
          BS_ASSERT (r_matrix);

          r_matrix->build_transpose (p_matrix, 0, 0, 0);

          a_matrix->triple_matrix_product (r_matrix, matrix, p_matrix, update);
          matrix = a_matrix;
          cop += matrix->get_n_non_zeros ();
          std::cout<<"OK\n";
        }

      set_n_levels (level);
      cop /= matrix->get_n_non_zeros ();
      set_cop (cop);

      std::cout<<"AMG setup OK. n_levels = "<<level<<" Cop = "<<cop<<"\n";

      std::cout<<"LU setup...";
      lu_fact->init_by_matrix (matrix);
      lu_solver->setup (lu_fact);
      std::cout<<"OK\n";
/*
      for (int level = 0; level < n_levels; ++level)
        {
          std::cout<<"AMG level = "<<level<<
          " S_n = "<<S[level]->get_n_rows ()<<
          " S_nnz = "<<S[level]->get_n_non_zeros ()<<"\n";
        }
*/
      return 0;
    }

    //////////////////////////////////////////////////////////////////////////
    BLUE_SKY_TYPE_STD_CREATE (amg_solver);
    BLUE_SKY_TYPE_STD_COPY (amg_solver);

    BLUE_SKY_TYPE_IMPL (amg_solver, amg_solver_iface, "amg_solver", "Algebraic Multigrid linear solver and preconditioner", "Algebraic Multigrid linear solver and preconditioner");

}  // blue_sky namespace
