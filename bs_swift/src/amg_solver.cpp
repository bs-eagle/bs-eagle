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
      //
      a.reserve (AMG_N_LEVELS_RESERVE);
      p.reserve (AMG_N_LEVELS_RESERVE);
      s.reserve (AMG_N_LEVELS_RESERVE);
      cf.reserve (AMG_N_LEVELS_RESERVE);
      smbuilder.reserve (AMG_N_LEVELS_RESERVE);
      coarser.reserve (AMG_N_LEVELS_RESERVE);
      pbuilder.reserve (AMG_N_LEVELS_RESERVE);

      // set default
      smbuilder.resize (1);
      coarser.resize (1);
      pbuilder.resize (1);
      smbuilder[0] = BS_KERNEL.create_object ("simple_smbuilder");
      coarser[0] = BS_KERNEL.create_object ("pmis2_coarse");
      pbuilder[0] = BS_KERNEL.create_object ("standart2_pbuild");

      lu_solver = BS_KERNEL.create_object ("blu_solver");
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
        // AMG
        prop->add_property_f (0.25, strength_threshold_idx,
                              std::string ("strength_threshold"));
        prop->add_property_f (0.01, max_row_sum_idx,
                              std::string ("max_row_sum"));
        prop->add_property_i (100,n_last_level_points_idx,
                              std::string ("n_last_level_points"));
      }

    int amg_solver::solve (sp_matrix_t matrix, spv_double sp_rhs, spv_double sp_sol)
    {
      BS_ASSERT (matrix);
      BS_ASSERT (sp_rhs->size ());
      BS_ASSERT (sp_sol->size ());
      BS_ASSERT (sp_rhs->size () == sp_sol->size ()) (sp_rhs->size ()) (sp_sol->size ());
      BS_ASSERT (prop);

      return 0;
    }

    int amg_solver::solve_prec (sp_matrix_t /*matrix*/, spv_double /*rhs*/, spv_double /*sol*/)
    {
      return 0;//prec->solve (matrix, rhs, sol);
    }

    /**
    * @brief setup for AMG
    *
    * @param matrix -- input matrix
    *
    * @return 0 if success
    */
    int amg_solver::setup (sp_matrix_t matrix_)
    {
      //init amg props
      BS_ASSERT (prop);
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

      t_long n = matrix->get_n_rows ();
      int n_levels;

      // matrix on first level
      a.push_back (matrix);

      for (int level = 0;;++level)
        {
          std::cout<<"AMG level = "<<level<<" n_rows = "<< n<<"\n";

          if (n <= n_last_level_points)
            {
              return lu_solver->setup (matrix);
            }

          //init tools
          sp_smbuild_t s_builder = get_smbuilder (level);
          sp_coarse_t  coarser   = get_coarser   (level);
          sp_pbuild_t  p_builder = get_pbuilder  (level);

          // build strength matrix (fill s_markers)
          spv_long sp_s_markers = BS_KERNEL.create_object (v_long::bs_type ());
          BS_ASSERT (sp_s_markers);
          s.push_back (sp_s_markers);

          s_builder->build (matrix, strength_threshold, max_row_sum, sp_s_markers);

          // coarse (fill cf_markers)
          spv_long sp_cf_markers = BS_KERNEL.create_object (v_long::bs_type ());
          BS_ASSERT (sp_cf_markers);
          cf.push_back (sp_cf_markers);
          spv_double sp_measure = BS_KERNEL.create_object (v_double::bs_type ());
          BS_ASSERT (sp_measure);

          t_long n_coarse_size = coarser->build (matrix, sp_measure, sp_cf_markers, sp_s_markers);

          if (n_coarse_size == n || n_coarse_size < 1)
            {
              n_levels = level;
              break;
            }

          // build prolongation (interpolation) matrix
          sp_bcsr_t sp_p_matrix = BS_KERNEL.create_object ("bcsr_matrix");
          BS_ASSERT (sp_p_matrix);
          p.push_back (sp_p_matrix);

          p_builder->build (matrix, n_coarse_size, max_connections,
                            sp_cf_markers, sp_s_markers, sp_p_matrix);

          // initialize next level matrix
          sp_bcsr_t a_matrix = BS_KERNEL.create_object ("bcsr_matrix");
          BS_ASSERT (a_matrix);
          a.push_back (a_matrix);

          sp_bcsr_t r_matrix = BS_KERNEL.create_object ("bcsr_matrix");
          BS_ASSERT (r_matrix);

          r_matrix->build_transpose (sp_p_matrix, 0, 0, 0);

          a_matrix->triple_matrix_product (r_matrix, matrix, sp_p_matrix, update);
          matrix = a_matrix;


        }
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
