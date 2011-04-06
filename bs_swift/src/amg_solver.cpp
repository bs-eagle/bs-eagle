/**
 * @file amg_solver.cpp
 * @brief implementation of AMG solver
 * @author
 * @version
 * @date 2011-03-30
 */

#include "amg_solver.h"
#include "simple_smbuilder.h"

#define AMG_N_LEVELS_RESERVE 10

namespace blue_sky
{
    //! constructor
    amg_solver::amg_solver (bs_type_ctor_param)
                : amg_solver_iface ()
    {
      strength_type = BS_KERNEL.create_object (v_int::bs_type ());
      coarse_type   = BS_KERNEL.create_object (v_int::bs_type ());
      interp_type   = BS_KERNEL.create_object (v_int::bs_type ());
      //
      A.reserve (AMG_N_LEVELS_RESERVE);
      P.reserve (AMG_N_LEVELS_RESERVE);
      s.reserve (AMG_N_LEVELS_RESERVE);
      cf.reserve (AMG_N_LEVELS_RESERVE);
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
      BS_ASSERT (prop);

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
      t_long nb = matrix->get_n_block_size ();
      int n_levels;

      // matrix on first level
      A.push_back (matrix);

      for (int level = 0;;++level)
        {
          spv_long sp_s_markers = BS_KERNEL.create_object (v_long::bs_type ());
          BS_ASSERT (sp_s_markers);

          sp_smbuild_t s_builder;
          if (get_strength_type (level) == 0)
            {
               //simple_smbuilder s_builder;
               s_builder = BS_KERNEL.create_object (simple_smbuilder::bs_type ());
            }
          t_double strength_threshold = 0.25;
          t_double max_row_sum = 0.01;

          s_builder->build (matrix, strength_threshold, max_row_sum, sp_s_markers);
          s.push_back (sp_s_markers);

          spv_long cf_markers = BS_KERNEL.create_object (v_long::bs_type ());
          cf.push_back (cf_markers);

          // initialize next level matrix
          sp_bcsr_t a_matrix = BS_KERNEL.create_object ("bcsr_matrix");
          BS_ASSERT (a_matrix);
          A.push_back (a_matrix);

          std::cout<<"AMG level = "<<level<<"\n";
          if (level > 2)
            {
              n_levels = level;
              break;
            }
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
