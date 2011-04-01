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
                : amg_solver_iface(),
                aver_cop (BS_KERNEL.create_object (v_long::bs_type ()))
    {
      //
      A.reserve (AMG_N_LEVELS_RESERVE);
      S.reserve (AMG_N_LEVELS_RESERVE);
      P.reserve (AMG_N_LEVELS_RESERVE);
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

    int amg_solver::solve_prec (sp_matrix_t matrix, spv_double rhs, spv_double sol)
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
      if (!matrix_)
        {
          bs_throw_exception ("AMG setup: Passed matrix is null");
        }

      sp_bcsr_matrix_t matrix (matrix_, bs_dynamic_cast ());

      if (!matrix)
        {
          bs_throw_exception ("AMG setup: Passed matrix is not BCSR");
        }
/*
      if (matrix->nb)
        {
          bs_throw_exception ("AMG setup: Passed matrix is not BCSR");
        }
       */
      BS_ASSERT (prop);

      // matrix on first level
      A.push_back (matrix);

      for (int level = 0;;++level)
        {
          //S.push_back (BS_KERNEL.create_object (bcsr_matrix_t::bs_type ()));

          //S[level]->init(10,10,1,10);

          // initialize next level matrix
          //A.push_back (BS_KERNEL.create_object (bcsr_matrix_t::bs_type ()));

          //A[level]->init(10,10,1,15);
          std::cout<<"AMG level = "<<level<<"\n";
          if (level > 2)
            {
              break;
            }
        }

      return 0;
    }

    //////////////////////////////////////////////////////////////////////////
    BLUE_SKY_TYPE_STD_CREATE (amg_solver);
    BLUE_SKY_TYPE_STD_COPY (amg_solver);

    BLUE_SKY_TYPE_IMPL (amg_solver, amg_solver_iface, "amg_solver", "Algebraic Multigrid linear solver and preconditioner", "Algebraic Multigrid linear solver and preconditioner");

}  // blue_sky namespace
