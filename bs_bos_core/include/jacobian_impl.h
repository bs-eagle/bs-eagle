/**
 * \file jacobian_impl.h
 * \brief impl of jacobian method
 * \author Sergey Miryanov
 * \date 23.01.2009
 * */
#ifndef BS_JACOBIAN_IMPL_H_
#define BS_JACOBIAN_IMPL_H_
#include "solve_helper.h"

namespace blue_sky {

  template <typename strategy_t>
  struct jacobian_impl
  {
    typedef typename strategy_t::item_t         item_t;
    typedef typename strategy_t::index_t        index_t;
    typedef typename strategy_t::matrix_t       matrix_t;
    typedef typename strategy_t::csr_matrix_t   bcsr_matrix_t;
    typedef typename strategy_t::item_array_t   item_array_t;
    typedef typename strategy_t::rhs_item_array_t rhs_item_array_t;

    typedef jacobian <strategy_t>               jacobian_t;
    typedef jacobian_matrix <strategy_t>        jmatrix_t;
    typedef linear_solver_base <strategy_t>     linear_solver_t;

    typedef smart_ptr <jacobian_t, true>        sp_jacobian_t;
    typedef smart_ptr <jmatrix_t, true>         sp_jmatrix_t;
    typedef smart_ptr <linear_solver_t, true>   sp_linear_solver_t;
    typedef smart_ptr <matrix_t>                sp_matrix_t;
    typedef smart_ptr <bcsr_matrix_t, true>     sp_bcsr_matrix_t;

    typedef sp_linear_solver_t                  sp_solver_t;

    jacobian_impl (sp_jacobian_t &jacobian, sp_jmatrix_t &jmatrix)
    : jacobian_ (jacobian),
    jmatrix_ (jmatrix),
    solver_ (jacobian_->get_solver ()),
    preconditioner_ (jacobian_->get_prec ()),
    rhs_ (jmatrix_->get_rhs ()),
    sol_ (jmatrix_->get_solution ()),
    regular_matrix_ (jmatrix_->get_regular_matrix ()),
    irregular_matrix_ (jmatrix_->get_irregular_matrix ())
    {
      if (!jmatrix_)
        {
          bs_throw_exception ("Jacobian matrix is not inited!");
        }
      if (!solver_)
        {
          bs_throw_exception ("Solver is not inited!");
        }
      if (!preconditioner_)
        {
          bs_throw_exception ("Preconditioner is not inited!");
        }
    }

    void
    setup_jacobian ()
    {
      full_matrix_print ();

      jmatrix_->prepare_matrix ();

      OMP_TIME_MEASURE_START(gmres_setup_timer);
      if (0 != solver_->setup (jmatrix_))
        {
          bs_throw_exception ("Internal error: cannot build preconditioner for linear solver.");
        }
      OMP_TIME_MEASURE_END (gmres_setup_timer);
    }

    item_t
    solve_jacobian (index_t &n_lin_iters)
    {
#ifdef _MPI
      BS_ASSERT (false && "MPI: NOT IMPL YET");
      //!ret_code = solver->solve (jm, jm->get_mpi_rhs(), jm->get_mpi_sol());
#else //_MPI
      sp_matrix_t working_matrix;
      //const blue_sky::setup_preconditioner<matrix_t> *setup_prec = dynamic_cast<const blue_sky::setup_preconditioner<matrix_t> *> ((objbase *)jmatrix_.get ());
      //if (setup_prec)
      //  {
      //    working_matrix = setup_prec->get_prepared_matrix ();
      //  }
      //else
      //  {
      //    working_matrix = jmatrix_.get ();
      //  }

      // TODO: BUG
      working_matrix = jmatrix_->get_prepared_matrix ();
      const sp_matrix_t &locked_working_mx (working_matrix);

      debug::print_memory_info ("-> jacobian_solver_solve");
      index_t ret_code = solve_helper (&*solver_, &(*locked_working_mx), rhs_, sol_);
      debug::print_memory_info ("<- jacobian_solver_solve");
#endif //_MPI

      jmatrix_->restore_sec_solution ();

      static size_t iter_counter = 0;
      ++iter_counter;
      ////jmatrix_->get_regular_matrix ()->ascii_write_in_csr_format (tools::string_formater ("j_reg_mx.bs.%d.txt", iter_counter).str);
      ////jmatrix_->get_irregular_matrix ()->ascii_write_in_csr_format (tools::string_formater ("j_irreg_mx.bs.%d.txt", iter_counter).str);
      //tools::save_seq_vector (tools::string_formater ("j_rhs.bs.%d.txt", iter_counter).str).save (rhs_);
      //tools::save_seq_vector (tools::string_formater ("j_rhs_flux.bs.%d.txt", iter_counter).str).save (jmatrix_->get_rhs_flux ());
      //tools::save_seq_vector (tools::string_formater ("j_sol.bs.%d.txt", iter_counter).str).save (sol_);
      //tools::save_seq_vector (tools::string_formater ("j_sec_sol.bs.%d.txt", iter_counter).str).save (jmatrix_->get_sec_solution ());

      if (ret_code < 0)
        {
          BOSERR (section::solvers, level::error) << "Linear solver failed with retcode = " << ret_code << bs_end;
          return 10.0f;
        }

#ifdef __FULL_MATRIX_PRINT__
      ++cur_iter;
#endif //__FULL_MATRIX_PRINT__

      n_lin_iters = solver_->get_prop()->get_iters();
      item_t tol = solver_->get_prop()->get_final_resid();
      if (tol > solver_->get_prop()->get_tolerance())
        {
          BOSERR (section::solvers, level::error)
            << "Linear solver failed with tolerance "
            << tol
            << bs_end;
        }
      else
        {
          BOSOUT (section::solvers, level::medium)
            << "Linear solver iterations " 
            << n_lin_iters 
            << ", tol = " << tol 
            << ", ret_code = " << ret_code 
            << bs_end;
        }

      return tol;
    }

    void
    full_matrix_print ()
    {
#ifdef __FULL_MATRIX_PRINT__
      if (jacobian_->cur_iter % jacobian_->sav_iters == 0)
        {
          char reg_filename [50];
          char irr_filename [50];

          sprintf (reg_filename, "regular_matrix%d.bs.b", jacobian_->cur_iter);
          sprintf (irr_filename, "irregular_matrix%d.bs.csr", jacobian_->cur_iter);

          regular_matrix_->ascii_write_in_csr_format(reg_filename);
          irregular_matrix_->ascii_write_in_csr_format(irr_filename);
        }
#endif
    }

  private:

    sp_jacobian_t     &jacobian_;
    sp_jmatrix_t      &jmatrix_;
    sp_solver_t       solver_;
    sp_solver_t       preconditioner_;

    rhs_item_array_t  &rhs_;
    item_array_t      &sol_;

    sp_matrix_t       regular_matrix_;
    sp_matrix_t       irregular_matrix_;
  };


} // namespace blue_sky


#endif  // #ifndef BS_JACOBIAN_IMPL_H_

