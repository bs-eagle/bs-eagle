/**
 *       \file  jacobian.cpp
 *      \brief  Implementation of Jacobian Matrix
 *     \author  Nikonov Max
 *       \date  27.06.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "jacobian.h"
#include BS_FORCE_PLUGIN_IMPORT()
#include "matrix_macroses.h"
#include BS_STOP_PLUGIN_IMPORT()

namespace blue_sky
  {

  /**
   * \class jacob_traits
   * \brief For sorting jacobian children
   * */
  struct jacob_traits : bs_node::sort_traits
    {
      struct jacob_key : bs_node::sort_traits::key_type
        {
          virtual bool sort_order (const key_ptr & ) const
            {
              return true;
            }
        };

      virtual const char * sort_name () const
        {
          return "jacob trait";
        };

      virtual key_ptr key_generator (const sp_link& /*l*/) const
        {
          return new jacob_key ();
        }

      virtual bool accepts (const sp_link& l)
      {
        return smart_ptr< jacobian, true >(l->data(), bs_dynamic_cast());
      }
    };

  jacobian::~jacobian ()
  {
  }

  /**
   * \brief  'default' ctor for jacobian
   * \param  param Additional params for ctor
   * */
  jacobian::jacobian (bs_type_ctor_param /*param*/)
  : bs_node(bs_node::create_node(new jacob_traits))
  {
  }

  /**
   * \brief  copy-ctor for jacobian
   * \param  src Instance of jacobian to be copied
   * */
  jacobian::jacobian (const jacobian& src)
  : bs_refcounter (src), bs_node(src)
  {
    *this = src;
  }

  void
  jacobian::setup_solver (const BS_SP (fi_params) &ts_params)
  {
    solver->get_prop()->set_i (max_iters_idx, ts_params->get_int(fi_params::LIN_ITERS_NUM));
    solver->get_prop()->set_f (tol_idx, ts_params->get_float(fi_params::LIN_SOLV_RESIDUAL));
    //solver->get_prop()->set_matbal_tolerance(ts_params->get_float(fi_params::LIN_SOLV_MATBAL_RESIDUAL));
  }

  void
  jacobian::setup_preconditioner (const BS_SP (fi_params) & /*ts_params*/)
  {
//    //set up preconditioner
//    if (prec_is_cpr_flag)
//      {
//#ifdef _MPI
//        BS_ASSERT (false && "MPI: NOT IMPL YET");
//        smart_ptr< /*!mpi_*/cpr_preconditioner<strategy_t> > cpr = ((static_cast< smart_ptr< two_stage_preconditioner<strategy_t> > >(preconditioner))->get_prec_1());
//#else // _MPI
//        smart_ptr< cpr_preconditioner<strategy_t> > cpr = ((static_cast< smart_ptr< two_stage_preconditioner<strategy_t> > >(preconditioner))->get_prec_1());
//#endif // _MPI
//
//        BS_ASSERT (cpr);
//        BS_ASSERT (cpr->get_prop ());
//
//        cpr->get_prop()->set_max_iters(ts_params->get_int(fi_params::AMG_LIN_ITERS_NUM));
//        cpr->get_prop()->set_matbal_tolerance(ts_params->get_float(fi_params::AMG_RESID));
//      }
  }

  int 
  jacobian::setup_solver_params (well_model_type /*model_type*/, int /*n_phases*/, const BS_SP (fi_params) &ts_params)
  {

    if (!ts_params)
      {
        throw bs_exception("Jacobian::setup_solver_params", "ts_params is not inited!");
      }

    if (!solver)
      {
        bs_throw_exception ("Solver is null");
      }
    if (!preconditioner)
      {
        bs_throw_exception ("Preconditioner is null");
      }

    if (solver_is_gmres_flag)
      {
        solver->get_prop ()->set_i (ortonorm_vlen, ts_params->get_int (fi_params::GMRES_ORTONORM_VLEN));
        //static_cast< smart_ptr< gmres_solver2<strategy_t> > >(solver)->m
        //= ts_params->get_int(fi_params::GMRES_ORTONORM_VLEN);
      }

    setup_solver (ts_params);
    setup_preconditioner (ts_params);

    solver->set_prec(preconditioner);
    return 0;
  }

  const BS_SP (lsolver_iface) &
  jacobian::get_solver () const
  {
    return solver;
  }

  const BS_SP (lsolver_iface) &
  jacobian::get_prec () const
  {
    return preconditioner;
  }

  void
  jacobian::init (t_long elements, t_long phases_, t_long secondary_)
  {
    secondary = secondary_;
    phases = phases_;

    matrix = BS_KERNEL.create_object ("mbcsr_matrix");
    flux_conn = BS_KERNEL.create_object ("flux_connections");

    BS_SP (bcsr_matrix_iface) flux = BS_KERNEL.create_object ("bcsr_matrix");
    BS_SP (bcsr_matrix_iface) facility = BS_KERNEL.create_object ("bcsr_matrix");
    BS_SP (bcsr_matrix_iface) accum = BS_KERNEL.create_object ("bcsr_matrix");

    matrix->add_matrix ("flux", flux);
    matrix->add_matrix ("facility", facility);
    matrix->add_matrix ("accum", accum);

    facility->alloc_rows_ptr (elements);

    accum->init (elements, elements, phases, elements);
    t_long *accum_rows = accum->get_rows_ptr ()->data ();
    t_long *accum_cols = accum->get_cols_ind ()->data ();
    t_float *accum_vals = accum->get_values ()->data ();
    // OPENMP
    for (t_long i = 0; i < elements + 1; ++i)
      {
        accum_rows[i] = i;
      }
    for (t_long i = 0; i < elements; ++i)
      {
        accum_cols[i] = i;
      }
    for (t_long i = 0; i < elements * phases * phases; ++i)
      {
        accum_vals[i] = 0;
      }

    rhs = BS_KERNEL.create_object (v_float::bs_type ());
    rhs_flux = BS_KERNEL.create_object (v_float::bs_type ());
    sec_rhs = BS_KERNEL.create_object (v_float::bs_type ());
    ss_diagonal = BS_KERNEL.create_object (v_float::bs_type ());
    sp_diagonal = BS_KERNEL.create_object (v_float::bs_type ());

    cfl_vector = BS_KERNEL.create_object (v_double::bs_type ());
    solution = BS_KERNEL.create_object (v_double::bs_type ());
    sec_solution = BS_KERNEL.create_object (v_double::bs_type ());

    boundary = BS_KERNEL.create_object (v_long::bs_type ());

    rhs->init (elements * phases, 0.0);
    rhs_flux->init (elements * phases, 0.0);
    cfl_vector->init (elements * phases, 0.0);
    solution->init (elements * phases, 0.0);

    if (secondary > 0)
      {
        sec_rhs->init (secondary * phases, 0.0);
        ss_diagonal->init (secondary * secondary * phases, 0.0);
        sp_diagonal->init (secondary * phases * phases, 0.0);
        sec_solution->init (secondary * phases, 0.0);
      }
    else
      {
        sec_rhs->init (1, 0.);
        ss_diagonal->init (1, 0.);
        sp_diagonal->init (1, 0.);
        sec_solution->init (1, 0.);
      }

    // FIXME: copied from reservoir_simulator
    //index_t N_block_size = cm->n_phases;
    //index_t N_blocks = hdm->get_mesh ()->get_n_active_elements ();
    //// FIXME: wrong init?
    //jmatrix->get_flux_matrix()->get_values ()->init (N_block_size, 0);
    //jmatrix->get_flux_matrix ()->get_values()->init (N_block_size * N_block_size * (2* hdm->get_mesh ()->get_n_connections() + hdm->get_mesh ()->get_n_active_elements()), item_t (0));

    //jmatrix->get_facility_matrix ()->init (N_blocks, N_blocks, N_block_size, 0);
  }

  BS_SP (bcsr_matrix_iface)
  jacobian::get_matrix (std::string const &name) const
  {
    BS_ASSERT (matrix);
    if (name == "trns") 
      return flux_conn->get_conn_trans ();

    return matrix->get_matrix (name);
  }

  BS_SP (mbcsr_matrix_iface)
  jacobian::get_matrix () const
  {
    BS_ASSERT (matrix);
    return matrix;
  }

  spv_float
  jacobian::get_ss_diagonal ()
  {
    return ss_diagonal;
  }

  spv_float
  jacobian::get_sp_diagonal ()
  {
    return sp_diagonal;
  }

  spv_float
  jacobian::get_sec_rhs ()
  {
    return sec_rhs;
  }

  spv_float
  jacobian::get_rhs ()
  {
    return rhs;
  }

  spv_double
  jacobian::get_rhs_dbl () const
  {
    spv_double t = BS_KERNEL.create_object(v_double::bs_type());
    t->resize(rhs->size());
    std::copy(rhs->begin(), rhs->end(), t->begin());
    return t;
  }

  spv_float
  jacobian::get_rhs_flux ()
  {
    return rhs_flux;
  }

  spv_double
  jacobian::get_cfl_vector ()
  {
    return cfl_vector;
  }

  spv_double
  jacobian::get_solution ()
  {
    return solution;
  }

  spv_double
  jacobian::get_sec_solution ()
  {
    return sec_solution;
  }

  spv_long
  jacobian::get_boundary ()
  {
    return boundary;
  }

  BS_SP (flux_connections_iface)
  jacobian::get_flux_connections ()
  {
    return flux_conn;
  }

  spv_long
  jacobian::get_m_array ()
  {
    BS_ASSERT (flux_conn);
    return flux_conn->get_matrix_block_idx_minus ();
  }

  spv_long
  jacobian::get_p_array ()
  {
    BS_ASSERT (flux_conn);
    return flux_conn->get_matrix_block_idx_plus ();
  }

  // FIXME:
  void
  jacobian::clear_solution ()
  {
    solution->assign (0);
  }

  // FIXME:
  void
  jacobian::summ_rhs ()
  {
    if (rhs->size () != rhs_flux->size ())
      {
        bs_throw_exception ("rhs and rhs_flux size mismatch");
      }

    t_float *rhs_ = rhs->data ();
    t_float *rhs_flux_ = rhs_flux->data ();

    // OPENMP:
    for (size_t i = 0, cnt = rhs->size (); i < cnt; ++i)
      {
        rhs_[i] += rhs_flux_[i];
      }
  }

  // FIXME:
  void
  jacobian::mult_flux_part (t_double mult)
  {
    spv_float facility  = matrix->get_matrix ("facility")->get_values ();
    spv_float flux      = matrix->get_matrix ("flux")->get_values ();
    t_float *facility_ = facility->data ();
    t_float *flux_     = flux->data ();
    t_float *rhs_flux_  = rhs_flux->data ();

    // FIXME: should we multiply facility?
    // OPENMP
    for (size_t i = 0, cnt = facility->size (); i < cnt; ++i)
      {
        facility_[i] *= mult;
      }
    for (size_t i = 0, cnt = flux->size (); i < cnt; ++i)
      {
        flux_[i] *= mult;
      }
    for (size_t i = 0, cnt = rhs_flux->size (); i < cnt; ++i)
      {
        rhs_flux_[i] *= mult;
      }
  }

  // FIXME:
  void
  jacobian::restore_sec_solution ()
  {
    if (secondary > 0)
      {
        if (phases != matrix->get_matrix ("flux")->get_n_block_size ())
          {
            bs_throw_exception ("phases and flux block_size mismatch");
          }

        spv_long rows     = matrix->get_matrix ("flux")->get_rows_ptr ();
        t_float *sp_diag  = sp_diagonal->data ();
        t_double *sol     = solution->data ();
        t_double *sec_sol = sec_solution->data ();
        t_float *sec_rhs_ = sec_rhs->data ();

        // OPENMP
        for (t_long i = 0, cnt = (t_long)rows->size (); i < cnt - 1; ++i)
          {
            t_float *sp_block   = &sp_diag[secondary * phases * i];
            t_double *sol_block = &sol[phases * i];
            sec_sol[i]          = sec_rhs_[i];

            VV_PROD_M (phases, sp_block, sol_block, sec_sol[i]);
          }
      }
  }

  BLUE_SKY_TYPE_STD_CREATE (jacobian);
  BLUE_SKY_TYPE_STD_COPY (jacobian);

  BLUE_SKY_TYPE_IMPL (jacobian, bs_node, "jacobian", "Jacobian", "Jacobian");
};

