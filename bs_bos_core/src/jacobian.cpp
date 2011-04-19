/**
 *       \file  jacobian.cpp
 *      \brief  Implementation of Jacobian Matrix
 *     \author  Nikonov Max
 *       \date  27.06.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"
#include "jacobian.h"

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
  jacobian::setup_preconditioner (const BS_SP (fi_params) &ts_params)
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
  jacobian::setup_solver_params (well_model_type model_type, int n_phases, const BS_SP (fi_params) &ts_params)
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
  jacobian::init (t_long elements, t_long phases, t_long secondary)
  {
    // FIXME: copied from reservoir_simulator
    //index_t N_block_size = cm->n_phases;
    //index_t N_blocks = hdm->get_mesh ()->get_n_active_elements ();
    //// FIXME: wrong init?
    //jmatrix->get_flux_matrix()->get_values ()->init (N_block_size, 0);
    //jmatrix->get_flux_matrix ()->get_values()->init (N_block_size * N_block_size * (2* hdm->get_mesh ()->get_n_connections() + hdm->get_mesh ()->get_n_active_elements()), item_t (0));

    //// FIXME: jmatrix->m_array, p_array, trns_matrix
    ////jmatrix->m_array = flux_conn->get_matrix_block_idx_minus();
    ////jmatrix->p_array = flux_conn->get_matrix_block_idx_plus ();
    ////jmatrix->trns_matrix = flux_conn->get_conn_trans();

    //jmatrix->get_facility_matrix ()->init (N_blocks, N_blocks, N_block_size, 0);
  }

  BS_SP (bcsr_matrix_iface)
  jacobian::get_matrix (std::string const &name) const
  {
    // FIXME:
    return BS_SP (bcsr_matrix_iface) ();
  }

  BS_SP (mbcsr_matrix_iface)
  jacobian::get_matrix () const
  {
    // FIXME:
    return BS_SP (mbcsr_matrix_iface) ();
  }

  BLUE_SKY_TYPE_STD_CREATE (jacobian);
  BLUE_SKY_TYPE_STD_COPY (jacobian);

  BLUE_SKY_TYPE_IMPL (jacobian, bs_node, "jacobian", "Jacobian", "Jacobian (double-long)");
};

