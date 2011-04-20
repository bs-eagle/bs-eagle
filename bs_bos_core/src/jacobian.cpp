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
// FIXME: 
//#include "two_stage_preconditioner.h"

namespace blue_sky
  {

    /**
     * \brief  'default' ctor for jacobian
     * \param  param Additional params for ctor
     * */
  jacobian::jacobian (bs_type_ctor_param /*param*/)
  : bs_node(bs_node::create_node(new jacob_traits))
  , jm (give_kernel::Instance().create_object(jac_matrix_iface::bs_type()))
  {
  }

  /**
   * \brief  copy-ctor for jacobian
   * \param  src Instance of jacobian to be copied
   * */
  jacobian::jacobian (const jacobian& src)
  : bs_refcounter (src), bs_node(src)
  , jm(src.jm)
  {
    *this = src;
  }

  void
  jacobian::setup_solver (const sp_fi_params &ts_params)
  {
    solver->get_prop()->set_i (max_iters_idx, ts_params->get_int(fi_params::LIN_ITERS_NUM));
    solver->get_prop()->set_f (tol_idx, ts_params->get_float(fi_params::LIN_SOLV_RESIDUAL));
    //solver->get_prop()->set_matbal_tolerance(ts_params->get_float(fi_params::LIN_SOLV_MATBAL_RESIDUAL));
  }

  void
  jacobian::setup_preconditioner (const sp_fi_params &ts_params)
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
  jacobian::setup_solver_params (well_model_type model_type, int n_phases, const sp_fi_params &ts_params)
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

  const jacobian::sp_lsolver &
  jacobian::get_solver () const
  {
    return solver;
  }

  const jacobian::sp_lsolver &
  jacobian::get_prec () const
  {
    return preconditioner;
  }

  //DEPRICATED
  void
  jacobian::begin ()
  {
    //this->get_jmatrix()->regular_matrix
    BS_ASSERT (false && "DEPRECATED");
  }
  //DEPRICATED
  void
  jacobian::end ()
  {
    BS_ASSERT (false && "DEPRECATED");
  }

  BLUE_SKY_TYPE_STD_CREATE (jacobian);
  BLUE_SKY_TYPE_STD_COPY (jacobian);

  BLUE_SKY_TYPE_IMPL (jacobian, bs_node, "jacobian", "Jacobian", "Jacobian (double-long)");
};

