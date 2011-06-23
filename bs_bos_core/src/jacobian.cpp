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
#include BS_FORCE_PLUGIN_IMPORT ()
#include "setup_preconditioner.h"
#include BS_STOP_PLUGIN_IMPORT ()
#include "two_stage_preconditioner.h"

namespace blue_sky
  {

    /**
     * \brief  'default' ctor for jacobian
     * \param  param Additional params for ctor
     * */
  jacobian::jacobian (bs_type_ctor_param /*param*/)
  : bs_node(bs_node::create_node(new jacob_traits))
  , jm (give_kernel::Instance().create_object(jacobian_matrix::bs_type()))
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

  //void
  //jacobian::create_solver (int n_phases, const sp_fi_params &ts_params)
  //{
  //  int p = ts_params->get_int (fi_params::LIN_SOLVER_TYPE);

  //  if (n_phases > 1)
  //    {
  //      if (p == FI_LIN_SOLVER_BICGSTAB)
  //        {
  //          solver = BS_KERNEL.create_object (bicgstab_solver::bs_type());
  //          solver_is_gmres_flag = 0;
  //        }
  //      else
  //        {
  //          #ifdef _MPI
  //          BS_ASSERT (false && "MPI: NOT IMPL YET");
  //          //!solver = give_kernel::Instance().create_object(mpi_gmres_solver::bs_type());
  //          #else //_MPI
  //          solver = BS_KERNEL.create_object (gmres_solver2::bs_type());
  //          #endif //_MPI
  //          solver_is_gmres_flag = 1;
  //        }
  //    }
  //  else
  //    {
  //      solver = BS_KERNEL.create_object (amg_solver::bs_type ());
  //      smart_ptr <amg_solver> amg (solver, bs_static_cast ());
  //      
  //      amg->amg_prop = BS_KERNEL.create_object (amg_properties::bs_type ());
  //      amg->set_coarse(BS_KERNEL.create_object (coarse_pmis_2::bs_type ()),true);
  //      amg->set_interpolator(BS_KERNEL.create_object (interpolator_standart_2::bs_type ()),true);

  //      BS_ASSERT (amg->amg_prop);
  //      amg->amg_prop->set_int (amg_properties::I_ADAPTIVE_TRESHOLD, 0);
  //      amg->amg_prop->set_int (amg_properties::I_UPDATE, 0);
  //      amg->amg_prop->set_float (amg_properties::FP_STRENGHT_THRESHOLD, 0.75);
  //    }

  //  if (!solver)
  //    {
  //      bs_throw_exception ("Can't create solver");
  //    }
  //}

//  void
//  jacobian::create_preconditioner (well_model_type model_type, int n_phases, const sp_fi_params &ts_params)
//  {
//    if (n_phases > 1)
//      {
//        int p = ts_params->get_int (fi_params::PREC_TYPE);
//        if (p == FI_LIN_PREC_ILU)
//          {
//#ifdef _MPI
//            BS_ASSERT (false && "MPI: NOT IMPL YET");
//            //!preconditioner = give_kernel::Instance().create_object(mpi_csr_ilu_prec::bs_type());
//#else //_MPI
//#ifdef ILU_PREC_PARALLEL
//            preconditioner = give_kernel::Instance().create_object(csr_pilu_prec::bs_type());
//#else //ILU_PREC_PARALLEL
//            preconditioner = give_kernel::Instance().create_object(csr_ilu_prec::bs_type());
//            //preconditioner = give_kernel::Instance().create_object(csr_ilu_cfl_prec::bs_type());
//#endif //ILU_PREC_PARALLEL
//#endif //_MPI
//            prec_is_cpr_flag = 0;
//          } // if (p == FI_LIN_PREC_ILU)
//        else if (p == FI_LIN_PREC_CPR_SOR)
//          {
//            preconditioner = give_kernel::Instance().create_object(two_stage_preconditioner::bs_type());
//            if (!preconditioner)
//              {
//                bs_throw_exception ("Can't create preconditioner");
//              }
//
//            prec_is_cpr_flag = 1;
//#ifdef _MPI
//            BS_ASSERT (false && "MPI: NOT IMPL YET");
//            if (model_type == BLACK_OIL)
//              {
//                //!static_cast< smart_ptr< two_stage_preconditioner> >(preconditioner)->
//                //!set_prec_1 (give_kernel::Instance().create_object(mpi_cpr_preconditioner::bs_type()));
//              }
//            else
//              {
//                throw bs_exception ("jacobian::create_preconditioner", "model type not supported");
//              }
//            //!static_cast< smart_ptr< two_stage_preconditioner> >(preconditioner)->
//            //!set_prec_2 (give_kernel::Instance().create_object(mpi_csr_ilu_prec::bs_type()));
//#else // _MPI
//            if ((model_type == BLACK_OIL) || (model_type == COMPOSIT))
//              {
//                sp_obj cpr_prec_raw (give_kernel::Instance().create_object(cpr_preconditioner::bs_type()));
//                smart_ptr <cpr_preconditioner> cpr_prec (cpr_prec_raw, bs_dynamic_cast ());
//                smart_ptr <two_stage_preconditioner> prec (preconditioner, bs_dynamic_cast ());
//                prec->set_prec_1 (cpr_prec.get ());
//              }
//            else
//              {
//                bs_throw_exception ("Model type not supported");
//              }
//            //////////!!!!!!!!!!!!!!!!!!!!!! sor_prec need for mpi_vector
//            //static_cast< smart_ptr< two_stage_preconditioner> >(preconditioner)->
//            //set_prec_2 (give_kernel::Instance().create_object(sor_prec::bs_type()));
//#endif // _MPI
//          } // if (p == FI_LIN_PREC_CPR_SOR)
//        else
//          {
//            BOSOUT (section::solvers, level::debug) << "create preconditioner" << bs_end;
//            preconditioner = give_kernel::Instance().create_object(two_stage_preconditioner::bs_type());
//            if (!preconditioner)
//              {
//                bs_throw_exception ("Can't create preconditioner");
//              }
//
//            prec_is_cpr_flag = 1;
//            smart_ptr <two_stage_preconditioner> prec (preconditioner, bs_dynamic_cast ());
//
//#ifdef _MPI
//            BS_ASSERT (false && "MPI: NOT IMPL YET");
//            if (model_type == BLACK_OIL)
//              {
//                //!static_cast< smart_ptr< two_stage_preconditioner> >(preconditioner)->
//                //!set_prec_1 (give_kernel::Instance().create_object(mpi_cpr_preconditioner::bs_type()));
//              }
//            else
//              {
//                throw bs_exception ("jacobian::create_preconditioner", "model type not supported");
//              }
//            //!static_cast< smart_ptr< two_stage_preconditioner> >(preconditioner)->
//            //!set_prec_2 (give_kernel::Instance().create_object(mpi_csr_ilu_prec::bs_type()));
//#else // _MPI
//            sp_obj cpr_prec_raw (give_kernel::Instance().create_object(cpr_preconditioner::bs_type()));
//            smart_ptr <cpr_preconditioner> cpr_prec (cpr_prec_raw, bs_static_cast ());
//            BS_ASSERT (cpr_prec_raw);
//            BS_ASSERT (cpr_prec);
//
//            prec->set_prec_1 (cpr_prec.get ());
//#ifdef ILU_PREC_PARALLEL
//            static_cast< smart_ptr< two_stage_preconditioner> >(preconditioner)->
//            set_prec_2 (give_kernel::Instance().create_object(csr_pilu_prec::bs_type()));
//#else // ILU_PREC_PARALLEL
//#ifdef BS_BOS_CORE_USE_CSR_ILU_CFL_PREC
//            sp_obj csr_prec_raw (give_kernel::Instance().create_object(csr_ilu_cfl_prec::bs_type()));
//            smart_ptr <csr_ilu_cfl_prec> csr_prec (csr_prec_raw, bs_static_cast ());
//#else
//            sp_obj csr_prec_raw (BS_KERNEL.create_object (csr_ilu_prec::bs_type ()));
//            smart_ptr <csr_ilu_prec> csr_prec (csr_prec_raw, bs_static_cast ());
//#endif
//            BS_ASSERT (csr_prec_raw);
//            BS_ASSERT (csr_prec);
//
//            prec->set_prec_2 (csr_prec.get ());
//#endif // ILU_PREC_PARALLEL
//#endif // _MPI
//          }
//      }
//    else   // if (n_phases > 1) - pressure system
//      {
//        int p = ts_params->get_int(fi_params::PREC_TYPE_ONE_PHASE);
//        if (p == FI_LIN_PREC_ILU)
//          {
//#ifdef _MPI
//            BS_ASSERT (false && "MPI: NOT IMPL YET");
//            //!preconditioner = give_kernel::Instance().create_object(mpi_csr_ilu_prec::bs_type());
//#else // _MPI
//#ifdef ILU_PREC_PARALLEL
//            preconditioner = give_kernel::Instance().create_object(csr_pilu_prec::bs_type());
//#else // ILU_PREC_PARALLEL
//            preconditioner = give_kernel::Instance().create_object(csr_ilu_prec::bs_type());
//#endif // ILU_PREC_PARALLEL
//#endif // _MPI
//            prec_is_cpr_flag = 0;
//          }
//        else if (p == FI_LIN_PREC_AMG)
//          {
//#ifdef _MPI
//            BS_ASSERT (false && "MPI: NOT IMPL YET");
//            //!preconditioner = give_kernel::Instance().create_object(mpi_amg_solver::bs_type());
//#else // _MPI
//            preconditioner = give_kernel::Instance().create_object(amg_solver::bs_type());
//
//#endif // _MPI
//          }
//        else
//          {
//            throw bs_exception ("jacobian::create_preconditioner", "this preconditioner not supported for 1 phase model");
//          }
//      }
//  }

  void
  jacobian::setup_solver (const sp_fi_params &ts_params)
  {
    solver->get_prop()->set_max_iters(ts_params->get_int(fi_params::LIN_ITERS_NUM));
    solver->get_prop()->set_tolerance(ts_params->get_float(fi_params::LIN_SOLV_RESIDUAL));
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
//        smart_ptr< /*!mpi_*/cpr_preconditioner> cpr = ((static_cast< smart_ptr< two_stage_preconditioner> >(preconditioner))->get_prec_1());
//#else // _MPI
//        smart_ptr< cpr_preconditioner> cpr = ((static_cast< smart_ptr< two_stage_preconditioner> >(preconditioner))->get_prec_1());
//#endif // _MPI
//
//        BS_ASSERT (cpr);
//        BS_ASSERT (cpr->get_prop ());
//
//        cpr->get_prop()->set_max_iters(ts_params->get_int(fi_params::AMG_LIN_ITERS_NUM));
//        cpr->get_prop()->set_matbal_tolerance(ts_params->get_float(fi_params::AMG_RESID));
//      }
  }

  int jacobian::setup_solver_params (well_model_type model_type, int n_phases, const sp_fi_params &ts_params)
  {

    if (!ts_params)
      {
        throw bs_exception("Jacobian::setup_solver_params", "ts_params is not inited!");
      }

    //if (ts_params->check_value (fi_params::LIN_SOLVER_TYPE) && !solver)
    //  {
    //    create_solver (n_phases, ts_params);
    //  }

    //if (ts_params->check_value (fi_params::PREC_TYPE) && !preconditioner)
    //  {
    //    create_preconditioner (model_type, n_phases, ts_params);
    //  }

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
        static_cast< smart_ptr< gmres_solver2> >(solver)->m
        = ts_params->get_int(fi_params::GMRES_ORTONORM_VLEN);
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
    BS_ASSERT (false && "NOT IMPL YET");
  }
  //DEPRICATED
  void
  jacobian::end ()
  {
    BS_ASSERT (false && "NOT IMPL YET");
  }

  BLUE_SKY_TYPE_STD_CREATE (jacobian);
  BLUE_SKY_TYPE_STD_COPY (jacobian);

  BLUE_SKY_TYPE_IMPL (jacobian, objbase, "jacobian", "jacobian", "jacobian");

};

