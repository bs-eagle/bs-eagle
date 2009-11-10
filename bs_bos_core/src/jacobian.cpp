/*! \file jacobian.cpp
		\brief jacobian class implementations
		\author Nikonov Max
*/
#include "stdafx.h"

#include "jacobian.h"
#include BS_FORCE_PLUGIN_IMPORT ()
#include "setup_preconditioner.h"
#include BS_STOP_PLUGIN_IMPORT ()
#include "two_stage_preconditioner.h"

// WTF??
#include "well_results_storage.h"
#include "fip_results_storage.h"

namespace blue_sky
  {

  template<class strategy_t>
  jacobian<strategy_t>::jacobian (bs_type_ctor_param /*param*/)
  : bs_node(bs_node::create_node(new jacob_traits))
  , jm (give_kernel::Instance().create_object(jacobian_matrix<strategy_t>::bs_type()))
  {
  }

  template<class strategy_t>
  jacobian<strategy_t>::jacobian (const jacobian<strategy_t>& src)
  : bs_refcounter (src), bs_node(src)
  , jm(src.jm)
  {
    *this = src;
  }

  //template <typename strategy_t>
  //void
  //jacobian <strategy_t>::create_solver (int n_phases, const sp_fi_params &ts_params)
  //{
  //  int p = ts_params->get_int (fi_params::LIN_SOLVER_TYPE);

  //  if (n_phases > 1)
  //    {
  //      if (p == FI_LIN_SOLVER_BICGSTAB)
  //        {
  //          solver = BS_KERNEL.create_object (bicgstab_solver<strategy_t>::bs_type());
  //          solver_is_gmres_flag = 0;
  //        }
  //      else
  //        {
  //          #ifdef _MPI
  //          BS_ASSERT (false && "MPI: NOT IMPL YET");
  //          //!solver = give_kernel::Instance().create_object(mpi_gmres_solver<strategy_t>::bs_type());
  //          #else //_MPI
  //          solver = BS_KERNEL.create_object (gmres_solver2<strategy_t>::bs_type());
  //          #endif //_MPI
  //          solver_is_gmres_flag = 1;
  //        }
  //    }
  //  else
  //    {
  //      solver = BS_KERNEL.create_object (amg_solver <strategy_t>::bs_type ());
  //      smart_ptr <amg_solver <strategy_t> > amg (solver, bs_static_cast ());
  //      
  //      amg->amg_prop = BS_KERNEL.create_object (amg_properties::bs_type ());
  //      amg->set_coarse(BS_KERNEL.create_object (coarse_pmis_2 <strategy_t> ::bs_type ()),true);
  //      amg->set_interpolator(BS_KERNEL.create_object (interpolator_standart_2 <strategy_t> ::bs_type ()),true);

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

//  template <typename strategy_t>
//  void
//  jacobian <strategy_t>::create_preconditioner (well_model_type model_type, int n_phases, const sp_fi_params &ts_params)
//  {
//    if (n_phases > 1)
//      {
//        int p = ts_params->get_int (fi_params::PREC_TYPE);
//        if (p == FI_LIN_PREC_ILU)
//          {
//#ifdef _MPI
//            BS_ASSERT (false && "MPI: NOT IMPL YET");
//            //!preconditioner = give_kernel::Instance().create_object(mpi_csr_ilu_prec<strategy_t>::bs_type());
//#else //_MPI
//#ifdef ILU_PREC_PARALLEL
//            preconditioner = give_kernel::Instance().create_object(csr_pilu_prec<strategy_t>::bs_type());
//#else //ILU_PREC_PARALLEL
//            preconditioner = give_kernel::Instance().create_object(csr_ilu_prec<strategy_t>::bs_type());
//            //preconditioner = give_kernel::Instance().create_object(csr_ilu_cfl_prec<strategy_t>::bs_type());
//#endif //ILU_PREC_PARALLEL
//#endif //_MPI
//            prec_is_cpr_flag = 0;
//          } // if (p == FI_LIN_PREC_ILU)
//        else if (p == FI_LIN_PREC_CPR_SOR)
//          {
//            preconditioner = give_kernel::Instance().create_object(two_stage_preconditioner<strategy_t>::bs_type());
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
//                //!static_cast< smart_ptr< two_stage_preconditioner<strategy_t> > >(preconditioner)->
//                //!set_prec_1 (give_kernel::Instance().create_object(mpi_cpr_preconditioner<strategy_t>::bs_type()));
//              }
//            else
//              {
//                throw bs_exception ("jacobian::create_preconditioner", "model type not supported");
//              }
//            //!static_cast< smart_ptr< two_stage_preconditioner<strategy_t> > >(preconditioner)->
//            //!set_prec_2 (give_kernel::Instance().create_object(mpi_csr_ilu_prec<strategy_t>::bs_type()));
//#else // _MPI
//            if ((model_type == BLACK_OIL) || (model_type == COMPOSIT))
//              {
//                sp_obj cpr_prec_raw (give_kernel::Instance().create_object(cpr_preconditioner<strategy_t>::bs_type()));
//                smart_ptr <cpr_preconditioner <strategy_t> > cpr_prec (cpr_prec_raw, bs_dynamic_cast ());
//                smart_ptr <two_stage_preconditioner <strategy_t> > prec (preconditioner, bs_dynamic_cast ());
//                prec->set_prec_1 (cpr_prec.get ());
//              }
//            else
//              {
//                bs_throw_exception ("Model type not supported");
//              }
//            //////////!!!!!!!!!!!!!!!!!!!!!! sor_prec need for mpi_vector
//            //static_cast< smart_ptr< two_stage_preconditioner<strategy_t> > >(preconditioner)->
//            //set_prec_2 (give_kernel::Instance().create_object(sor_prec<strategy_t>::bs_type()));
//#endif // _MPI
//          } // if (p == FI_LIN_PREC_CPR_SOR)
//        else
//          {
//            BOSOUT (section::solvers, level::debug) << "create preconditioner" << bs_end;
//            preconditioner = give_kernel::Instance().create_object(two_stage_preconditioner<strategy_t>::bs_type());
//            if (!preconditioner)
//              {
//                bs_throw_exception ("Can't create preconditioner");
//              }
//
//            prec_is_cpr_flag = 1;
//            smart_ptr <two_stage_preconditioner <strategy_t> > prec (preconditioner, bs_dynamic_cast ());
//
//#ifdef _MPI
//            BS_ASSERT (false && "MPI: NOT IMPL YET");
//            if (model_type == BLACK_OIL)
//              {
//                //!static_cast< smart_ptr< two_stage_preconditioner<strategy_t> > >(preconditioner)->
//                //!set_prec_1 (give_kernel::Instance().create_object(mpi_cpr_preconditioner<strategy_t>::bs_type()));
//              }
//            else
//              {
//                throw bs_exception ("jacobian::create_preconditioner", "model type not supported");
//              }
//            //!static_cast< smart_ptr< two_stage_preconditioner<strategy_t> > >(preconditioner)->
//            //!set_prec_2 (give_kernel::Instance().create_object(mpi_csr_ilu_prec<strategy_t>::bs_type()));
//#else // _MPI
//            sp_obj cpr_prec_raw (give_kernel::Instance().create_object(cpr_preconditioner<strategy_t>::bs_type()));
//            smart_ptr <cpr_preconditioner <strategy_t> > cpr_prec (cpr_prec_raw, bs_static_cast ());
//            BS_ASSERT (cpr_prec_raw);
//            BS_ASSERT (cpr_prec);
//
//            prec->set_prec_1 (cpr_prec.get ());
//#ifdef ILU_PREC_PARALLEL
//            static_cast< smart_ptr< two_stage_preconditioner<strategy_t> > >(preconditioner)->
//            set_prec_2 (give_kernel::Instance().create_object(csr_pilu_prec<strategy_t>::bs_type()));
//#else // ILU_PREC_PARALLEL
//#ifdef BS_BOS_CORE_USE_CSR_ILU_CFL_PREC
//            sp_obj csr_prec_raw (give_kernel::Instance().create_object(csr_ilu_cfl_prec<strategy_t>::bs_type()));
//            smart_ptr <csr_ilu_cfl_prec <strategy_t> > csr_prec (csr_prec_raw, bs_static_cast ());
//#else
//            sp_obj csr_prec_raw (BS_KERNEL.create_object (csr_ilu_prec <strategy_t>::bs_type ()));
//            smart_ptr <csr_ilu_prec <strategy_t> > csr_prec (csr_prec_raw, bs_static_cast ());
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
//            //!preconditioner = give_kernel::Instance().create_object(mpi_csr_ilu_prec<strategy_t>::bs_type());
//#else // _MPI
//#ifdef ILU_PREC_PARALLEL
//            preconditioner = give_kernel::Instance().create_object(csr_pilu_prec<strategy_t>::bs_type());
//#else // ILU_PREC_PARALLEL
//            preconditioner = give_kernel::Instance().create_object(csr_ilu_prec<strategy_t>::bs_type());
//#endif // ILU_PREC_PARALLEL
//#endif // _MPI
//            prec_is_cpr_flag = 0;
//          }
//        else if (p == FI_LIN_PREC_AMG)
//          {
//#ifdef _MPI
//            BS_ASSERT (false && "MPI: NOT IMPL YET");
//            //!preconditioner = give_kernel::Instance().create_object(mpi_amg_solver<strategy_t>::bs_type());
//#else // _MPI
//            preconditioner = give_kernel::Instance().create_object(amg_solver<strategy_t>::bs_type());
//
//#endif // _MPI
//          }
//        else
//          {
//            throw bs_exception ("jacobian::create_preconditioner", "this preconditioner not supported for 1 phase model");
//          }
//      }
//  }

  template <typename strategy_t>
  void
  jacobian <strategy_t>::setup_solver (const sp_fi_params &ts_params)
  {
    solver->get_prop()->set_max_iters(ts_params->get_int(fi_params::LIN_ITERS_NUM));
    solver->get_prop()->set_tolerance(ts_params->get_float(fi_params::LIN_SOLV_RESIDUAL));
    //solver->get_prop()->set_matbal_tolerance(ts_params->get_float(fi_params::LIN_SOLV_MATBAL_RESIDUAL));
  }

  template <typename strategy_t>
  void
  jacobian <strategy_t>::setup_preconditioner (const sp_fi_params &ts_params)
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

  template<class strategy_t>
  int jacobian<strategy_t>::setup_solver_params (well_model_type model_type, int n_phases, const sp_fi_params &ts_params)
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
        static_cast< smart_ptr< gmres_solver2<strategy_t> > >(solver)->m
        = ts_params->get_int(fi_params::GMRES_ORTONORM_VLEN);
      }

    setup_solver (ts_params);
    setup_preconditioner (ts_params);

    solver->set_prec(preconditioner);
    return 0;
  }

  template <typename strategy_t>
  const typename jacobian <strategy_t>::sp_lsolver &
  jacobian <strategy_t>::get_solver () const
  {
    return solver;
  }
  template <typename strategy_t>
  const typename jacobian <strategy_t>::sp_lsolver &
  jacobian <strategy_t>::get_prec () const
  {
    return preconditioner;
  }

  //DEPRICATED
  template <typename strategy_t>
  void
  jacobian<strategy_t>::begin ()
  {
    //this->get_jmatrix()->regular_matrix
    BS_ASSERT (false && "NOT IMPL YET");
  }
  //DEPRICATED
  template <typename strategy_t>
  void
  jacobian<strategy_t>::end ()
  {
    BS_ASSERT (false && "NOT IMPL YET");
  }

  BLUE_SKY_TYPE_STD_CREATE_T_DEF(jacobian, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF(jacobian, (class));

  BLUE_SKY_TYPE_IMPL_T_EXT(1, (jacobian<base_strategy_fi>) , 1, (objbase), "jacobian_fi", "Jacobian-float-int", "Jacobian with float items and integer indexes", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (jacobian<base_strategy_di>) , 1, (objbase), "jacobian_di", "Jacobian-double-int", "Jacobian with double items and integer indexes", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (jacobian<base_strategy_mixi>) , 1, (objbase), "jacobian_mixi", "Jacobian-mixi", "Jacobian mixi", false);

};

