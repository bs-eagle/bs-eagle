/** 
 * @file main.cpp
 * @brief 
 * @author Oleg Borschuk
 * @date 2009-08-28
 */
#include "bs_lsolvers_stdafx.h"
//#include "bs_kernel.h"

#include "prop.h"
#include "cgs_solver.h"
#include "gmres_solver.h"
#include "bicgstab_solver.h"
#include "tfqmr_solver.h"
#include "bcsr_ilu_prec.h"
#include "two_stage_prec.h"
#include "blu_solver.h"
#include "py_prop.h"
#include "py_iface.h"

using namespace blue_sky;
using namespace blue_sky::python;
using namespace boost::python;

namespace blue_sky {
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("lsolvers", "1.0.0", "Linear solvers for BS", "Linear solvers for BS", "lsolvers")

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    //const plugin_descriptor &pd = *bs_init.pd_;

    bool res = true;

    //res &= BS_KERNEL.register_type(*bs_init.pd_, matrix_base<seq_vector<float>, seq_vector<int> >::bs_type()); BS_ASSERT (res);
    //res &= BS_KERNEL.register_type(*bs_init.pd_, matrix_base<seq_vector<double>, seq_vector<int> >::bs_type()); BS_ASSERT (res);

    //res &= BS_KERNEL.register_type(*bs_init.pd_, bcsr_matrix<seq_vector<float>, seq_vector<int> >::bs_type()); BS_ASSERT (res);
    //res &= BS_KERNEL.register_type(*bs_init.pd_, bcsr_matrix<seq_vector<double>, seq_vector<int> >::bs_type()); BS_ASSERT (res);

    res &= BS_KERNEL.register_type(*bs_init.pd_, prop<float, int, std::string, bool>::bs_type()); BS_ASSERT (res);

    res &= BS_KERNEL.register_type(*bs_init.pd_, cgs_solver<base_strategy_fif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, cgs_solver<base_strategy_did>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, cgs_solver<base_strategy_dif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, cgs_solver<base_strategy_flf>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, cgs_solver<base_strategy_dld>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, cgs_solver<base_strategy_dlf>::bs_type()); BS_ASSERT (res);

    res &= BS_KERNEL.register_type(*bs_init.pd_, gmres_solver<base_strategy_fif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, gmres_solver<base_strategy_did>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, gmres_solver<base_strategy_dif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, gmres_solver<base_strategy_flf>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, gmres_solver<base_strategy_dld>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, gmres_solver<base_strategy_dlf>::bs_type()); BS_ASSERT (res);

    res &= BS_KERNEL.register_type(*bs_init.pd_, bicgstab_solver<base_strategy_fif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bicgstab_solver<base_strategy_did>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bicgstab_solver<base_strategy_dif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bicgstab_solver<base_strategy_flf>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bicgstab_solver<base_strategy_dld>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bicgstab_solver<base_strategy_dlf>::bs_type()); BS_ASSERT (res);

    res &= BS_KERNEL.register_type(*bs_init.pd_, tfqmr_solver<base_strategy_fif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, tfqmr_solver<base_strategy_did>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, tfqmr_solver<base_strategy_dif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, tfqmr_solver<base_strategy_flf>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, tfqmr_solver<base_strategy_dld>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, tfqmr_solver<base_strategy_dlf>::bs_type()); BS_ASSERT (res);

    res &= BS_KERNEL.register_type(*bs_init.pd_, bcsr_ilu_prec<base_strategy_fif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bcsr_ilu_prec<base_strategy_did>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bcsr_ilu_prec<base_strategy_dif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bcsr_ilu_prec<base_strategy_flf>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bcsr_ilu_prec<base_strategy_dld>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bcsr_ilu_prec<base_strategy_dlf>::bs_type()); BS_ASSERT (res);

    res &= BS_KERNEL.register_type(*bs_init.pd_, blu_solver<base_strategy_fif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, blu_solver<base_strategy_did>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, blu_solver<base_strategy_dif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, blu_solver<base_strategy_flf>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, blu_solver<base_strategy_dld>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, blu_solver<base_strategy_dlf>::bs_type()); BS_ASSERT (res);

    res &= BS_KERNEL.register_type(*bs_init.pd_, two_stage_prec<base_strategy_fif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, two_stage_prec<base_strategy_did>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, two_stage_prec<base_strategy_dif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, two_stage_prec<base_strategy_flf>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, two_stage_prec<base_strategy_dld>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, two_stage_prec<base_strategy_dlf>::bs_type()); BS_ASSERT (res);

    return res;
  }
}

#ifdef BSPY_EXPORTING_PLUGIN
BLUE_SKY_INIT_PY_FUN
{
  using namespace boost::python;

  python::py_export_prop ();
  python::py_export_lsolvers ();
  
}
#endif
