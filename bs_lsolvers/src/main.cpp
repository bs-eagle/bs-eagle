/**
 * @file main.cpp
 * @brief
 * @author Oleg Borschuk
 * @date 2009-08-28
 */
#include "bs_lsolvers_stdafx.h"

//#include "prop.h"
#include "cgs_solver.h"
#include "gmres_solver.h"
#include "bicgstab_solver.h"
#include "tfqmr_solver.h"
#include "gs_solver.h"
#include "bcsr_ilu_prec.h"
#include "two_stage_prec.h"
#include "blu_solver.h"
//#include "py_prop.h"
#include "py_lsolver_iface.h"

using namespace blue_sky;
using namespace blue_sky::python;
using namespace boost::python;

namespace blue_sky {
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("lsolvers", "1.0.0", "Linear solvers for BS", "Linear solvers for BS", "lsolvers");

  namespace {
  bool
  register_types (const plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, cgs_solver::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, gmres_solver::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, bicgstab_solver::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, tfqmr_solver::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, gs_solver::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, bcsr_ilu_prec::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, blu_solver::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, two_stage_prec::bs_type()); BS_ASSERT (res);

    return res;
  }
  }

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    return register_types (*bs_init.pd_);
  }
}

#ifdef BSPY_EXPORTING_PLUGIN
namespace {
  void 
  init_py_subsystem ()
  {
    using namespace boost::python;

    //python::py_export_prop ();
    python::py_export_lsolvers ();
  }
}
BLUE_SKY_INIT_PY_FUN
{
  init_py_subsystem ();
}
#ifdef _DEBUG
BOOST_PYTHON_MODULE (bs_lsolvers_d)
#else
BOOST_PYTHON_MODULE (bs_lsolvers)
#endif
{
  init_py_subsystem ();
  std::cout << &BS_KERNEL << std::endl;
  bool res = register_types (*blue_sky::bs_get_plugin_descriptor ());
  if (!res)
    throw "Can't register lsolver types";
}
#endif
