// bs_base_linear_solvers.cpp : Defines the entry point for the DLL application.
//

#include "bs_base_linear_solvers_stdafx.h"

#include "linear_solvers.h"
#include "py_linear_solvers.h"

namespace blue_sky {
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("bs_base_linear_solvers", "1.0.0", "BS_BASE_LINEAR_SOLVERS", "BS_BASE_LINEAR_SOLVERS", "bs_base_linear_solvers")

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    const plugin_descriptor &pd = *bs_init.pd_;

    bool res = true;

    res &= blue_sky::linear_solver_prop_register_type (pd); BS_ASSERT (res);
    res &= blue_sky::linear_solvers_register_type (pd); BS_ASSERT (res);

    return res;
  }
}

#ifdef BSPY_EXPORTING_PLUGIN
BLUE_SKY_INIT_PY_FUN
{
  using namespace blue_sky::python;

  py_export_linear_solver_prop ();
  py_export_linear_solvers ();
}
#endif
