// bs_base_linear_solvers.cpp : Defines the entry point for the DLL application.
//

#include "bs_base_linear_solvers_stdafx.h"

#include "linear_solvers.h"
#include "py_linear_solvers.h"
#include "throw_exception.h"

using namespace blue_sky;

namespace blue_sky {
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("bs_base_linear_solvers", "1.0.0", "BS_BASE_LINEAR_SOLVERS", "BS_BASE_LINEAR_SOLVERS", "bs_base_linear_solvers")

  namespace
  {
    bool
    register_types (plugin_descriptor const &pd)
    {
      bool res = true;

      res &= blue_sky::linear_solver_prop_register_type (pd); BS_ASSERT (res);
      res &= blue_sky::linear_solvers_register_type (pd); BS_ASSERT (res);

      return res;
    }
  }

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    return register_types (*bs_init.pd_);
  }
}

#ifdef BSPY_EXPORTING_PLUGIN
namespace 
{
  void
  init_py_subsystem ()
  {
    using namespace blue_sky::python;

    py_export_linear_solver_prop ();
    py_export_linear_solvers ();
  }
}
BLUE_SKY_INIT_PY_FUN
{
  init_py_subsystem ();
}

#ifdef _DEBUG
BOOST_PYTHON_MODULE (bs_base_linear_solvers_d)
#else
BOOST_PYTHON_MODULE (bs_base_linear_solvers)
#endif
{
  init_py_subsystem ();
  if (!blue_sky::register_types (*blue_sky::bs_get_plugin_descriptor ()))
    bs_throw_exception ("Can't register types");
}

#endif
