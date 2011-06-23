// bs_csr_ilu_prec.cpp : Defines the entry point for the DLL application.
//

#include "bs_csr_ilu_prec_stdafx.h"

#include "csr_ilu_prec.h"
#include "py_csr_ilu_prec.h"

using namespace blue_sky;

namespace blue_sky {
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("bs_csr_ilu_prec", "1.0.0", "BS_CSR_ILU_PREC", "BS_CSR_ILU_PREC", "bs_csr_ilu_prec")

  namespace
  {
    bool
    register_types (plugin_descriptor const &pd)
    {
      bool res = true;

      res &= blue_sky::csr_ilu_prec_register_type (pd); BS_ASSERT (res);

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
    blue_sky::python::py_export_csr_ilu_prec ();
  }
}
BLUE_SKY_INIT_PY_FUN
{
  init_py_subsystem ();
}
#ifdef _DEBUG
BOOST_PYTHON_MODULE (bs_csr_ilu_prec_d)
#else
BOOST_PYTHON_MODULE (bs_csr_ilu_prec)
#endif
{
  init_py_subsystem ();
  if (!blue_sky::register_types (*blue_sky::bs_get_plugin_descriptor ()))
    bs_throw_exception ("Can't register types");
}

#endif
