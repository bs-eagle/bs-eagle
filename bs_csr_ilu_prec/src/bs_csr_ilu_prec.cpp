// bs_csr_ilu_prec.cpp : Defines the entry point for the DLL application.
//

#include "bs_csr_ilu_prec_stdafx.h"

#include "csr_ilu_prec.h"
#include "py_csr_ilu_prec.h"

namespace blue_sky {
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("bs_csr_ilu_prec", "1.0.0", "BS_CSR_ILU_PREC", "BS_CSR_ILU_PREC", "bs_csr_ilu_prec")

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    const plugin_descriptor &pd = *bs_init.pd_;

    bool res = true;

    res &= blue_sky::csr_ilu_prec_register_type (pd); BS_ASSERT (res);

    return res;
  }
}

#ifdef BSPY_EXPORTING_PLUGIN
BLUE_SKY_INIT_PY_FUN
{
  blue_sky::python::py_export_csr_ilu_prec ();
}
#endif
