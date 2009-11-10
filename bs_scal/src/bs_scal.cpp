// bs_scal.cpp : Defines the entry point for the DLL application.
//

#include "bs_scal_stdafx.h"
#include "scal_3p.h"
#include "py_scal_wrapper.h"

namespace blue_sky
{
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("bs_scal", "1.0.0", "BS_SCAL", "BS_SCAL", "bs_scal")

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    const plugin_descriptor & pd = *bs_init.pd_;

    bool res = true;

    res &= blue_sky::scal_register_types (pd); BS_ASSERT (res);

    return res;
  }
}//bs

#ifdef BSPY_EXPORTING_PLUGIN
BLUE_SKY_INIT_PY_FUN
{
  blue_sky::python::py_export_scal ();
}
#endif //BSPY_EXPORT_PLUGIN

