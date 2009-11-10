// bs_pvt.cpp : Defines the entry point for the DLL application.
//

#include "bs_pvt_stdafx.h"

#include "pvt_base.h"
#include "py_pvt.h"

namespace blue_sky
{
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("bs_pvt", "1.0.0", "BS_PVT", "BS_PVT", "bs_pvt")

    BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    const plugin_descriptor & pd = *bs_init.pd_;

    bool res = true;

    res &= blue_sky::pvt_register_types (pd); BS_ASSERT (res);

    return res;
  }
}//bs

#ifdef BSPY_EXPORTING_PLUGIN
BLUE_SKY_INIT_PY_FUN
{
  blue_sky::python::py_export_pvt ();
}
#endif //BSPY_EXPORT_PLUGIN
