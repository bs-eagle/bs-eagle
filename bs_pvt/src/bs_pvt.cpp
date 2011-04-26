// bs_pvt.cpp : Defines the entry point for the DLL application.
//

#include "bs_pvt_stdafx.h"

#include "pvt_base.h"
#include "py_pvt.h"

namespace blue_sky
{
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("bs_pvt", "1.0.0", "BS_PVT", "BS_PVT", "bs_pvt");

  bool
  register_types (const plugin_descriptor &pd)
  {
    bool res = true;

    res &= blue_sky::pvt_register_types (pd); BS_ASSERT (res);

    return res;
  }

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    return register_types (*bs_init.pd_);
  }
}//bs

#ifdef BSPY_EXPORTING_PLUGIN
BLUE_SKY_INIT_PY_FUN
{
  blue_sky::python::py_export_pvt ();
}
#ifdef _DEBUG
BOOST_PYTHON_MODULE (bs_pvt_d)
#else
BOOST_PYTHON_MODULE (bs_pvt)
#endif
{
  bs_init_py_subsystem ();
  std::cout << &BS_KERNEL << std::endl;
  bool res = blue_sky::pvt_register_types (*blue_sky::bs_get_plugin_descriptor ());
  if (!res)
    throw "Can't register pvt types";
}
#endif //BSPY_EXPORT_PLUGIN
