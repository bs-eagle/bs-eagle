// bs_scal.cpp : Defines the entry point for the DLL application.
//

#include "bs_scal_stdafx.h"
#include "scal_3p.h"
#include "py_scal_wrapper.h"

namespace blue_sky
{
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("bs_scal", "1.0.0", "BS_SCAL", "BS_SCAL", "bs_scal");

  namespace {
  bool
  register_types (const plugin_descriptor &pd)
  {
    bool res = true;

    res &= blue_sky::scal_register_types (pd); BS_ASSERT (res);

    return res;
  }
  }


  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    return register_types (*bs_init.pd_);
  }
}//bs

#ifdef BSPY_EXPORTING_PLUGIN
BLUE_SKY_INIT_PY_FUN
{
  blue_sky::python::py_export_scal ();
}
#ifdef _DEBUG
BOOST_PYTHON_MODULE (bs_scal_d)
#else
BOOST_PYTHON_MODULE (bs_scal)
#endif
{
  bs_init_py_subsystem ();
  std::cout << &BS_KERNEL << std::endl;
  bool res = blue_sky::scal_register_types (*blue_sky::bs_get_plugin_descriptor ());
  if (!res)
    throw "Can't register scal types";
}
#endif //BSPY_EXPORT_PLUGIN

