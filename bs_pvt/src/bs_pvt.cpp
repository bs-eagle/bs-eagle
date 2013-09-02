// bs_pvt.cpp : Defines the entry point for the DLL application.
//

#include "bs_pvt_stdafx.h"

#include "pvt_base.h"
#include "pvt_3p_iface.h"
#include "pvt_3p.h"
#include "py_pvt.h"
#include "pvt_keywords.hpp"
#include "bs_serialize_decl.h"

namespace blue_sky
{
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("bs_pvt", "1.0.0", "BS_PVT", "BS_PVT", "bs_pvt");

  namespace {
  bool
  register_types (const plugin_descriptor &pd)
  {
    bool res = true;

    res &= blue_sky::pvt_register_types (pd); BS_ASSERT (res); 
    res &= BS_KERNEL.register_type (pd, blue_sky::pvt_keywords::bs_type ()); BS_ASSERT (res);

    // force serialization typeinfo registering
    //serialize_register_eti< pvt_3p >();

    return res;
  }
  }

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    return register_types (*bs_init.pd_);
  }
}//bs

#ifdef BSPY_EXPORTING_PLUGIN
namespace {
  void
  init_py_subsystem ()
  {
    blue_sky::python::py_export_pvt ();
  }
}
BLUE_SKY_INIT_PY_FUN
{
  init_py_subsystem ();
}
#ifdef _DEBUG
BOOST_PYTHON_MODULE (bs_pvt_d)
#else
BOOST_PYTHON_MODULE (bs_pvt)
#endif
{
  init_py_subsystem ();
  std::cout << &BS_KERNEL << std::endl;
  bool res = blue_sky::register_types (*blue_sky::bs_get_plugin_descriptor ());
  if (!res)
    throw "Can't register pvt types";
}
#endif //BSPY_EXPORT_PLUGIN
