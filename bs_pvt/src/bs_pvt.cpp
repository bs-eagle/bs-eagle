// bs_pvt.cpp : Defines the entry point for the DLL application.
//

#include "bs_pvt_stdafx.h"

#include "pvt_base.h"
#include "py_pvt.h"

using namespace blue_sky;

namespace blue_sky
{
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("bs_pvt", "1.0.0", "BS_PVT", "BS_PVT", "bs_pvt")

  namespace
  {
    bool
    register_types (plugin_descriptor const &pd)
    {
      bool res = true;

      res &= blue_sky::pvt_register_types (pd); BS_ASSERT (res);

      return res;
    }
  }

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    return register_types (*bs_init.pd_);
  }
}//bs

#ifdef BSPY_EXPORTING_PLUGIN
namespace 
{
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
  if (!blue_sky::register_types (*blue_sky::bs_get_plugin_descriptor ()))
    bs_throw_exception ("Can't register types");
}

#endif //BSPY_EXPORT_PLUGIN
