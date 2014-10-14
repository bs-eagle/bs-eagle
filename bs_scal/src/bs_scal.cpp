// bs_scal.cpp : Defines the entry point for the DLL application.
//

#include "scal_3p.h"
#include "py_scal_wrapper.h"
#include "scal_keywords.hpp"
#include "bs_serialize_decl.h"

namespace blue_sky
{
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("bs_scal", "1.0.0", "BS_SCAL", "BS_SCAL", "bs_scal");

  namespace {
  bool
  register_types (const plugin_descriptor &pd)
  {
    bool res = true;

    res &= blue_sky::scal_register_types (pd); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, scal_keywords::bs_type ()); BS_ASSERT (res);

    // force serialization typeinfo registering
    //serialize_register_eti< scal_3p >();

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
    blue_sky::python::py_export_scal ();
  }
}
BLUE_SKY_INIT_PY_FUN
{
  init_py_subsystem ();
}
#ifdef _DEBUG
BOOST_PYTHON_MODULE (bs_scal_d)
#else
BOOST_PYTHON_MODULE (bs_scal)
#endif
{
  init_py_subsystem ();
  std::cout << &BS_KERNEL << std::endl;
  bool res = blue_sky::register_types (*blue_sky::bs_get_plugin_descriptor ());
  if (!res)
    throw "Can't register scal types";
}
#endif //BSPY_EXPORT_PLUGIN

