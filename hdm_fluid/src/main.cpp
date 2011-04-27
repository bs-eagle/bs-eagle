/**
 * \file main.cpp
 * \brief
 * \author Sergey Miryanov
 * \date 23.03.2009
 * */

#include "pvt_dead_oil.h"
#include "fluids.h"

using namespace blue_sky;
using namespace blue_sky::python;
using namespace boost::python;

namespace blue_sky {
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("hdm_fluid", "1.0.0", "Fluid representation for hydrodynamic simulator", "Fluid representation for hydrodynamic simulator", "hdm_fluid");

  namespace {
  bool
  register_types (plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type(pd, pvt_dead_oil::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(pd, fluids::bs_type()); BS_ASSERT (res);

    return res;
  }
  }

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    return register_types (*bs_init.pd_);
  }
}
//#if 0
#ifdef BSPY_EXPORTING_PLUGIN
BLUE_SKY_INIT_PY_FUN
{
  using namespace boost::python;

}
#ifdef _DEBUG
BOOST_PYTHON_MODULE (hdm_fluid_d)
#else
BOOST_PYTHON_MODULE (hdm_fluid)
#endif
{
  bs_init_py_subsystem ();
  std::cout << &BS_KERNEL << std::endl;
  bool res = blue_sky::register_types (*blue_sky::bs_get_plugin_descriptor ());
  if (!res)
    throw "Can't register hdm_fluid types";
}
#endif

//#endif //0
