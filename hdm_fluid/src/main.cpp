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
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("hdm_fluid", "1.0.0", "Fluid representation for hydrodynamic simulator", "Fluid representation for hydrodynamic simulator", "hdm_fluid")

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    //const plugin_descriptor &pd = *bs_init.pd_;

    bool res = true;

    res &= BS_KERNEL.register_type(*bs_init.pd_, pvt_dead_oil::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, fluids::bs_type()); BS_ASSERT (res);

    return res;
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
}
#endif

//#endif //0
