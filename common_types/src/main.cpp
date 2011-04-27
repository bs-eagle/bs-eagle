/** 
 * @file main.cpp
 * @brief 
 * @author Oleg Borschuk
 * @date 2009-08-28
 */

#include "prop.h"
#include "table.h"
#include "h5_pool.h"
#include "py_prop.h"
#include "py_table.h"
#include "py_pool.h"

using namespace blue_sky;
using namespace blue_sky::python;
using namespace boost::python;

namespace blue_sky {
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("comm", "1.0.0", "Common data types", "Common data types", "comm")

  namespace {
  bool 
  register_types (const plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, prop::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, table::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, h5_pool::bs_type()); BS_ASSERT (res);

    return res;
  }
  }

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    return register_types (*bs_init.pd_);
  }
}

#ifdef BSPY_EXPORTING_PLUGIN
BLUE_SKY_INIT_PY_FUN
{
  using namespace boost::python;

  python::py_export_prop ();
  python::py_export_table ();
  python::py_export_pool ();
}
#ifdef _DEBUG
BOOST_PYTHON_MODULE (common_types_d)
#else
BOOST_PYTHON_MODULE (common_types)
#endif
{
  bs_init_py_subsystem ();
  std::cout << &BS_KERNEL << std::endl;
  bool res = blue_sky::register_types (*blue_sky::bs_get_plugin_descriptor ());
  if (!res)
    throw "Can't register common types";

}
#endif
