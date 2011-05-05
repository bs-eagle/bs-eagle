/**
 * @file main.cpp
 * @brief
 * @author Oleg Borschuk
 * @date 2009-08-28
 */

#include "prop.h"
#include "table.h"
#include "py_prop.h"
#include "py_table.h"
#include "vartype_table.h"

using namespace blue_sky;
using namespace blue_sky::python;
using namespace boost::python;

#define REG_TYPE(S)                     \
    res &= BS_KERNEL.register_type(pd, S::bs_type()); BS_ASSERT (res);

namespace blue_sky {
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("comm", "1.0.0", "Common data types", "Common data types", "comm")

  namespace {
  bool
  register_types (const plugin_descriptor &pd)
  {
    bool res = true;

    REG_TYPE (prop)
    REG_TYPE (table)
    res &= BS_KERNEL.register_type (pd, vartype_table <t_float>::bs_type ());

    return res;
  }
  }

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    return register_types (*bs_init.pd_);
  }
}

#ifdef BSPY_EXPORTING_PLUGIN
namespace {
  void
  init_py_subsystem ()
  {
    using namespace boost::python;

    python::py_export_prop ();
    python::py_export_table ();
    // FIXME: export vartype table
  }
}
BLUE_SKY_INIT_PY_FUN
{
  init_py_subsystem ();
}
#ifdef _DEBUG
BOOST_PYTHON_MODULE (common_types_d)
#else
BOOST_PYTHON_MODULE (common_types)
#endif
{
  init_py_subsystem ();
  std::cout << &BS_KERNEL << std::endl;
  bool res = register_types (*blue_sky::bs_get_plugin_descriptor ());
  if (!res)
    throw "Can't register common types";

}
#endif