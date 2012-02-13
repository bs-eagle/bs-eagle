/**
 * @file main.cpp
 * @brief
 * @author Oleg Borschuk
 * @date 2009-08-28
 */
#ifdef BSPY_EXPORTING_PLUGIN
#include <boost/python.hpp>
#endif

#include "prop.h"
#include "table.h"
#include "gis.h"
#include "frac.h"
#include "perf.h"
#include "traj.h"
#include "py_prop.h"
#include "py_table.h"
#include "py_gis.h"
#include "py_frac.h"
#include "py_perf.h"
#include "py_traj.h"
#include "vartype_table.h"
#include "bs_serialize_decl.h"

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
    REG_TYPE (gis)
    REG_TYPE (frac)
    REG_TYPE (perf)
    REG_TYPE (traj)
    res &= BS_KERNEL.register_type (pd, vartype_table <t_float>::bs_type ());

    // force serialization typeinfo registering
    //serialize_register_eti< table >();
    //serialize_register_eti< prop >();
    //serialize_register_eti< vartype_table< t_float > >();

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
    python::py_export_gis ();
    python::py_export_frac ();
    python::py_export_perf ();
    python::py_export_traj ();
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
