/**
 * @file main.cpp
 * @brief
 * @author Oleg Borschuk
 * @date 2009-08-28
 */
#ifdef BSPY_EXPORTING_PLUGIN
#include <boost/python.hpp>
#endif

#include "sql_well.h" 
#include "py_sql_well.h"
#include "well_keywords.hpp"
#include "bs_serialize_decl.h"

using namespace blue_sky;
using namespace blue_sky::python;
using namespace boost::python;

#define REG_TYPE(S)                     \
    res &= BS_KERNEL.register_type(pd, S::bs_type()); BS_ASSERT (res);

namespace blue_sky {
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("sql_well", "1.0.0", "Well storage using sqlite", "Well storage using sqlite", "sql_well")

  namespace {
  bool
  register_types (const plugin_descriptor &pd)
  {
    bool res = true;
    setlocale(LC_NUMERIC, "C");

    REG_TYPE (sql_well);
    res &= BS_KERNEL.register_type (pd, blue_sky::well_keywords::bs_type ()); BS_ASSERT (res);
//res &= BS_KERNEL.register_type (pd, vartype_table <t_float>::bs_type ());
    
    // force serialization typeinfo registering
    //serialize_register_eti< sql_well >();

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
    py_export_sql_well ();
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
