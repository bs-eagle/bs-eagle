/**
 * @file main.cpp
 * @brief
 * @author Oleg Borschuk
 * @date 2009-08-28
 */
#ifdef BSPY_EXPORTING_PLUGIN
#include <boost/python.hpp>
#endif
#include "bos_reader.h"
#include "dt_tools.h"
#include "py_bos_reader.h"
#include "py_dt_tools.h"
#include "pool_iface.h"

using namespace blue_sky;
using namespace blue_sky::python;
using namespace boost::python;

#define REG_TYPE(S)                     \
    res &= BS_KERNEL.register_type(pd, S::bs_type()); BS_ASSERT (res);

namespace blue_sky {
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("alg", "1.0.0", "Common algorithms", "Common algorithms", "alg")

  namespace {
  bool
  register_types (const plugin_descriptor &pd)
  {
    bool res = true;

    REG_TYPE (bos_reader)
    REG_TYPE (dt_tools)
    //REG_TYPE (prop)
    //REG_TYPE (table)
    //REG_TYPE (gis)
    //REG_TYPE (frac)
    //REG_TYPE (perf)
    //REG_TYPE (traj)
    //res &= BS_KERNEL.register_type (pd, vartype_table <t_float>::bs_type ());

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

  // forward declaration to export to python
  void read_grdecl(const std::string& fname, const std::string dir, ulong nx, ulong ny, ulong nz,
    smart_ptr< h5_pool_iface > pool);

  spv_float wlog_mean_projection(spv_float wlog_data, spv_float wlog_dept, spv_float dest_grid);

}

#ifdef BSPY_EXPORTING_PLUGIN
namespace {
  void
  init_py_subsystem ()
  {
    using namespace boost::python;

    python::py_export_bos_reader ();
    python::py_export_dt_tools ();

    def("read_grdecl", &blue_sky::read_grdecl);
    def("wlog_mean_projection", &blue_sky::wlog_mean_projection);
    //python::py_export_table ();
    //python::py_export_gis ();
    //python::py_export_frac ();
    //python::py_export_perf ();
    //python::py_export_traj ();
    // FIXME: export vartype table
  }
}
BLUE_SKY_INIT_PY_FUN
{
  init_py_subsystem ();
}
#ifdef _DEBUG
BOOST_PYTHON_MODULE (common_alg_d)
#else
BOOST_PYTHON_MODULE (common_alg)
#endif
{
  init_py_subsystem ();
  std::cout << &BS_KERNEL << std::endl;
  bool res = register_types (*blue_sky::bs_get_plugin_descriptor ());
  if (!res)
    throw "Can't register common types";

}
#endif
