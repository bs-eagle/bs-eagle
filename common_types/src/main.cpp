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

using namespace blue_sky;
using namespace blue_sky::python;
using namespace boost::python;

namespace blue_sky {
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("common_types", "1.0.0", "Common data types", "Common data types", "common_types")

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    //const plugin_descriptor &pd = *bs_init.pd_;

    bool res = true;

    res &= BS_KERNEL.register_type(*bs_init.pd_, prop<float, int, std::string, bool>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, table::bs_type()); BS_ASSERT (res);

    return res;
  }
}

#ifdef BSPY_EXPORTING_PLUGIN
BLUE_SKY_INIT_PY_FUN
{
  using namespace boost::python;

  python::py_export_prop ();
  python::py_export_table ();
}
#endif
