// bs_bos_core_data_storage.cpp : Defines the entry point for the DLL application.
//

#include "bs_bos_core_data_storage_stdafx.h"

#include "data_manager.h"
#include "data_class.h"
#include "read_class.h"

#include "py_data_manager.h"
#include "py_data_class.h"
#include "py_bs_pool.h"
#include "py_read_class.h"

using namespace blue_sky;

namespace blue_sky
{
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("bs_bos_core_data_storage", "1.0.0", "BS_BOS_CORE_DATA_STORAGE", "BS_BOS_CORE_DATA_STORAGE", "bs_bos_core_data_storage")

  namespace
  {
    bool
    register_types (plugin_descriptor const &pd)
    {
      bool res = true;

      res &= BS_KERNEL.register_type (pd, data_manager::bs_type ()); BS_ASSERT (res);

      res &= BS_KERNEL.register_type (pd, idata::bs_type ()); BS_ASSERT (res);
      res &= BS_KERNEL.register_type(pd, FRead::bs_type()); BS_ASSERT (res);

      return res;
    }
  }

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    return register_types (*bs_init.pd_);
  }
}

#ifdef BSPY_EXPORTING_PLUGIN
namespace
{
  void
  init_py_subsystem ()
  {
    using namespace blue_sky;

    python::export_FRead ();
    //python::py_export_array_maps ();
    python::py_export_data_manager ();
    python::py_export_idata ();
  }
}
BLUE_SKY_INIT_PY_FUN
{
  init_py_subsystem ();
}
#ifdef _DEBUG
BOOST_PYTHON_MODULE (bs_bos_core_data_storage_d)
#else
BOOST_PYTHON_MODULE (bs_bos_core_data_storage)
#endif
{
  init_py_subsystem ();
  if (!blue_sky::register_types (*blue_sky::bs_get_plugin_descriptor ()))
    bs_throw_exception ("Can't register types");
}

#endif //BSPY_EXPORT_PLUGIN
