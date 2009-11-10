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

namespace blue_sky
{
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("bs_bos_core_data_storage", "1.0.0", "BS_BOS_CORE_DATA_STORAGE", "BS_BOS_CORE_DATA_STORAGE", "bs_bos_core_data_storage")

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    //bool res = true;
    const plugin_descriptor & pd = *bs_init.pd_;

    bool res = true;

    res &= BS_KERNEL.register_type (pd, data_manager <base_strategy_fi>::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, data_manager <base_strategy_di>::bs_type ()); BS_ASSERT (res);

    res &= BS_KERNEL.register_type (pd, idata::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(pd, FRead::bs_type()); BS_ASSERT (res);

    return res;
  }
}

#ifdef BSPY_EXPORTING_PLUGIN
BLUE_SKY_INIT_PY_FUN
{
  using namespace blue_sky;

  python::export_FRead ();
  //python::py_export_array_maps ();
  python::py_export_data_manager ();
  python::py_export_idata ();
}
#endif //BSPY_EXPORT_PLUGIN
