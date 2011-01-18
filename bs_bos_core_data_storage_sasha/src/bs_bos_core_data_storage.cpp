// bs_bos_core_data_storage.cpp : Defines the entry point for the DLL application.
//

#include "bs_bos_core_data_storage_stdafx.h"

#include "data_manager.h"
#include "data_class.h"
#include "read_class.h"
#include "pool_treeish.h"

#include "py_data_manager.h"
#include "py_data_class.h"
#include "py_bs_pool.h"
#include "py_pool_treeish.h"

namespace blue_sky
{
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("bs_bos_core_data_storage", "1.0.0", "BS_BOS_CORE_DATA_STORAGE", "BS_BOS_CORE_DATA_STORAGE", "bs_bos_core_data_storage")

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    //bool res = true;
    const plugin_descriptor & pd = *bs_init.pd_;

    bool res = true;

    res &= BS_KERNEL.register_type (pd, data_manager <base_strategy_fif>::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, data_manager <base_strategy_did>::bs_type ()); BS_ASSERT (res);

    res &= BS_KERNEL.register_type (pd, idata <base_strategy_fif>::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, idata <base_strategy_did>::bs_type ()); BS_ASSERT (res);

    res &= BS_KERNEL.register_type(pd, FRead::bs_type()); BS_ASSERT (res);

	res &= register_bs_pool_node(pd);
    return res;
  }
}

#ifdef BSPY_EXPORTING_PLUGIN
BLUE_SKY_INIT_PY_FUN
{
  using namespace blue_sky;

  //python::py_export_array_maps ();
  python::py_export_data_manager ();
  python::py_export_idata ();
  python::py_export_bs_pool_node();

}
#endif //BSPY_EXPORT_PLUGIN
