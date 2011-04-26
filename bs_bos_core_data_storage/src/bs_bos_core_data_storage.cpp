// bs_bos_core_data_storage.cpp : Defines the entry point for the DLL application.
//

#include "bs_bos_core_data_storage_stdafx.h"

#include "hydrodynamic_model.h"
#include "data_class.h"
#include "read_class.h"

#include "py_hydrodynamic_model.h"
#include "py_keyword_manager.h"
#include "py_data_class.h"
#include "py_read_class.h"

namespace blue_sky
{
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("bs_bos_core_data_storage", "1.0.0", "BS_BOS_CORE_DATA_STORAGE", "BS_BOS_CORE_DATA_STORAGE", "bs_bos_core_data_storage");

  bool
  register_types (const plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, hydrodynamic_model::bs_type ()); BS_ASSERT (res);

    res &= BS_KERNEL.register_type (pd, idata::bs_type ()); BS_ASSERT (res);
    
    res &= BS_KERNEL.register_type(pd, FRead::bs_type()); BS_ASSERT (res);

    return res;
  }

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    return register_types (*bs_init.pd_);
  }
}

#ifdef BSPY_EXPORTING_PLUGIN
BLUE_SKY_INIT_PY_FUN
{
  using namespace blue_sky;

  python::export_FRead ();
  python::py_export_hydrodynamic_model ();
  python::py_export_idata ();
  python::export_keyword_manager();
}
#ifdef _DEBUG
BOOST_PYTHON_MODULE (bs_bos_core_data_storage_d)
#else
BOOST_PYTHON_MODULE (bs_bos_core_data_storage)
#endif
{
  bs_init_py_subsystem ();
  std::cout << &BS_KERNEL << std::endl;
  bool res = blue_sky::register_types (*blue_sky::bs_get_plugin_descriptor ());
  if (!res)
    throw "Can't register bs-bos-core data_storage types";
}
#endif //BSPY_EXPORT_PLUGIN
