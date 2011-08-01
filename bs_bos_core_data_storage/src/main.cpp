// bs_bos_core_data_storage.cpp : Defines the entry point for the DLL application.
//

#include "bs_bos_core_data_storage_stdafx.h"

#include "hdm.h"
#include "well.h"
#include "data_class.h"
#include "read_class.h"

#include "py_hdm.h"
#include "py_well.h"
#include "py_well_storage.h"
#include "py_keyword_manager.h"
#include "py_data_class.h"
#include "py_read_class.h"
#include "equil_keywords.hpp"
#include "explicit_keywords.hpp"

namespace blue_sky
{
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("bs_bos_core_data_storage", "1.0.0", "BS_BOS_CORE_DATA_STORAGE", "BS_BOS_CORE_DATA_STORAGE", "bs_bos_core_data_storage");

  namespace {
  bool
  register_types (const plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, hdm::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, well_obj::bs_type ()); BS_ASSERT (res);

    res &= BS_KERNEL.register_type (pd, idata::bs_type ()); BS_ASSERT (res);
    
    res &= BS_KERNEL.register_type(pd, FRead::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, equil_keywords::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, explicit_keywords::bs_type ()); BS_ASSERT (res);
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
    using namespace blue_sky;

    python::export_FRead ();
    python::py_export_hdm ();
    python::py_export_idata ();
    python::export_keyword_manager();
    python::py_export_well_storage();
    python::py_export_well();
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
  std::cout << &BS_KERNEL << std::endl;
  bool res = blue_sky::register_types (*blue_sky::bs_get_plugin_descriptor ());
  if (!res)
    throw "Can't register bs-bos-core data_storage types";
}
#endif //BSPY_EXPORT_PLUGIN
