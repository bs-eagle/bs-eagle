/**
 * \file bs_bos_core_base.cpp
 * \brief
 * \author Sergey Miryanov
 * \date 24.03.2009
 * */
#include "bs_bos_core_base_stdafx.h"

#include "py_bs_assert.h"
#include "named_pbase_access.h"
#include "data_dimens.h"

using namespace blue_sky;

#ifdef BSPY_EXPORTING_PLUGIN
using namespace blue_sky::python;
using namespace boost::python;
#endif

namespace blue_sky {
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("bs_bos_core_base", "1.0.0", "BS_BOS_CORE_BASE", "BS_BOS_CORE_BASE", "bs_bos_core_base")

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    init_bos_logs();
   
    const plugin_descriptor &pd = *bs_init.pd_;

    bool res = true;
    res &= BS_KERNEL.register_type (pd, property_base::bs_type ());
    res &= BS_KERNEL.register_type (pd, named_pbase::bs_type ());

    return res;
  }
}

#ifdef BSPY_EXPORTING_PLUGIN
BLUE_SKY_INIT_PY_FUN
{
  using namespace boost::python;

  //python::py_export_assert ();
  python::py_export_named_pbase ("named_pbase");
  python::export_data_dimens ();
}
#endif 
