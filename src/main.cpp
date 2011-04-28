#include "bs_hdf5_storage.h"
#include "bs_kernel.h"
#include <stdlib.h>
#include "py_bs_hdf5_storage.h"
#include "bos_report.h"
#include "h5_pool_ex.hpp"

using namespace blue_sky;
using namespace blue_sky::python;
using namespace boost::python;

namespace blue_sky {

  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("bs_hdf5_storage", "1.0.0", "Blue Sky HDF5 storage plugin", "Blue Sky HDF5 storage plugin", "bs_hdf5_storage")

  namespace {
    bool
    register_types (plugin_descriptor const &pd)
    {
      bool res = BLUE_SKY_REGISTER_TYPE(pd, bs_hdf5_storage); BS_ASSERT (res);
      res &= BLUE_SKY_REGISTER_TYPE (pd, h5_pool_ex); BS_ASSERT (res);

      return res;
    }
  }

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    return register_types (*bs_init.pd_);
  }

#ifdef BSPY_EXPORTING_PLUGIN
  namespace {
    void
    init_py_subsystem ()
    {
      py_export_bs_hdf5_storage ();
    }
  }
  BLUE_SKY_INIT_PY_FUN
  {
    init_py_subsystem ();
  }
#ifdef _DEBUG
BOOST_PYTHON_MODULE (bs_hdf5_storage_d)
#else
BOOST_PYTHON_MODULE (bs_hdf5_storage)
#endif
{
  init_py_subsystem ();
  std::cout << &BS_KERNEL << std::endl;
  bool res = blue_sky::register_types (*blue_sky::bs_get_plugin_descriptor ());
  if (!res)
    throw "Can't register hdf5_storage types";
}
#endif //BSPY_EXPORTING_PLUGIN

}
