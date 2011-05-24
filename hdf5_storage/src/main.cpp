#ifdef BSPY_EXPORTING_PLUGIN
#include <boost/python.hpp>
#endif

#include "bs_kernel.h"
#include <stdlib.h>
#include "bos_report.h"
#include "h5_pool.hpp"
#include "py_pool.h"

#include "hdf5_group_impl.hpp"

using namespace blue_sky;
using namespace blue_sky::python;
using namespace boost::python;

namespace blue_sky {

  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("bs_hdf5_storage", "1.0.0", "Blue Sky HDF5 storage plugin", "Blue Sky HDF5 storage plugin", "bs_hdf5_storage")

  namespace {
    bool
    register_types (plugin_descriptor const &pd)
    {
      bool res = true; 
      res &= BLUE_SKY_REGISTER_TYPE (pd, h5_pool); BS_ASSERT (res);
      res &= BLUE_SKY_REGISTER_TYPE (pd, hdf5_group_impl); BS_ASSERT (res);

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
      python::py_export_pool ();
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
