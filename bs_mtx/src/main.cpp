/**
 * \file main.cpp
 * \brief
 * \author Sergey Miryanov
 * \date 23.03.2009
 * */

#include "bcsr.h"
#include "bdiag_matrix.h"
#include "py_matrix_iface.h"
#include "py_bcsr.h"
#include "py_dens.h"
#include "bcsr_matrix_tools.h"
#include "dens_matrix_tools.h"
#include "dens_matrix.h"
#include "mbcsr_matrix.h"
#include "py_mbcsr_matrix.h"

using namespace blue_sky;
using namespace blue_sky::python;
using namespace boost::python;

#define REG_TYPE(S)                     \
    res &= BS_KERNEL.register_type(pd, S::bs_type()); BS_ASSERT (res);

namespace {
  using namespace blue_sky;

  bool
  register_types (const plugin_descriptor &pd)
  {
    bool res = true;

    REG_TYPE (dens_matrix)
    REG_TYPE (dens_matrix_tools)
    REG_TYPE (bcsr)
    REG_TYPE (bcsr_matrix_tools)
    REG_TYPE (bdiag_matrix)
    REG_TYPE (mbcsr_matrix)

    return res;
  }
}

namespace blue_sky {
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("mx", "1.0.0", "Base matrixes for blue_sky", "Base matrixes for blue_sky", "mx")

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    return register_types (*bs_init.pd_);
  }
}
//#if 0
#ifdef BSPY_EXPORTING_PLUGIN
BLUE_SKY_INIT_PY_FUN
{
  using namespace boost::python;

  python::py_export_matrices ();
  python::py_export_bcsr_matrices ();
  python::py_export_mbcsr_matrices ();
  python::py_export_dens_matrices ();

}

#ifdef _DEBUG
BOOST_PYTHON_MODULE (bs_mtx_d)
#else
BOOST_PYTHON_MODULE (bs_mtx)
#endif
{
  bs_init_py_subsystem ();
  std::cout << &BS_KERNEL << std::endl;
  bool res = register_types (*blue_sky::bs_get_plugin_descriptor ());
  if (!res)
    throw "Can't register mtx types";
}
#endif

//#endif //0
