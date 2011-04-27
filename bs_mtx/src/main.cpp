/**
 * \file main.cpp
 * \brief
 * \author Sergey Miryanov
 * \date 23.03.2009
 * */

#include "bcsr.h"
#include "bdiag_matrix.h"
#include "jac_matrix.h"
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

namespace blue_sky {
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("mx", "1.0.0", "Base matrixes for blue_sky", "Base matrixes for blue_sky", "mx")

  namespace {
  bool
  register_types (const plugin_descriptor &pd)
  {
    //const plugin_descriptor &pd = *bs_init.pd_;

    bool res = true;

    res &= BS_KERNEL.register_type (pd, dens_matrix::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, dens_matrix_tools::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, bcsr::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, bcsr_matrix_tools::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, bdiag_matrix::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, jac_matrix::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, mbcsr_matrix::bs_type()); BS_ASSERT (res);

    return res;
  }
  }
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
  bool res = blue_sky::register_types (*blue_sky::bs_get_plugin_descriptor ());
  if (!res)
    throw "Can't register mtx types";
}
#endif

//#endif //0
