/**
 * \file main.cpp
 * \brief
 * \author Sergey Miryanov
 * \date 23.03.2009
 * */
#include "bs_matrix_stdafx.h"

#include "bcsr.h"
#include "bdiag_matrix.h"
#include "jac_matrix.h"
#include "py_matrix_iface.h"
#include "py_bcsr.h"
#include "py_dens.h"
#include "bcsr_matrix_tools.h"
#include "dens_matrix_tools.h"
#include "dens_matrix.h"
//#include "bsvector.h"
//#include "py_bsvector.h"

#include "py_naive_file_reader.h"
#include "py_save_seq_vector.h"

using namespace blue_sky;
using namespace blue_sky::python;
using namespace boost::python;

namespace blue_sky {
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("mx", "1.0.0", "Base matrixes for blue_sky", "Base matrixes for blue_sky", "mx")

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    //const plugin_descriptor &pd = *bs_init.pd_;

    bool res = true;

    //res &= BS_KERNEL.register_type(*bs_init.pd_, ndarray<float>::bs_type()); BS_ASSERT (res);
    //res &= BS_KERNEL.register_type(*bs_init.pd_, ndarray<double>::bs_type()); BS_ASSERT (res);
    //res &= BS_KERNEL.register_type(*bs_init.pd_, ndarray<int>::bs_type()); BS_ASSERT (res);
    //res &= BS_KERNEL.register_type(*bs_init.pd_, ndarray<long>::bs_type()); BS_ASSERT (res);
    //res &= BS_KERNEL.register_type(*bs_init.pd_, ndarray<char>::bs_type()); BS_ASSERT (res);

    //res &= BS_KERNEL.register_type(*bs_init.pd_, bsvector<float>::bs_type()); BS_ASSERT (res);
    //res &= BS_KERNEL.register_type(*bs_init.pd_, bsvector<double>::bs_type()); BS_ASSERT (res);
    //res &= BS_KERNEL.register_type(*bs_init.pd_, bsvector<int>::bs_type()); BS_ASSERT (res);
    //res &= BS_KERNEL.register_type(*bs_init.pd_, bsvector<long>::bs_type()); BS_ASSERT (res);
    //res &= BS_KERNEL.register_type(*bs_init.pd_, bsvector<char>::bs_type()); BS_ASSERT (res);

    //res &= BS_KERNEL.register_type(*bs_init.pd_, bsmatrix<float>::bs_type()); BS_ASSERT (res);
    //res &= BS_KERNEL.register_type(*bs_init.pd_, bsmatrix<double>::bs_type()); BS_ASSERT (res);
    //res &= BS_KERNEL.register_type(*bs_init.pd_, bsmatrix<int>::bs_type()); BS_ASSERT (res);
    //res &= BS_KERNEL.register_type(*bs_init.pd_, bsmatrix<long>::bs_type()); BS_ASSERT (res);
    //res &= BS_KERNEL.register_type(*bs_init.pd_, bsmatrix<char>::bs_type()); BS_ASSERT (res);

    res &= BS_KERNEL.register_type(*bs_init.pd_, dens_matrix<base_strategy_fif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, dens_matrix<base_strategy_did>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, dens_matrix<base_strategy_dif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, dens_matrix<base_strategy_flf>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, dens_matrix<base_strategy_dld>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, dens_matrix<base_strategy_dlf>::bs_type()); BS_ASSERT (res);

    res &= BS_KERNEL.register_type(*bs_init.pd_, dens_matrix_tools<base_strategy_fif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, dens_matrix_tools<base_strategy_did>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, dens_matrix_tools<base_strategy_dif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, dens_matrix_tools<base_strategy_flf>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, dens_matrix_tools<base_strategy_dld>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, dens_matrix_tools<base_strategy_dlf>::bs_type()); BS_ASSERT (res);

    res &= BS_KERNEL.register_type(*bs_init.pd_, bcsr<base_strategy_fif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bcsr<base_strategy_did>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bcsr<base_strategy_dif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bcsr<base_strategy_flf>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bcsr<base_strategy_dld>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bcsr<base_strategy_dlf>::bs_type()); BS_ASSERT (res);

    res &= BS_KERNEL.register_type(*bs_init.pd_, bcsr_matrix_tools<base_strategy_fif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bcsr_matrix_tools<base_strategy_did>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bcsr_matrix_tools<base_strategy_dif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bcsr_matrix_tools<base_strategy_flf>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bcsr_matrix_tools<base_strategy_dld>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bcsr_matrix_tools<base_strategy_dlf>::bs_type()); BS_ASSERT (res);
    
    res &= BS_KERNEL.register_type(*bs_init.pd_, bdiag_matrix<base_strategy_fif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bdiag_matrix<base_strategy_did>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bdiag_matrix<base_strategy_dif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bdiag_matrix<base_strategy_flf>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bdiag_matrix<base_strategy_dld>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, bdiag_matrix<base_strategy_dlf>::bs_type()); BS_ASSERT (res);

    res &= BS_KERNEL.register_type(*bs_init.pd_, jac_matrix<base_strategy_fif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, jac_matrix<base_strategy_did>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, jac_matrix<base_strategy_dif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, jac_matrix<base_strategy_flf>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, jac_matrix<base_strategy_dld>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(*bs_init.pd_, jac_matrix<base_strategy_dlf>::bs_type()); BS_ASSERT (res);

    //res &= blue_sky::jacobian_matrix_register_type (pd); BS_ASSERT (res);

    return res;
  }
}
//#if 0
#ifdef BSPY_EXPORTING_PLUGIN
BLUE_SKY_INIT_PY_FUN
{
  using namespace boost::python;

  //Python vector types
  python::py_export_naive_file_reader ();
  python::py_export_save_shared_vector ();

  python::py_export_matrices ();
  python::py_export_bcsr_matrices ();
  python::py_export_dens_matrices ();
  //python::py_export_bsvector ();
  //python::py_export_jacobian_matrix ();

  //python::py_export_read_b_matrix ();

  
}
#endif

//#endif //0
