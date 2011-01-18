/**
* \file   py_bcsr.cpp
* \brief  Python wrapper for BCSR matrices
* \date 2008-04-04
*/
#include "bs_mtx_stdafx.h"
#include "py_bcsr.h"
#include "bs_object_base.h"

using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //! export matrices to python
  void py_export_bcsr_matrices ()
  {
    using namespace boost::python;


    strategy_exporter::export_class <bcsr_matrix_iface, matrix_iface, py_bcsr_matrix_iface_exporter> ("bcsr_matrix_iface");
    strategy_exporter::export_class <bcsr_amg_matrix_iface, bcsr_matrix_iface, py_bcsr_amg_matrix_iface_exporter> ("bcsr_amg_matrix_iface");
    strategy_exporter::export_class <bcsr, bcsr_amg_matrix_iface, py_bcsr_amg_matrix_iface_exporter> ("bcsr_matrix");
    strategy_exporter::export_base <bcsr_matrix_tools, py_matrix_bcsr_tools_exporter> ("bcsr_matrix_tools");

  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
