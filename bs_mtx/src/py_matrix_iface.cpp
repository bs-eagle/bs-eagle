/**
* \file   py_matrix_iface.cpp
* \brief  Python wrapper for linear solvers
* \author Miryanov Sergey
* \date 2008-04-04
*/
#include "bs_matrix_stdafx.h"
#include "py_matrix_iface.h"

using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //! export matrices to python
  void py_export_matrices ()
  {
    using namespace boost::python;

    strategy_exporter::export_base <matrix_iface, py_matrix_iface_exporter> ("matrix_iface");
    strategy_exporter::export_class <bdiag_matrix_iface, matrix_iface, py_bdiag_matrix_iface_exporter> ("bdiag_matrix_iface");
    strategy_exporter::export_class <jac_matrix_iface, matrix_iface, py_jac_matrix_iface_exporter> ("jac_matrix_iface");
    strategy_exporter::export_class <bdiag_matrix, bdiag_matrix_iface, py_bdiag_matrix_iface_exporter> ("bdiag_matrix");
    strategy_exporter::export_class <jac_matrix, jac_matrix_iface, py_jac_matrix_iface_exporter> ("jac_matrix");

  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN

