/** 
 * @file py_dens.cpp
 * @brief 
 * @date 2009-12-09
 */

#include "bs_matrix_stdafx.h"
#include "py_dens.h"

using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //! export matrices to python
  void py_export_dens_matrices ()
  {
    using namespace boost::python;


    strategy_exporter::export_class <dens_matrix_iface, matrix_iface, py_dens_matrix_iface_exporter> ("dens_matrix_iface");
    strategy_exporter::export_class <dens_matrix, dens_matrix_iface, py_dens_matrix_iface_exporter> ("dens_matrix");
    strategy_exporter::export_base  <dens_matrix_tools, py_matrix_dens_tools_exporter> ("dens_matrix_tools");
  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
