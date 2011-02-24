/**
* \file   py_bcsr.cpp
* \brief  Python wrapper for BCSR matrices
* \date 2008-04-04
*/
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


    class_exporter<bcsr_matrix_iface, matrix_iface, py_bcsr_matrix_iface_exporter>::export_class ("bcsr_matrix_iface");
    class_exporter<bcsr_amg_matrix_iface, bcsr_matrix_iface, py_bcsr_amg_matrix_iface_exporter>::export_class ("bcsr_amg_matrix_iface");
    class_exporter<bcsr, bcsr_amg_matrix_iface, py_bcsr_amg_matrix_iface_exporter>::export_class ("bcsr_matrix");
    base_exporter<bcsr_matrix_tools, py_matrix_bcsr_tools_exporter>::export_class ("bcsr_matrix_tools");

  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
