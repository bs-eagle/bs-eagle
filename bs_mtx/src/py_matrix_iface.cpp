/**
* \file   py_matrix_iface.cpp
* \brief  Python wrapper for linear solvers
* \author Miryanov Sergey
* \date 2008-04-04
*/
#ifdef BSPY_EXPORTING_PLUGIN
#include <boost/python.hpp>
#endif
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

    base_exporter<matrix_iface, py_matrix_iface_exporter>::export_class ("matrix_iface");
    class_exporter<bdiag_matrix_iface, matrix_iface, py_bdiag_matrix_iface_exporter>::export_class  ("bdiag_matrix_iface");
    class_exporter<bdiag_matrix, bdiag_matrix_iface, py_bdiag_matrix_iface_exporter>::export_class ("bdiag_matrix");

  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN

