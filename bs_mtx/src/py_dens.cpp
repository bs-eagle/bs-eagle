/** 
 * @file py_dens.cpp
 * @brief 
 * @date 2009-12-09
 */
#ifdef BSPY_EXPORTING_PLUGIN
#include <boost/python.hpp>
#endif

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


    class_exporter<dens_matrix_iface, matrix_iface, py_dens_matrix_iface_exporter>::export_class ("dens_matrix_iface");
    class_exporter<dens_matrix, dens_matrix_iface, py_dens_matrix_iface_exporter>::export_class ("dens_matrix");
    base_exporter<dens_matrix_tools, py_matrix_dens_tools_exporter>::export_class ("dens_matrix_tools");
  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
