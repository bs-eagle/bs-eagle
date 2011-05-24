/** 
 * @file py_mbcsr_matrix.cpp
 * @brief python wrapper for #mbcsr_matrix
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-03-02
 */
#ifdef BSPY_EXPORTING_PLUGIN
#include <boost/python.hpp>
#endif

#include "py_mbcsr_matrix.h"
#include "bs_object_base.h"

using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //! export matrices to python
  void py_export_mbcsr_matrices ()
  {
    using namespace boost::python;


    class_exporter<mbcsr_matrix_iface, matrix_iface, py_mbcsr_matrix_iface_exporter>::export_class ("mbcsr_matrix_iface");
    class_exporter<mbcsr_matrix, mbcsr_matrix_iface, py_mbcsr_matrix_iface_exporter>::export_class ("mbcsr_matrix");

  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
