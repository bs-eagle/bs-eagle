/** 
 * @file py_frac.cpp
 * @brief python wraper for well fracture storage
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-29
 */

#include "py_frac.h"
#include "frac.h"

using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //! export matrices to python
  void py_export_frac ()
  {
    using namespace boost::python;

    base_exporter <frac_iface, py_frac_exporter>::export_class ("frac_iface");

    class_exporter <frac, frac_iface, py_frac_exporter>::export_class ("frac");

  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
