/** 
 * @file py_well.cpp
 * @brief python wraper for well wellture storage
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-29
 */

#include "py_well.h"
#include "well.h"

using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //! export matrices to python
  void py_export_well ()
  {
    using namespace boost::python;

    base_exporter <well_iface, py_well_exporter>::export_class ("well_iface");

    class_exporter <well, well_iface, py_well_exporter>::export_class ("well");

  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
