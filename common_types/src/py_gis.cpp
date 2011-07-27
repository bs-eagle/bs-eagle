/** 
 * @file py_gis.cpp
 * @brief python interface to WELL GIS
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-26
 */

#include "py_gis.h"
#include "gis.h"

using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //! export matrices to python
  void py_export_gis ()
  {
    using namespace boost::python;

    base_exporter <gis_iface, py_gis_exporter>::export_class ("gis_iface");

    class_exporter <gis, gis_iface, py_gis_exporter>::export_class ("gis");

  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
