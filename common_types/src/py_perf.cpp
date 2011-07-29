/** 
 * @file py_perf.cpp
 * @brief python wraper for well perf storage
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-29
 */

#include "py_perf.h"
#include "perf.h"

using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //! export matrices to python
  void py_export_perf ()
  {
    using namespace boost::python;

    base_exporter <perf_iface, py_perf_exporter>::export_class ("perf_iface");

    class_exporter <perf, perf_iface, py_perf_exporter>::export_class ("perf");

  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
