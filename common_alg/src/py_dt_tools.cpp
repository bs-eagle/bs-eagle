/** 
 * @file py_bos_reader.cpp
 * @brief python wraper for BOS ascii file reader
 * @author Oleg Borschuk
 * @version 
 * @date 2012-03-01
 */

#include "py_dt_tools.h"
#include "dt_tools.h"

using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {

  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //! export matrices to python
  void py_export_dt_tools ()
  {
    using namespace boost::python;

    base_exporter <dt_tools_iface, py_dt_tools_exporter>::export_class ("dt_tools_iface");
    class_exporter <dt_tools, dt_tools_iface, py_dt_tools_exporter>::export_class ("dt_tools");


    // register implicit conversion to interface
    implicitly_convertible<
      smart_ptr< dt_tools >,
      smart_ptr< dt_tools_iface >
    >();
  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
