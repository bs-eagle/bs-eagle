/** 
 * @file py_upsc.cpp
 * @brief python wraper for hdm upscaling
 * @author ALina Yapparova
 * @version 
 * @date 2012-20-02
 */
#include "bs_mesh_stdafx.h"
#include "py_upsc.h"
#include "upsc.h"

using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //! export matrices to python
  void py_export_upsc ()
  {
    using namespace boost::python;

    base_exporter <upsc_iface, py_upsc_exporter>::export_class ("upsc_iface");

    class_exporter <upsc, upsc_iface, py_upsc_exporter>::export_class ("upsc");

  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
