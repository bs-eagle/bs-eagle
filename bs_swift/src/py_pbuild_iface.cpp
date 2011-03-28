/**
* \file   py_pbuild_iface.cpp
* \brief  Python wrapper for Prolangation matrix builder
*/
//#include "amg_stdafx.h"
#include "py_pbuild_iface.h"
#include "direct_pbuild.h"
#include "standart_pbuild.h"
#include "standart2_pbuild.h"
#include "standart3_pbuild.h"

using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //! export matrices to python
  void py_export_pbuild ()
  {
    using namespace boost::python;

    base_exporter<amg_pbuild_iface, py_pbuild_exporter>::export_class ("amg_pbuild_iface");
    class_exporter<direct_pbuild, amg_pbuild_iface, py_pbuild_exporter>::export_class ("direct_pbuild");
    class_exporter<standart_pbuild, amg_pbuild_iface, py_pbuild_exporter>::export_class ("standart_pbuild");
    class_exporter<standart2_pbuild, amg_pbuild_iface, py_pbuild_exporter>::export_class ("standart2_pbuild");
    class_exporter<standart3_pbuild, amg_pbuild_iface, py_pbuild_exporter>::export_class ("standart3_pbuild");

  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN

