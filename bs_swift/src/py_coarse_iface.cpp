/**
* \file   py_matrix_iface.cpp
* \brief  Python wrapper for linear solvers
* \author Miryanov Sergey
* \date 2008-04-04
*/
//#include "amg_stdafx.h"
#include "py_coarse_iface.h"
#include "ruge_coarse.h"
#include "pmis_coarse.h"
#include "pmis2_coarse.h"
#include "cljp_coarse.h"
#include "cljp2_coarse.h"

using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //! export matrices to python
  void py_export_coarse ()
  {
    using namespace boost::python;

    base_exporter<amg_coarse_iface, py_coarse_exporter>::export_class ("amg_coarse_iface");
    class_exporter<ruge_coarse, amg_coarse_iface, py_coarse_exporter>::export_class ("ruge_coarse");
    class_exporter<pmis_coarse, amg_coarse_iface, py_coarse_exporter>::export_class ("pmis_coarse");
    class_exporter<pmis2_coarse, amg_coarse_iface, py_coarse_exporter>::export_class ("pmis2_coarse");
    class_exporter<cljp_coarse, amg_coarse_iface, py_coarse_exporter>::export_class ("cljp_coarse");
    class_exporter<cljp2_coarse, amg_coarse_iface, py_coarse_exporter>::export_class ("cljp2_coarse");

  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN

