#include "bs_mesh_stdafx.h"

#include "py_mesh_grdecl.h"
#include "export_python_wrapper.h"

#ifdef BSPY_EXPORTING_PLUGIN
using namespace boost::python;

namespace blue_sky {
namespace python {

    void py_export_mesh_grdecl ()
    {
      strategy_exporter::export_class <bs_mesh_grdecl, rs_mesh_iface, mesh_grdecl_exporter> ("mesh_grdecl");
    }

} // namespace python
} // namespace blue_sky
#endif