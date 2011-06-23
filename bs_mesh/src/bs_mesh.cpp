// bs_mesh.cpp : Defines the entry point for the DLL application.
//

#include "bs_mesh_stdafx.h"

#include "bs_mesh_ijk.h"
#include "bs_mesh_grdecl.h"
#include "bs_flux_connections.h"

#include "mesh_ijk_keywords.h"
#include "mesh_grdecl_keywords.h"

#include "py_rs_mesh.h"
#include "py_mesh_grdecl.h"
#include "py_flux_connections.h"

namespace blue_sky
{
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("bs_mesh", "1.0.0", "BS_MESH", "BS_MESH", "bs_mesh")

  namespace 
  {
    bool
    register_types (const plugin_descriptor &pd)
    {
      bool res = true;

      //res &= BS_KERNEL.register_type(*bs_init.pd_, bs_flux_connections::bs_type()); BS_ASSERT (res);
      //res &= blue_sky::give_kernel::Instance().register_type(*bs_init.pd_, mesh_rs::bs_type()); BS_ASSERT (res);
      //res &= blue_sky::give_kernel::Instance().register_type(*bs_init.pd_, mesh_ijk::bs_type()); BS_ASSERT (res);
      
      res &= BS_KERNEL.register_type(pd, mesh_ijk_keywords::bs_type()); BS_ASSERT (res);
      res &= BS_KERNEL.register_type(pd, bs_mesh_grdecl::bs_type()); BS_ASSERT (res);
      res &= BS_KERNEL.register_type(pd, mesh_grdecl_keywords::bs_type()); BS_ASSERT (res);
      
      //mpi_mesh_grdecl
#ifdef _MPI_MY
      res &= BS_KERNEL.register_type(pd, mpi_mesh_grdecl::bs_type()); BS_ASSERT (res);
#endif

      // flux_connections
      res &= BS_KERNEL.register_type(pd, bs_flux_connections::bs_type()); BS_ASSERT (res);
      
      return res;
    }
  }

  BLUE_SKY_REGISTER_PLUGIN_FUN
  {
    return register_types (*bs_init.pd_);
  }
}//bs

#ifdef BSPY_EXPORTING_PLUGIN
namespace 
{
  void
  init_py_subsystem ()
  {
    using namespace blue_sky;

    python::py_export_mesh ();
    python::py_export_mesh_grdecl ();
    //python::py_export_flux_connections ();
  }
}
BLUE_SKY_INIT_PY_FUN
{
  init_py_subsystem ();
}

#ifdef _DEBUG
BOOST_PYTHON_MODULE (bs_mesh_d)
#else
BOOST_PYTHON_MODULE (bs_mesh)
#endif
{
  init_py_subsystem ();
  if (!blue_sky::register_types (*blue_sky::bs_get_plugin_descriptor ()))
    bs_throw_exception ("Can't register types");
}

#endif //BSPY_EXPORT_PLUGIN

