// bs_mesh.cpp : Defines the entry point for the DLL application.
//

#include "bs_mesh_stdafx.h"

#include "bs_mesh_ijk.h"
#include "bs_mesh_grdecl.h"
#include "bs_flux_connections.h"

#include "mesh_ijk_keywords.h"
#include "mesh_grdecl_keywords.h"

#include "py_rs_mesh.h"

namespace blue_sky
{
  BLUE_SKY_PLUGIN_DESCRIPTOR_EXT ("bs_mesh", "1.0.0", "BS_MESH", "BS_MESH", "bs_mesh");

  namespace {
  bool
  register_types (const plugin_descriptor &pd)
  {
    bool res = true;

    //res &= BS_KERNEL.register_type(*bs_init.pd_, bs_flux_connections<base_strategy_fif>::bs_type()); BS_ASSERT (res);
    //res &= BS_KERNEL.register_type(*bs_init.pd_, bs_flux_connections<base_strategy_did>::bs_type()); BS_ASSERT (res);
    //res &= BS_KERNEL.register_type(*bs_init.pd_, bs_flux_connections<base_strategy_dif>::bs_type()); BS_ASSERT (res);

    ////mesh_rs(virtual)
    //res &= blue_sky::give_kernel::Instance().register_type(*bs_init.pd_, mesh_rs<base_strategy_fif>::bs_type()); BS_ASSERT (res);
    //res &= blue_sky::give_kernel::Instance().register_type(*bs_init.pd_, mesh_rs<base_strategy_did>::bs_type()); BS_ASSERT (res);

    //  //mesh_ijk
    //res &= blue_sky::give_kernel::Instance().register_type(*bs_init.pd_, mesh_ijk<base_strategy_fif>::bs_type()); BS_ASSERT (res);
    //res &= blue_sky::give_kernel::Instance().register_type(*bs_init.pd_, mesh_ijk<base_strategy_did>::bs_type()); BS_ASSERT (res);
    
    res &= BS_KERNEL.register_type(pd, mesh_ijk_keywords::bs_type()); BS_ASSERT (res);
  
    res &= BS_KERNEL.register_type(pd, bs_mesh_grdecl::bs_type()); BS_ASSERT (res);
    
    res &= BS_KERNEL.register_type(pd, mesh_grdecl_keywords::bs_type()); BS_ASSERT (res);
    
    //mpi_mesh_grdecl
#ifdef _MPI_MY
    res &= BS_KERNEL.register_type(pd, mpi_mesh_grdecl<base_strategy_fif>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(pd, mpi_mesh_grdecl<base_strategy_did>::bs_type()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type(pd, mpi_mesh_grdecl<base_strategy_dif>::bs_type()); BS_ASSERT (res);
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
/*-----------------------------------------------------------------
 * forward declaration of exporting functions
 *----------------------------------------------------------------*/
namespace blue_sky { namespace python {

void py_export_mesh_grdecl();
void py_export_flux_connections();

}}

/*-----------------------------------------------------------------
 * callback that make Python bindings
 *----------------------------------------------------------------*/
BLUE_SKY_INIT_PY_FUN
{
  using namespace blue_sky;

  python::py_export_mesh ();
  python::py_export_mesh_grdecl ();
  python::py_export_flux_connections ();
}
#ifdef _DEBUG
BOOST_PYTHON_MODULE (bs_mesh_d)
#else
BOOST_PYTHON_MODULE (bs_mesh)
#endif
{
  bs_init_py_subsystem ();
  std::cout << &BS_KERNEL << std::endl;
  bool res = blue_sky::register_types (*blue_sky::bs_get_plugin_descriptor ());
  if (!res)
    throw "Can't register mesh types";
}
#endif //BSPY_EXPORT_PLUGIN

