/**
 * \file create_mesh.cpp
 * \brief helper to create mesh of specified type 
 * \author Sergey Miryanov
 * \date 22.06.2009
 * */
#include "bs_mesh_stdafx.h"
#include "strategy_name.h"
#include "create_mesh.h"

#include "bs_mesh_grdecl.h"
//#include "bs_mesh_grdecl_mpfa.h"
#include "bs_mesh_ijk.h"
#include "bs_flux_connections.h"

namespace blue_sky {
namespace mesh_detail {

  //smart_ptr <create_mesh::mesh_iface_t, true>
  //create_mesh::create (const std::string &mesh_name, const smart_ptr <create_mesh::mesh_iface_t, true> &mesh)
  //{
  //  const std::vector <type_tuple> &types = BS_KERNEL.registered_types ();
  //  
  //  for (size_t i = 0, cnt = types.size (); i < cnt; ++i)
  //    {
  //      const type_descriptor &td = types[i].td_;
  //      if (td.stype_ == mesh_name)
  //        {
  //          smart_ptr <create_mesh::mesh_iface_t, true> new_mesh (BS_KERNEL.create_object (td), bs_dynamic_cast ());
  //          if (!new_mesh)
  //            {
  //              bs_throw_exception (boost::format ("Type (%s) not registered") % mesh_name);
  //            }

  //          return new_mesh;
  //        }
  //    }

  //  bs_throw_exception ("Unknown mesh typename");
  //    

  //
  //  
  //    
  //  /*
  //  if (!mesh)
  //    {
  //     
  //        
  //      smart_ptr <typename create_mesh <mesh_t>::mesh_iface_t, true> new_mesh = BS_KERNEL.create_object (mesh_t::bs_type ());
  //      if (!new_mesh)
  //        {
  //          bs_throw_exception ("Can't create mesh");
  //        }

  //      return new_mesh;
  //    }
  //  else
  //    {
  //      smart_ptr <typename create_mesh <mesh_t>::mesh_iface_t, true> p_mesh (mesh, bs_dynamic_cast ());
  //      if (p_mesh)
  //        {
  //          return p_mesh;
  //        }
  //      else
  //        {
  //          bs_throw_exception (boost::format ("Can't cast passed mesh (type: %s)") % bs::type_name (mesh));
  //        }
  //    }
  //   */
  //}

  //smart_ptr <create_mesh::smesh_iface_t, true>
  //create_mesh::get_structured (const smart_ptr <mesh_iface_t, true> &mesh)
  //{
  //  smart_ptr <smesh_iface_t, true> structured_mesh (mesh, bs_dynamic_cast ());
  //  if (!structured_mesh)
  //    {
  //      bs_throw_exception ("Passed mesh is not structured");
  //    }

  //  return structured_mesh;
  //}

  //smart_ptr <create_flux_conn::flux_conn_iface_t, true>
  //create_flux_conn <flux_conn_t>::create ()
  //{
  //  smart_ptr <flux_conn_t, true> new_flux_conn (BS_KERNEL.create_object (flux_conn_t::bs_type ()), bs_dynamic_cast ());
  //  if (!new_flux_conn)
  //    {
  //      bs_throw_exception ("Can't create flux connections!");
  //    }

  //  return new_flux_conn;
  //}

} // namespace mesh_detail
} // namespace blue_sky


