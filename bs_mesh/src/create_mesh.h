/**
 * \file create_mesh.h
 * \brief helper to create mesh of specified type
 * \author Sergey Miryanov
 * \date 22.06.2009
 * */
#ifndef BS_MESH_CREATE_MESH_ADAPTER_H_
#define BS_MESH_CREATE_MESH_ADAPTER_H_

#include "rs_mesh_iface.h"
#include "rs_smesh_iface.h"
#include "flux_connections_iface.h"

namespace blue_sky {
namespace mesh_detail {

  //class BS_API_PLUGIN create_mesh
  //  {
  //  public:
  //    typedef rs_mesh_iface mesh_iface_t;
  //    typedef rs_smesh_iface smesh_iface_t;

  //    static smart_ptr <mesh_iface_t, true>  
  //    create (const std::string &mesh_name, const smart_ptr <mesh_iface_t, true> &mesh);

  //    static smart_ptr <smesh_iface_t, true>
  //    get_structured (const smart_ptr <mesh_iface_t, true> &mesh);
  //  };
  //
  //template <typename flux_conn_t>
  //class BS_API_PLUGIN create_flux_conn
  //{
  //public:
  //  typedef typename flux_conn_t::strategy_type strategy_t;
  //  typedef flux_connections_iface flux_conn_iface_t;
  //  
  //  static smart_ptr <flux_conn_iface_t, true>  
  //  create ();
  //};

} // namespace mesh_detail
} // namespace blue_sky

#endif // #ifndef BS_MESH_CREATE_MESH_ADAPTER_H_

