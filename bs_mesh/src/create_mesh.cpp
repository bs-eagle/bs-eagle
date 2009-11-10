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

  template <typename strategy_t>
  smart_ptr <typename create_mesh <strategy_t>::mesh_iface_t, true>
  create_mesh <strategy_t>::create (const std::string &mesh_name, const smart_ptr <typename create_mesh <strategy_t>::mesh_iface_t, true> &mesh)
  {
    const std::vector <type_tuple> &types = BS_KERNEL.registered_types ();
    
    std::string mesh_strat_name = mesh_name + "_" + tools::strategy_name <strategy_t>::name ();
    for (size_t i = 0, cnt = types.size (); i < cnt; ++i)
      {
        const type_descriptor &td = types[i].td_;
        if (td.stype_ == mesh_strat_name)
          {
            smart_ptr <typename create_mesh <strategy_t>::mesh_iface_t, true> new_mesh (BS_KERNEL.create_object (td), bs_dynamic_cast ());
            if (!new_mesh)
              {
                bs_throw_exception (boost::format ("Type (%s) not registered") % mesh_strat_name);
              }

            return new_mesh;
          }
      }

    bs_throw_exception ("Unknown mesh typename");
      

  
    
      
    /*
    if (!mesh)
      {
       
          
        smart_ptr <typename create_mesh <mesh_t>::mesh_iface_t, true> new_mesh = BS_KERNEL.create_object (mesh_t::bs_type ());
        if (!new_mesh)
          {
            bs_throw_exception ("Can't create mesh");
          }

        return new_mesh;
      }
    else
      {
        smart_ptr <typename create_mesh <mesh_t>::mesh_iface_t, true> p_mesh (mesh, bs_dynamic_cast ());
        if (p_mesh)
          {
            return p_mesh;
          }
        else
          {
            bs_throw_exception (boost::format ("Can't cast passed mesh (type: %s)") % bs::type_name (mesh));
          }
      }
     */
  }

  template <typename strategy_t>
  smart_ptr <typename create_mesh <strategy_t>::smesh_iface_t, true>
  create_mesh <strategy_t>::get_structured (const smart_ptr <mesh_iface_t, true> &mesh)
  {
    smart_ptr <smesh_iface_t, true> structured_mesh (mesh, bs_dynamic_cast ());
    if (!structured_mesh)
      {
        bs_throw_exception ("Passed mesh is not structured");
      }

    return structured_mesh;
  }

  template class create_mesh <base_strategy_fi >;
  template class create_mesh <base_strategy_di >;
  template class create_mesh <base_strategy_mixi>; 
  
  
  template <typename flux_conn_t>
  smart_ptr <typename create_flux_conn <flux_conn_t>::flux_conn_iface_t, true>
  create_flux_conn <flux_conn_t>::create ()
  {
    smart_ptr <flux_conn_t, true> new_flux_conn (BS_KERNEL.create_object (flux_conn_t::bs_type ()), bs_dynamic_cast ());
    if (!new_flux_conn)
      {
        bs_throw_exception ("Can't create flux connections!");
      }

    return new_flux_conn;
  }

#define SPEC_FLUX_CONN(T)                                                                 \
  template class create_flux_conn <T <base_strategy_fi> >;  \
  template class create_flux_conn <T <base_strategy_di> >;  \
  template class create_flux_conn <T <base_strategy_mixi> >; 

  SPEC_FLUX_CONN (bs_flux_connections);

} // namespace mesh_detail
} // namespace blue_sky


