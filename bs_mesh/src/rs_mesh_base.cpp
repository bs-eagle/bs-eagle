/*!
	\file rs_mesh_base.cpp
  \brief This file implements base class for reservoir simulation meshes
  \author Mark Khait
  \date 2009-07-20
  */
  
#include "bs_mesh_stdafx.h"
#include "rs_mesh_base.h"

using namespace blue_sky;

void
rs_mesh_base::init_props (const sp_idata_t &idata)
{
  //base_t::init_props (idata);
  
  minpv = idata->minimal_pore_volume;
  minsv = idata->minimal_splice_volume;
  max_thickness = idata->maximum_splice_thickness;

  sp_actnum = idata->get_int_non_empty_array("ACTNUM");
  sp_poro = idata->get_float_non_empty_array("PORO");
  sp_ntg = idata->get_float_array("NTG");
  sp_multpv = idata->get_float_array("MULTPV");
  
  base_t::n_elements = static_cast <index_t> (sp_actnum.size ());
  base_t::n_active_elements =  std::accumulate(sp_actnum.begin(), sp_actnum.end(),0);
}


int rs_mesh_base::init_int_to_ext()
{
  if (base_t::ext_to_int.size() == 0)
    return -1;
    
  base_t::int_to_ext.resize (base_t::n_active_elements, 0);

  for (size_t i = 0; i < base_t::ext_to_int.size(); i++)
    {
      if (base_t::ext_to_int[i] != -1)
        {
          base_t::int_to_ext [base_t::ext_to_int[i]] = (index_t)i;
        }
    }
  return 0;  
}

void rs_mesh_base::check_data() const
{
  base_t::check_data ();
  
  if (minpv < 0)
    bs_throw_exception (boost::format ("minpv = %d is out of range")% minpv);
  if (minsv < 0)
    bs_throw_exception (boost::format ("minsv = %d is out of range")% minsv);
  if (max_thickness < 0)
    bs_throw_exception (boost::format ("max_thickness = %d is out of range")% max_thickness);
    
  if (!sp_actnum.size ())
    bs_throw_exception ("ACTNUM array is not initialized");
  if (!sp_poro.size ())
    bs_throw_exception ("PORO array is not initialized");
  if (!depths.size ())
    bs_throw_exception ("depths array is not initialized");
}

hdf5_group_v2 &
rs_mesh_base::save_info (hdf5_group_v2 &group) const
{
  // FIXME: dims
  //  smart_ptr <rs_smesh_iface> smesh (mesh, bs_dynamic_cast ());
  //  if (!smesh)
  //    {
  //      bs_throw_exception (boost::format ("Can't cast mesh to structired mesh, %s") % mesh->bs_resolve_type ().stype_);
  //    }
  //  const rs_smesh_iface::index_point3d_t &dims = smesh->get_dimens ();

  //  hdf5_group_v2 group = file_["/mesh"];

  //  group.write ("initial_data", 
  //    hdf5_pod (dims[0]) << dims[1] << dims[2]
  //      << smesh->get_n_active_elements () << smesh->get_n_active_elements ());

  //  group.write ("original_elements", smesh->get_ext_to_int ());
  //  group.write ("original_planes",   smesh->get_ext_to_int ());
  //
  group.write ("initial_data", hdf5_pod (0) << 0 << 0 << this->get_n_active_elements () << this->get_n_active_elements ());
  group.write ("original_elements", this->get_ext_to_int ());
  group.write ("original_planes", this->get_ext_to_int ());
  group.write ("params", hdf5_pod (this->minpv) << this->minsv << this->max_thickness);

  return group;
}

hdf5_group_v2 &
rs_mesh_base::save_data (hdf5_group_v2 &group) const
{
  group
    .write ("actnum",   this->sp_actnum)
    .write ("poro",     this->sp_poro)
    .write ("ntg",      this->sp_ntg)
    .write ("multpv",   this->sp_multpv)
    .write ("depths",   this->depths)
    .write ("volumes",  this->volumes)
    .write ("ext_to_int", this->get_ext_to_int ())
    .write ("int_to_ext", this->get_int_to_ext ())
    ;

  return group;
}
