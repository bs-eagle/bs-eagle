/*!
	\file rs_mesh_base.cpp
  \brief This file implements base class for reservoir simulation meshes
  \author Mark Khait
  \date 2009-07-20
  */
  
#include "bs_mesh_stdafx.h"

#ifndef PURE_MESH
  using namespace blue_sky;
#endif

#include "rs_mesh_base.h"


rs_mesh_base ::rs_mesh_base ()
#ifndef PURE_MESH
: depths (give_kernel::Instance().create_object(v_float::bs_type ()))
#endif
{
}


void
rs_mesh_base ::init_props (const sp_hdm_t hdm)
{
#ifndef PURE_MESH  
  well_pool = hdm->get_well_pool();
  minpv = hdm->get_prop ()->get_f (L"minimal_pore_volume");
  minsv = hdm->get_prop ()->get_f (L"minimal_splice_volume");
  max_thickness = hdm->get_prop ()->get_f (L"maximum_splice_thickness");
  darcy_constant = hdm->get_darcy_constant ();
  
  actnum_array = hdm->get_pool ()->get_i_data("ACTNUM");
  poro_array = hdm->get_pool ()->get_fp_data("PORO");
  ntg_array = hdm->get_pool ()->get_fp_data("NTG");
  multpv_array = hdm->get_pool ()->get_fp_data("MULTPV");
  
  n_elements = static_cast <t_long> (actnum_array->size ());
  n_active_elements =  std::accumulate (actnum_array->begin (), actnum_array->end (),0);
#else
  minpv = hdm->minpv;
  minsv = hdm->minsv;
  max_thickness = hdm->max_thickness;
  darcy_constant = hdm->darcy_constant;
  
  actnum_array = hdm->actnum_array;
  poro_array = hdm->poro_array;
  ntg_array = hdm->ntg_array;
  multpv_array = hdm->multpv_array;
  
  n_elements = hdm->n_elements;
  //n_active_elements =  std::accumulate (actnum_array->begin (), actnum_array->end (),0);
#endif

}

int rs_mesh_base::init_int_to_ext()
{
#ifndef PURE_MESH
  t_long n_ext = t_long(ext_to_int->size());
  if (n_ext == 0)
    return -1;

  //int_to_ext->init (n_active_elements, 0);
  int_to_ext->resize (n_active_elements);
  
  t_long *int_to_ext_data = int_to_ext->data ();
  t_long *ext_to_int_data = ext_to_int->data ();
  
  stdv_long wrong;
  stdv_long wrong_idx;

#else
  t_long n_ext = n_elements;
  if (!ext_to_int)
    return -1;

  spv_long int_to_ext_data = int_to_ext;
  spv_long ext_to_int_data = ext_to_int;
#endif
    
  
  for (t_long i = 0; i < n_ext; i++)
    {
      if (ext_to_int_data[i] != -1)
        {
#ifndef PURE_MESH          
          if (ext_to_int_data[i] >= n_active_elements)
            {
              wrong.push_back (ext_to_int_data[i]);
              wrong_idx.push_back (i);
            }
          else
#endif
            {
              int_to_ext_data[ext_to_int_data[i]] = i;
            }
        }
    }

#ifndef PURE_MESH
  if (wrong.size ())
    {
      for (size_t i = 0, cnt = wrong.size (); i < cnt; ++i)
        {
          BOSERR (section::mesh, level::error) 
            << boost::format ("ext_to_int[%d] == %s >= %d") % wrong_idx[i] % wrong[i] % n_active_elements 
            << bs_end;
        }

      bs_throw_exception ("ext_to_int out of n_active_elements");
    }
#endif
  return 0;  
}


void rs_mesh_base::check_data() const
{
  base_t::check_data ();
#ifndef PURE_MESH 
  if (minpv < 0)
    bs_throw_exception (boost::format ("minpv = %d is out of range")% minpv);
  if (minsv < 0)
    bs_throw_exception (boost::format ("minsv = %d is out of range")% minsv);
  if (max_thickness < 0)
    bs_throw_exception (boost::format ("max_thickness = %d is out of range")% max_thickness);
    
  // FIXME: get_fp_data, get_i_data in init_props will raise exception if no such array
  // maybe we can remove this checks, or don't raise exceptions in init_props
  if (!actnum_array)
    bs_throw_exception ("ACTNUM array is not initialized");
  if (!poro_array)
    bs_throw_exception ("PORO array is not initialized");
  if (!depths->size ())
    bs_throw_exception ("depths array is not initialized");
#endif
}


//BS_INST_STRAT(rs_mesh_base);
