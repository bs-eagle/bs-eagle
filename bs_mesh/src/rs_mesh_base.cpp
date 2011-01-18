/*!
	\file rs_mesh_base.cpp
  \brief This file implements base class for reservoir simulation meshes
  \author Mark Khait
  \date 2009-07-20
  */
  
#include "bs_mesh_stdafx.h"
#include "rs_mesh_base.h"

using namespace blue_sky;

template <typename strategy_t>
rs_mesh_base <strategy_t>::rs_mesh_base ()
{
  depths = give_kernel::Instance().create_object(bs_array<fp_type_t>::bs_type());
  actnum_array = 0; 
  poro_array = 0;
  ntg_array = 0;
  multpv_array = 0;
  darcy_constant = ph_const.darcy_constant;
}

template <typename strategy_t>
void
rs_mesh_base <strategy_t>::init_props (const sp_idata_t &idata)
{
  //base_t::init_props (idata);
  
  minpv = idata->minimal_pore_volume;
  minsv = idata->minimal_splice_volume;
  max_thickness = idata->maximum_splice_thickness;
  
  sp_fp_storage_array_t data_array;
  sp_i_array_t sp_actnum_array;
  
  sp_actnum_array = idata->get_int_non_empty_array("ACTNUM");
  if (sp_actnum_array->size()) actnum_array = &(*sp_actnum_array)[0];
  
  base_t::n_elements = static_cast <i_type_t> (sp_actnum_array->size ());
  base_t::n_active_elements =  std::accumulate(sp_actnum_array->begin(), sp_actnum_array->end(),0);
  
  data_array = idata->get_fp_non_empty_array("PORO");
  if (data_array->size()) poro_array = &(*data_array)[0];
  
  data_array = idata->get_fp_array("NTG");
  if (data_array && data_array->size()) ntg_array = &(*data_array)[0];
  
  data_array = idata->get_fp_array("MULTPV");
  if (data_array && data_array->size()) multpv_array = &(*data_array)[0];
}


template<class strategy_t>
int rs_mesh_base<strategy_t>::init_int_to_ext()
{
  i_type_t *ext_to_int_data, *int_to_ext_data; 
  i_type_t n_ext = base_t::ext_to_int->size();
  if (n_ext == 0)
    return -1;
    
  base_t::int_to_ext->resize (base_t::n_active_elements);
  base_t::int_to_ext->assign(0);
  
  int_to_ext_data = &(*base_t::int_to_ext)[0];
  ext_to_int_data = &(*base_t::ext_to_int)[0];
  
  for (i_type_t i = 0; i < n_ext; i++)
    {
      if (ext_to_int_data[i] != -1)
        {
          int_to_ext_data[ext_to_int_data[i]] = i;
        }
    }
  return 0;  
}

template<class strategy_t>
void rs_mesh_base<strategy_t>::check_data() const
{
  base_t::check_data ();
  
  if (minpv < 0)
    bs_throw_exception (boost::format ("minpv = %d is out of range")% minpv);
  if (minsv < 0)
    bs_throw_exception (boost::format ("minsv = %d is out of range")% minsv);
  if (max_thickness < 0)
    bs_throw_exception (boost::format ("max_thickness = %d is out of range")% max_thickness);
    
  if (!actnum_array)
    bs_throw_exception ("ACTNUM array is not initialized");
  if (!poro_array)
    bs_throw_exception ("PORO array is not initialized");
  if (!depths->size ())
    bs_throw_exception ("depths array is not initialized");
}


BS_INST_STRAT(rs_mesh_base);
