/// @file scal_2p_data_holder.cpp
/// @brief scal_2p_data_holder implementation
/// @author uentity
/// @version 
/// @date 26.01.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "scal_3p.h"
#include "scale_array_holder.h"
#include "scal_region_info.h"
#include "scal_region.h"

#include "scal_data_source.h"
#include "scal_data_vector.h"
#include "scale_arrays_placement_strategies.h"
#include "scal_data_placement_strategies.h"

#include "scal_2p_data_holder.h"

namespace blue_sky {
  //////////////////////////////////////////////////////////////////////////

  scal_2p_data_holder::scal_2p_data_holder (bs_type_ctor_param param /* = NULL */)
  {
    data_ = BS_KERNEL.create_object (item_array_t::bs_type ());
    scal_table_array.clear ();
  }

  scal_2p_data_holder::scal_2p_data_holder (const this_t& s)
  : bs_refcounter (s), scal_2p_data_holder_iface (s)
  {
    BS_ASSERT (false && "NOT IMPL YET");
  }

  void
  scal_2p_data_holder::precalc (scal_region_info_t &info, const scal_region_t &region, scal::data_placement::scal_data_placement_type type, bool is_water)
  {
    // replace calc_sorp and calc_spr
    precalc_min_index_and_spr (region.Krp, region.Sp, info.Krp_min_greater_zero, info.spr);
    precalc_min_index_and_spr (region.Krop, region.So, info.Krop_min_greater_zero, info.sorp);

#ifdef _DEBUG
    BOSOUT (section::scal, level::debug) << boost::format ("krp_min: {%d} krop: {%d}") % info.Krp_min_greater_zero % info.Krop_min_greater_zero << bs_end;
#endif

    if (is_water)  // water-oil system
      {
        info.pcp_max = fabs (region.Pcp.front ());
      }
    else // gas-oil system
      {
        info.pcp_max = fabs (region.Pcp.back ());
      }
       
  }

  void
  scal_2p_data_holder::add_spof (sp_array_item_t const &data, t_long region_index, bool is_water)
  {
    BS_ASSERT ((data->size () % 4) == 0) (data->size ());

    scal_region_info_t info;
    info.So_count   = (int)data->size () / 4;
    info.Sp_count   = info.So_count;
    info.sp_offset  = (int)data_->size ();
    info.so_offset  = (int)data_->size ();

    data_->resize (data_->size () + info.So_count * 5);
    placement_info_.type = scal::data_placement::spof;

#ifdef _DEBUG
    item_t *raw_data = &(*data_)[0];
    raw_data;
#endif

    scal::data_placement::all_regions_t::place_spof_data (data_, placement_info_, data, is_water);
    region_.push_back (info);
    scal_region_t const &new_region = get_region_from_info (region_.size () - 1);
    precalc (region_.back (), new_region, scal::data_placement::spof, is_water);
    if (!check (new_region, is_water))
      {
        //region_.pop_back ();
        //data_.resize (data_.size () - info.So_count * 5);
        throw bs_exception ("scal_2p_data_holder::add_spof", "Could not compute residual saturation");
        // TODO: LOG
      }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    if (scal_table_array.size ())
      {
        BS_ASSERT (region_index >= 0 && region_index < scal_table_array.size ());
        t_int table_len = (t_int) data->size () / 4;
        t_float const *data_array = data->data ();
        
        spv_double sp = BS_KERNEL.create_object (v_double::bs_type ()); 
        spv_double so = BS_KERNEL.create_object (v_double::bs_type ()); 
        spv_double krp = BS_KERNEL.create_object (v_double::bs_type ()); 
        spv_double krop = BS_KERNEL.create_object (v_double::bs_type ()); 
        spv_double pcp = BS_KERNEL.create_object (v_double::bs_type ()); 
        sp->resize (table_len);
        so->resize (table_len);
        krp->resize (table_len);
        krop->resize (table_len);
        pcp->resize (table_len);
        
        t_double *sp_ = &(*sp)[0];
        t_double *so_ = &(*so)[0];
        t_double *krp_ = &(*krp)[0];
        t_double *krop_ = &(*krop)[0];
        t_double *pcp_ = &(*pcp)[0];
        
        for (t_long i = 0; i < table_len; i++)
          {
            sp_[i]   = data_array[i*4 + 0];
            so_[i]   = 1.0 - sp_[i];
            krp_[i]  = data_array[i*4 + 1];
            krop_[i] = data_array[i*4 + 2];
            pcp_[i]  = is_water ? -data_array[i*4 + 3] : data_array[i*4 + 3];
          }
        
        sp_table scal_table = scal_table_array[region_index]; 
        if (is_water)
          {
            scal_table->add_col_vector (SCAL_TABLE_SP, L"SW", sp);
            scal_table->add_col_vector (SCAL_TABLE_SO, L"SO", so);
            scal_table->add_col_vector (SCAL_TABLE_KRP, L"KRW", krp);
            scal_table->add_col_vector (SCAL_TABLE_KROP, L"KROW", krop);  
            scal_table->add_col_vector (SCAL_TABLE_PCP, L"PCW", pcp);
          }
        else
          {
            scal_table->add_col_vector (SCAL_TABLE_SP, L"SG", sp);
            scal_table->add_col_vector (SCAL_TABLE_SO, L"SO", so);
            scal_table->add_col_vector (SCAL_TABLE_KRP, L"KRG", krp);
            scal_table->add_col_vector (SCAL_TABLE_KROP, L"KROG", krop);  
            scal_table->add_col_vector (SCAL_TABLE_PCP, L"PCG", pcp);
          }      
      }
    ///////////////////////////////////////////////////////////////////////////////    
  }

  void
  scal_2p_data_holder::add_spfn (sp_array_item_t const &data, t_long region_index, bool is_water)
  {
    BS_ASSERT ((data->size () % 3) == 0) (data->size ());

    if (placement_info_.type == scal::data_placement::scal_data_placement_null ||
        placement_info_.type == scal::data_placement::spfn_sofX)
      {
        scal_region_info_t info;
        info.So_count   = -1;
        info.Sp_count   = (int)data->size () / 3;
        info.so_offset  = -1;
        info.sp_offset  = (int)data_->size ();

        data_->resize (data_->size () + data->size ());
        placement_info_.type = scal::data_placement::spfn_sofX;

        scal::data_placement::all_regions_t::place_spfn_data (data_, placement_info_, data, is_water);
        region_.push_back (info);
      }
    else if (placement_info_.type == scal::data_placement::sofX_spfn)
      {
        BS_ASSERT (region_index < static_cast <t_long> (region_.size ())) (region_index) (region_.size ());
        scal_region_info_t &info = region_[region_index];

        info.Sp_count   = (int)data->size () / 3;
        info.sp_offset  = (int)data_->size ();

        data_->resize (data_->size () + data->size ());
        scal::data_placement::all_regions_t::place_spfn_data (data_, placement_info_, data, is_water);
        scal_region_t const &new_region = get_region_from_info (region_index);
        precalc (info, new_region, scal::data_placement::sofX_spfn, is_water);
        if (!check (new_region, is_water))
          {
            //region_.pop_back ();
            //data_.resize (data_.size () - info.So_count * 5);
            throw bs_exception ("scal_2p_data_holder::add_spfn", "Could not compute residual saturation");
            // TODO: LOG
          }
      }
    else
      {
        BS_ASSERT (false && "UNSUPPORTED DATA PLACEMENT TYPE") (placement_info_.type);
      }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////
    if (scal_table_array.size ())
      {
        BS_ASSERT (region_index >= 0 && region_index < scal_table_array.size ());
        t_int table_len = (t_int) data->size () / 3;
        t_float const *data_array = data->data ();
        
        spv_double sp = BS_KERNEL.create_object (v_double::bs_type ()); 
        spv_double krp = BS_KERNEL.create_object (v_double::bs_type ()); 
        spv_double pcp = BS_KERNEL.create_object (v_double::bs_type ()); 
        sp->resize (table_len);
        krp->resize (table_len);
        pcp->resize (table_len);
        
        t_double *sp_ = &(*sp)[0];
        t_double *krp_ = &(*krp)[0];
        t_double *pcp_ = &(*pcp)[0];
        
        for (t_long i = 0; i < table_len; i++)
          {
            sp_[i]   = data_array[i*3 + 0];
            krp_[i]  = data_array[i*3 + 1];
            pcp_[i]  = is_water ? -data_array[i*3 + 2] : data_array[i*3 + 2];
          }
        
        sp_table scal_table = scal_table_array[region_index]; 
        if (is_water)
          {
            scal_table->add_col_vector (SCAL_TABLE_SP, L"SW", sp);
            scal_table->add_col_vector (SCAL_TABLE_KRP, L"KRW", krp);
            scal_table->add_col_vector (SCAL_TABLE_PCP, L"PCW", pcp);
          }
        else
          {
            scal_table->add_col_vector (SCAL_TABLE_SP, L"SG", sp);
            scal_table->add_col_vector (SCAL_TABLE_KRP, L"KRG", krp);
            scal_table->add_col_vector (SCAL_TABLE_PCP, L"PCG", pcp);
          }      
      }
    ///////////////////////////////////////////////////////////////////////////////    
      
  }

  void
  scal_2p_data_holder::add_sof3 (sp_array_item_t const &data, t_long region_index, bool is_water)
  {
    typedef t_int   index_t;

    BS_ASSERT ((data->size () % 3) == 0) (data->size ());

    if (placement_info_.type == scal::data_placement::scal_data_placement_null ||
        placement_info_.type == scal::data_placement::sofX_spfn)
      {
        scal_region_info_t info;
        info.So_count   = (int)data->size () / 3;
        info.Sp_count   = -1;
        info.so_offset  = (int)data_->size ();
        info.sp_offset  = -1;

        data_->resize (data_->size () + info.So_count * 2);
        placement_info_.type = scal::data_placement::sofX_spfn;

        scal::data_placement::all_regions_t::place_sof3_data (data_, placement_info_, data, is_water);
        region_.push_back (info);
      }
    else if (placement_info_.type == scal::data_placement::spfn_sofX)
      {
        BS_ASSERT (region_index < static_cast <t_long> (region_.size ())) (region_index) (region_.size ());
        scal_region_info_t &info = region_[region_index];

        info.So_count   = (int)data->size () / 3;
        info.so_offset  = (int)data_->size ();

        data_->resize (data_->size () + info.So_count * 2);
        scal::data_placement::all_regions_t::place_sof3_data (data_, placement_info_, data, is_water);
        scal_region_t const &new_region = get_region_from_info (region_index);
        precalc (info, new_region, scal::data_placement::spfn_sofX, is_water);
        if (!check (new_region, is_water))
          {
            //region_.pop_back ();
            //data_.resize (data_.size () - info.So_count * 5);
            throw bs_exception ("scal_2p_data_holder::add_sof3", "Could not compute residual saturation");
            // TODO: LOG
          }
      }
    else
      {
        BS_ASSERT (false && "UNSUPPORTED DATA PLACEMENT TYPE") (placement_info_.type);
      }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    if (scal_table_array.size ())
      {
        BS_ASSERT (region_index >= 0 && region_index < scal_table_array.size ());
        t_int table_len = (t_int) data->size () / 3;
        t_float const *data_array = data->data ();
        
        spv_double so = BS_KERNEL.create_object (v_double::bs_type ()); 
        spv_double krop = BS_KERNEL.create_object (v_double::bs_type ()); 
        so->resize (table_len);
        krop->resize (table_len);
        
        t_double *so_ = &(*so)[0];
        t_double *krop_ = &(*krop)[0];
        
        for (t_long i = 0; i < table_len; i++)
          {
            so_[i]   = data_array[i*3 + 0];
            krop_[i] = is_water ? data_array[i*3 + 1] : data_array[i*3 + 2];
          }
        
        sp_table scal_table = scal_table_array[region_index]; 
        if (is_water)
          {
            scal_table->add_col_vector (SCAL_TABLE_SO, L"SO", so);
            scal_table->add_col_vector (SCAL_TABLE_KROP, L"KROW", krop);  
          }
        else
          {
            scal_table->add_col_vector (SCAL_TABLE_SO, L"SO", so);
            scal_table->add_col_vector (SCAL_TABLE_KROP, L"KROG", krop);  
          }      
      }
    ///////////////////////////////////////////////////////////////////////////////    
      
  }

  void
  scal_2p_data_holder::add_sof2 (sp_array_item_t const &data, t_long region_index, bool is_water)
  {
    typedef t_int   index_t;

    BS_ASSERT ((data->size () % 2) == 0) (data->size ());

    if (placement_info_.type == scal::data_placement::scal_data_placement_null ||
        placement_info_.type == scal::data_placement::sofX_spfn)
      {
        scal_region_info_t info;
        info.So_count   = (int)data->size () / 2;
        info.Sp_count   = -1;
        info.so_offset  = (int)data_->size ();
        info.sp_offset  = -1;

        data_->resize (data_->size () + data->size ());
        placement_info_.type = scal::data_placement::sofX_spfn;

        scal::data_placement::all_regions_t::place_sof2_data (data_, placement_info_, data);
        region_.push_back (info);
      }
    else if (placement_info_.type == scal::data_placement::spfn_sofX)
      {
        BS_ASSERT (region_index < region_.size ()) (region_index) (region_.size ());
        scal_region_info_t &info = region_[region_index];

        info.So_count   = (int)data->size () / 2;
        info.so_offset  = (int)data_->size ();

        data_->resize (data_->size () + data->size ());
        scal::data_placement::all_regions_t::place_sof2_data (data_, placement_info_, data);
        scal_region_t const &new_region = get_region_from_info (region_index);
        precalc (info, new_region, scal::data_placement::spfn_sofX, is_water);
        if (!check (new_region, is_water))
          {
            //region_.pop_back ();
            //data_.resize (data_.size () - info.So_count * 5);
            throw bs_exception ("scal_2p_data_holder::add_sof2", "Could not compute residual saturation");
            // TODO: LOG
          }
      }
    else
      {
        BS_ASSERT (false && "UNSUPPORTED DATA PLACEMENT TYPE") (placement_info_.type);
      }
    /////////////////////////////////////////////////////////////////////////////////////////////////
    if (scal_table_array.size ())
      {
        BS_ASSERT (region_index >= 0 && region_index < scal_table_array.size ());
        t_int table_len = (t_int) data->size () / 2;
        t_float const *data_array = data->data ();
        
        spv_double so = BS_KERNEL.create_object (v_double::bs_type ()); 
        spv_double krop = BS_KERNEL.create_object (v_double::bs_type ()); 
        so->resize (table_len);
        krop->resize (table_len);
        
        t_double *so_ = &(*so)[0];
        t_double *krop_ = &(*krop)[0];
        
        for (t_long i = 0; i < table_len; i++)
          {
            so_[i]   = data_array[i*2 + 0];
            krop_[i] = data_array[i*2 + 1];
          }
        
        sp_table scal_table = scal_table_array[region_index]; 
        if (is_water)
          {
            scal_table->add_col_vector (SCAL_TABLE_SO, L"SO", so);
            scal_table->add_col_vector (SCAL_TABLE_KROP, L"KROW", krop);  
          }
        else
          {
            scal_table->add_col_vector (SCAL_TABLE_SO, L"SO", so);
            scal_table->add_col_vector (SCAL_TABLE_KROP, L"KROG", krop);  
          }      
      }
    ///////////////////////////////////////////////////////////////////////////////    
      
  }

  bool
  scal_2p_data_holder::check (const scal_region_t &region, bool is_water)
  {
    if (is_water)
      {
        for (size_t i = 0, cnt = region.Krp.size (); i < cnt; ++i)
          {
            if (fabs (region.Krp[(int)i]) > EPS_DIFF)
              {
#if 0
                if (i == region.Krp.size() - 1)
                  return false;
                else
#endif 
                  return true;
              }
          }
      }
    else
      {
        size_t k = region.Krp.size () - 1;
        while (region.Krp[(int)k] == region.Krp[(int)k - 1] && k > 0)
          --k;

        if (k == 0)
          return false;
      }

    return true;
  }

  void
  scal_2p_data_holder::precalc_min_index_and_spr (const data_vector_t &main, const data_vector_t &slave, int &min_index, item_t &spr)
  {
    for (int i = 0, found = 0, cnt = (int)main.size (); i < cnt; ++i)
      {
        // store index of Krop that we will use in interpolation
        if (!found && main[i] > 0.0)
          {
            min_index = i;
            found = true;
          }

        // calc sorp
        if (main[i] < EPS_DIFF)  // BUG: for spr (not for sorp) this condition is main[i] == 0.0
          {
            spr = (slave[i] > 0.0 ? slave[i] : 0.0);
          }
      }
  }

  void
  scal_2p_data_holder::update_so (const sp_scal_data_t &water_data)
  {
    BS_ERROR (water_data->get_region_info ().size () == region_.size (), "update_so");// (water_data->get_region_info ().size ()) (region_.size ());

    for (size_t i = 0, cnt = region_.size (); i < cnt; ++i)
      {
        scal_region_t gas_region = get_region ((int)i);
        const scal_region_t &water_region = water_data->get_region ((int)i);

        item_t swc = water_region.Sp[0];
        for (size_t j = 0, jcnt = gas_region.So.size (); j < jcnt; ++j)
          {
            gas_region.So[(int)j] -= swc; // TODO: BUG:
          }
      }
    ////////////////////////////////////////////////////////////////////
    if (scal_table_array.size ())
      {
        BS_ASSERT (scal_table_array.size () == water_data->scal_table_array.size ());
        t_long n_regions = scal_table_array.size ();
        for (t_long i = 0; i < n_regions; i++)
          {
            sp_table water_scal_table = water_data->scal_table_array[i];
            item_t swc = water_scal_table->get_value (0, SCAL_TABLE_SP);
            
            sp_table gas_scal_table = scal_table_array[i];
            t_int so_table_len = gas_scal_table->get_n_rows (SCAL_TABLE_SO);
            BS_ASSERT (so_table_len > 0);
            t_double *so_ = gas_scal_table->get_col_ptr (SCAL_TABLE_SO);
            for (t_int j = 0; j < so_table_len; j++)
              {
                so_[i] -= swc;
              }
          }
      }
    ////////////////////////////////////////////////////////////////////
      
  }

  BLUE_SKY_TYPE_STD_CREATE (scal_2p_data_holder);
  BLUE_SKY_TYPE_STD_COPY (scal_2p_data_holder);
  BLUE_SKY_TYPE_IMPL(scal_2p_data_holder, scal_2p_data_holder_iface, "scal_2p_data_holder", "scal_2p_data_holder calculation class", "scal_2p_data_holder calculation");

} // eof blue_sky
