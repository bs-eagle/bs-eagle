/**
 * \file scal_3p.cpp
 * \brief
 * \author Sergey Miryanov
 * \date 19.05.2008
 * */
#include "bs_scal_stdafx.h"

#include "scal_3p.h"
#include "scal_data_source.h"
#include "scal_data_vector.h"
#include "scale_arrays_placement_strategies.h"
#include "scal_data_placement_strategies.h"

#include "scale_array_holder.h"
#include "scal_region_info.h"
#include "scal_region.h"
#include "scal_2p_data_holder.h"

#include "jfunction.h"

namespace blue_sky
  {
    BS_TYPE_IMPL_T_EXT_MEM(bs_array, 2, (calc_model_data, bs_array_shared));
    BS_TYPE_IMPL_T_EXT_MEM(bs_array, 2, (calc_model_data, vector_traits));

  //////////////////////////////////////////////////////////////////////////

  scale_array_holder::scale_array_holder (bs_type_ctor_param param /* = NULL */)
  {
    data_pool = BS_KERNEL.create_object ("float_var_table");
    BS_ASSERT (data_pool);

    data_pool->init ( scale_array_name_total);
    
    data_pool->set_col_name (socr, "SOCR");
    data_pool->set_col_name (scr, "SCR");
    data_pool->set_col_name (su, "SU");
    data_pool->set_col_name (sl, "SL");
    data_pool->set_col_name (pcp, "PCP");
    data_pool->set_col_name (krp, "KRP");
    data_pool->set_col_name (krop, "KROP");
    data_pool->set_col_name (krpr, "KRPR");
    data_pool->set_col_name (krorp, "KRORP");
    BOOST_STATIC_ASSERT (krorp == scale_array_name_total - 1);
  }

  scale_array_holder::scale_array_holder (const scale_array_holder& s)
  : bs_refcounter (s), scale_array_holder_iface (s)
  {
    BS_ASSERT (false && "NOT IMPL YET");
  }
  //////////////////////////////////////////////////////////////////////////

  scal_2p_data_holder::scal_2p_data_holder (bs_type_ctor_param param /* = NULL */)
  {
    data_ = BS_KERNEL.create_object (item_array_t::bs_type ());
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
  scal_2p_data_holder::add_spof (sp_array_item_t const &data, bool is_water)
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
    precalc (region_.back (), get_region_internal ((int)region_.size () - 1), scal::data_placement::spof, is_water);
    if (!check (get_region_internal ((int)region_.size () - 1), is_water))
      {
        //region_.pop_back ();
        //data_.resize (data_.size () - info.So_count * 5);
        //throw bs_exception ("scal_2p_data_holder::add_spof", "Could not compute residual saturation");
        // TODO: LOG
      }
  }

  void
  scal_2p_data_holder::add_spfn (sp_array_item_t const &data, t_long region_index, bool is_water)
  {
    BS_ASSERT ((data->size () % 3) == 0) (data->size ());

    if (placement_info_.type == scal::data_placement::scal_data_placement_null ||
        placement_info_.type == scal::data_placement::spfn_sof3)
      {
        scal_region_info_t info;
        info.So_count   = -1;
        info.Sp_count   = (int)data->size () / 3;
        info.so_offset  = -1;
        info.sp_offset  = (int)data_->size ();

        data_->resize (data_->size () + data->size ());
        placement_info_.type = scal::data_placement::spfn_sof3;

        scal::data_placement::all_regions_t::place_spfn_data (data_, placement_info_, data, is_water);
        region_.push_back (info);
      }
    else if (placement_info_.type == scal::data_placement::sof3_spfn)
      {
        BS_ASSERT (region_index < region_.size ()) (region_index) (region_.size ());
        scal_region_info_t &info = region_[region_index];

        info.Sp_count   = (int)data->size () / 3;
        info.sp_offset  = (int)data_->size ();

        data_->resize (data_->size () + data->size ());
        scal::data_placement::all_regions_t::place_spfn_data (data_, placement_info_, data, is_water);
        precalc (info, get_region_internal ((int)region_index), scal::data_placement::sof3_spfn, is_water);
        if (!check (get_region_internal ((int)region_index), is_water))
          {
            //region_.pop_back ();
            //data_.resize (data_.size () - info.So_count * 5);
            //throw bs_exception ("scal_2p_data_holder::add_spfn", "Could not compute residual saturation");
            // TODO: LOG
          }
      }
    else
      {
        BS_ASSERT (false && "UNSUPPORTED DATA PLACEMENT TYPE") (placement_info_.type);
      }
  }

  void
  scal_2p_data_holder::add_sof3 (sp_array_item_t const &data, t_long region_index, bool is_water)
  {
    typedef t_int   index_t;

    BS_ASSERT ((data->size () % 3) == 0) (data->size ());

    if (placement_info_.type == scal::data_placement::scal_data_placement_null ||
        placement_info_.type == scal::data_placement::sof3_spfn)
      {
        scal_region_info_t info;
        info.So_count   = (int)data->size () / 3;
        info.Sp_count   = -1;
        info.so_offset  = (int)data_->size ();
        info.sp_offset  = -1;

        data_->resize (data_->size () + info.So_count * 2);
        placement_info_.type = scal::data_placement::sof3_spfn;

        scal::data_placement::all_regions_t::place_sof3_data (data_, placement_info_, data, is_water);
        region_.push_back (info);
      }
    else if (placement_info_.type == scal::data_placement::spfn_sof3)
      {
        BS_ASSERT (region_index < region_.size ()) (region_index) (region_.size ());
        scal_region_info_t &info = region_[region_index];

        info.So_count   = (int)data->size () / 3;
        info.so_offset  = (int)data_->size ();

        data_->resize (data_->size () + info.So_count * 2);
        scal::data_placement::all_regions_t::place_sof3_data (data_, placement_info_, data, is_water);
        precalc (info, get_region_internal ((index_t)region_index), scal::data_placement::spfn_sof3, is_water);
        if (!check (get_region_internal ((index_t)region_index), is_water))
          {
            //region_.pop_back ();
            //data_.resize (data_.size () - info.So_count * 5);
            //throw bs_exception ("scal_2p_data_holder::add_sof3", "Could not compute residual saturation");
            // TODO: LOG
          }
      }
    else
      {
        BS_ASSERT (false && "UNSUPPORTED DATA PLACEMENT TYPE") (placement_info_.type);
      }
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
                if (i == region.Krp.size() - 1)
                  return false;
                else
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
        scal_region_t gas_region = get_region_internal ((int)i);
        const scal_region_t &water_region = water_data->get_region_internal ((int)i);

        item_t swc = water_region.Sp[0];
        for (size_t j = 0, jcnt = gas_region.Sp.size (); j < jcnt; ++j)
          {
            gas_region.So[(int)j] -= swc; // TODO: BUG:
          }
      }
  }

  //////////////////////////////////////////////////////////////////////////
  struct scal_3p::scal_3p_impl_base
  {
    typedef t_double                          item_t;
    typedef t_long                            index_t;
    typedef scal_3p                           scal_3p_t;

    virtual ~scal_3p_impl_base () {}

    virtual void
    get_relative_perm (index_t cell_index,
      const sp_array_item_t saturation,
      const sp_array_index_t sat_regions,
      sp_array_item_t relative_perm,
      sp_array_item_t s_deriv_relative_perm) const = 0;

    virtual void
    get_capillary (index_t cell_index,
      const sp_array_item_t saturation,
      const sp_array_index_t sat_regions,
      const sp_array_item_t perm,
      const sp_array_item_t poro,
      sp_array_item_t cap,
      sp_array_item_t s_deriv_cap) const = 0;

    virtual void
    process (const sp_array_item_t &saturation,
      const sp_array_index_t &sat_regions,
      const stdv_float &perm,
      const stdv_float &poro,
      data_array_t &data) const = 0;

    virtual void
    process_init (index_t i, const item_t *pressure, index_t sat_reg, const item_t *perm_array, item_t poro, item_t *sat, item_t *pc_limit) const = 0;

    virtual void
    calc_pcp (index_t cell_index, item_t sat, index_t sat_reg, item_t cap, item_t &pcp) const = 0;

    virtual void
    calc_gas_water_zone (index_t cell_index, index_t sat_reg, const item_t *perm_array, item_t poro, item_t pcgw, item_t &sw, item_t &sg) const = 0;

    virtual bool
    is_water () const = 0;

    virtual bool
    is_gas () const = 0;

    virtual bool
    is_oil () const = 0;
  };

  template <bool is_w, bool is_g, bool is_o, RPO_MODEL_ENUM rpo_model>
  struct scal_3p_impl : scal_3p::scal_3p_impl_base
  {
    typedef t_double                                      item_t;
    typedef t_long                                        index_t;
    typedef v_long                                        index_array_t;
    typedef v_double                                      item_array_t;
    typedef scal_3p                                       scal_3p_t;
    typedef scal_3p_t::phase_d_t                          phase_d_t;
    typedef scal_3p_t::data_array_t                       data_array_t;
    typedef scal_3p_t::data_t                             data_t;
    typedef scal_3p_t::sp_array_index_t                   sp_array_index_t;
    typedef scal_3p_t::sp_array_item_t                    sp_array_item_t;


    typedef scal_3p_t::scale_array_holder_t               scale_array_holder_t;
    typedef scal_3p_t::scal_2p_data_holder_t              scal_2p_data_holder_t;
    typedef scal_3p_t::scal_region_t                      scal_region_t;

    typedef scal_3p_t::sp_scale_array_holder_t            sp_scale_array_holder_t;
    typedef scal_3p_t::sp_scal_2p_data_holder_t           sp_scal_2p_data_holder_t;
    typedef scal_3p_t::sp_jfunction_t                     sp_jfunction_t;

    enum {
      n_phases = is_w + is_g + is_o,
    };

  public:
    scal_3p_impl (const sp_scale_array_holder_t &water_scale, const sp_scale_array_holder_t &gas_scale, 
      const sp_scal_2p_data_holder_t &water_data, const sp_scal_2p_data_holder_t &gas_data, 
      const sp_jfunction_t &water_jfunc, const sp_jfunction_t &gas_jfunc,
      const phase_d_t &phase_d, const phase_d_t &sat_d, const bool is_scalecrs_)
    : water_scale (water_scale),
    gas_scale (gas_scale),
    water_data (water_data),
    gas_data (gas_data),
    water_jfunc (water_jfunc),
    gas_jfunc (gas_jfunc),
    is_scalecrs (is_scalecrs_)
    {
      BOOST_STATIC_ASSERT (phase_d_t::static_size == (size_t)FI_PHASE_TOT);

      i_w = phase_d[FI_PHASE_WATER];
      i_g = phase_d[FI_PHASE_GAS];
      i_o = phase_d[FI_PHASE_OIL];

      i_w_w = i_w * n_phases + i_w;
      i_w_g = i_w * n_phases + i_g;
      i_w_o = i_w * n_phases + i_o;
      i_g_w = i_g * n_phases + i_w;
      i_g_g = i_g * n_phases + i_g;
      i_g_o = i_g * n_phases + i_o;
      i_o_w = i_o * n_phases + i_w;
      i_o_g = i_o * n_phases + i_g;
      i_o_o = i_o * n_phases + i_o;

      i_s_w = sat_d[FI_PHASE_WATER];
      i_s_g = sat_d[FI_PHASE_GAS];
      
    }

    void
    get_relative_perm (index_t cell_index,
      const sp_array_item_t saturation,
      const sp_array_index_t sat_regions,
      sp_array_item_t relative_perm, 
      sp_array_item_t s_deriv_relative_perm) const
    {
      if (cell_index >= (index_t)sat_regions->size ())
        {
          bs_throw_exception (boost::format ("Sat_regions size (%d) smaller than cell_index (%d)") % (index_t)sat_regions->size () % cell_index);
        }

      if (cell_index * n_phases >= (index_t)saturation->size ())
        {
          bs_throw_exception (boost::format ("Saturation size (%d) smaller than cell_index (%d) * n_phases (%d)") % (index_t)saturation->size () % cell_index % (index_t)n_phases);
        }

      if ((index_t)relative_perm->size () != n_phases)
        {
          bs_throw_exception (boost::format ("Size of relative_perm (%d) should be equal to n_phases (%d)") % (index_t)relative_perm->size () % (index_t)n_phases);
        }
      if ((index_t)s_deriv_relative_perm->size () != (n_phases * n_phases))
        {
          bs_throw_exception (boost::format ("Size of s_deriv_relative_perm (%d) should be equal to n_phases * n_phases (%d)") % (index_t)s_deriv_relative_perm->size () % (index_t)(n_phases * n_phases));
        }

      index_t sat_reg = (*sat_regions)[cell_index];
      process (cell_index, &(*saturation)[cell_index * n_phases], sat_reg, &(*relative_perm)[0], &(*s_deriv_relative_perm)[0]);
    }

    void
    get_capillary (index_t cell_index,
      const sp_array_item_t saturation,
      const sp_array_index_t sat_regions,
      const sp_array_item_t perm,
      const sp_array_item_t poro,
      sp_array_item_t cap,
      sp_array_item_t s_deriv_cap) const
    {
      if (cell_index >= (index_t)sat_regions->size ())
        {
          bs_throw_exception (boost::format ("Sat_regions size (%d) smaller than cell_index (%d)") % (index_t)sat_regions->size () % cell_index);
        }

      if (cell_index * n_phases >= (index_t)saturation->size ())
        {
          bs_throw_exception (boost::format ("Saturation size (%d) smaller than cell_index (%d) * n_phases (%d)") % (index_t)saturation->size () % cell_index % (index_t)n_phases);
        }

      if (cell_index * PLANE_ORIENTATION_TOTAL >= (index_t)perm->size ())
        {
          bs_throw_exception (boost::format ("Perm size (%d) smaller than cell_index (%d) * PLANE_ORIENTATION_TOTAL (%d)") % (index_t)perm->size () % cell_index % (index_t)PLANE_ORIENTATION_TOTAL);
        }
      if (cell_index >= (index_t)poro->size ())
        {
          bs_throw_exception (boost::format ("Poro size (%d) smaller than cell_index (%d)") % (index_t)poro->size () % cell_index);
        }

      if ((index_t)cap->size () != (n_phases - 1))
        {
          bs_throw_exception (boost::format ("Size of cap (%d) should be equal to n_phases - 1 (%d)") % (index_t)cap->size () % (index_t)(n_phases - 1));
        }
      if ((index_t)s_deriv_cap->size () != (n_phases - 1))
        {
          bs_throw_exception (boost::format ("Size of s_deriv_cap (%d) should be equal to n_phases - 1 (%d)") % (index_t)s_deriv_cap->size () % (index_t)(n_phases - 1));
        }

      index_t sat_reg = (*sat_regions)[cell_index];
      process_capillary (cell_index, &(*saturation)[cell_index * n_phases], sat_reg, &(*perm)[cell_index * PLANE_ORIENTATION_TOTAL], (*poro)[cell_index], &(*cap)[0], &(*s_deriv_cap)[0]);
    }

    void
    process (const spv_double &saturation,
      const spv_long &sat_regions,
      const stdv_float &perm,
      const stdv_float &poro,
      data_array_t &data_) const
    {
      for (index_t i = 0, cnt = (index_t) data_.size (); i < cnt; ++i)
        {
          data_t &data_i = data_[i];
          index_t sat_reg = (*sat_regions)[i];

          process (i, &(*saturation)[i * n_phases], sat_reg, &data_i.relative_perm[0], &data_i.s_deriv_relative_perm[0]);
          process_capillary (i, &(*saturation)[i * n_phases], sat_reg, &perm[i * PLANE_ORIENTATION_TOTAL], poro[i], &data_i.cap_pressure[0], &data_i.s_deriv_cap_pressure[0]);
        }
    }

    void
    process_init (index_t cell_index, const item_t *pressure, index_t sat_reg, const item_t *perm_array,
      item_t poro, item_t *sat, item_t *pc_limit) const
    {
      item_t cap, cap_min, cap_max;
      item_t mult;

      BS_ASSERT (pressure);
      BS_ASSERT (perm_array);
      BS_ASSERT (sat);
      BS_ASSERT (pc_limit);

      if (is_w)
        {
          BS_ASSERT (water_jfunc);
          if (!water_jfunc)
            throw bs_exception ("scal_3p::process_capillary", "water_jfunc is null");

          cap = pressure[i_w] - pressure[i_o];

          scal_region_t w_region = water_data->get_region (sat_reg);

          mult = 1;
          if (!water_jfunc->valid ())
            {
              item_t pcw_max    = w_region.get_pcp_max ();
              item_t pcw        = water_scale->get (pcp, pcw_max) [cell_index];

              if ((pcw_max * pcw) > EPS_DIFF)
                {
                  mult = pcw / pcw_max;
                  cap /= mult;
                }
            }

          process_init (cell_index, cap, *water_scale, w_region, perm_array, poro, water_jfunc, sat[i_s_w], cap_max, cap_min);

          pc_limit[i_s_w] = cap_min * mult;
        }
      if (is_g)
        {
          BS_ASSERT (gas_jfunc);
          if (!gas_jfunc)
            throw bs_exception ("scal_3p::process_capillary", "gas_jfunc is null");

          cap = pressure[i_g] - pressure[i_o];

          process_init (cell_index, cap, *gas_scale, gas_data->get_region (sat_reg), perm_array, poro, gas_jfunc, sat[i_s_g], cap_min, cap_max);

          pc_limit[i_s_g] = cap_max;
        }
    }

    void
    calc_pcp (index_t cell_index, item_t sat, index_t sat_reg, item_t cap, item_t &pcp) const
    {
      if (is_w)
        {
          item_t cap_table;

          scal_region_t w_region = water_data->get_region (sat_reg);

          if (w_region.get_pcp_max () < EPS_DIFF)
            {
              pcp = 0;
              return ;
            }

          item_t s_max  = w_region.get_phase_sat_max ();
          item_t s_min  = w_region.get_phase_sat_min ();

          item_t su     = water_scale->get (blue_sky::su, s_max) [cell_index];
          item_t sl     = water_scale->get (blue_sky::sl, s_min) [cell_index];

          item_t s      = scale_table (sl, s_min, su, s_max, sat);
          //item_t mult   = (s_max - s_min) / (su - sl);

          interpolate (w_region.Sp, w_region.Pcp, s, cap_table, std::less <item_t> ());

          if (-cap_table > EPS_DIFF)
            pcp = cap * w_region.get_pcp_max () / cap_table;
          else
            pcp = 0;
        }
      else if (is_g)
        {
          bs_throw_exception ("Error: calculating PCW for gas");
        }
    }

    void
    calc_gas_water_zone (index_t cell_index, index_t sat_reg, const item_t *perm_array, item_t poro, item_t pcgw, item_t &sw, item_t &sg) const
    {
      std::vector <item_t> pcgw_table;
      size_t n_table, i_table;
      item_t sat_cell[2], cap_cell[2];
      item_t sw_max, sw_min, swu, swl, s;

      if (!(is_w && is_g))
        {
          bs_throw_exception ("Error: calculating gas-water transition zone");
        }

      const scal_region_t w_region = water_data->get_region (sat_reg);
      const scal_region_t g_region = gas_data->get_region (sat_reg);

      n_table = w_region.Sp.size ();
      pcgw_table.resize (n_table, 0);

      sw_max  = w_region.get_phase_sat_max ();
      sw_min  = w_region.get_phase_sat_min ();
      swu     = water_scale->get (blue_sky::su, sw_max) [cell_index];
      swl     = water_scale->get (blue_sky::sl, sw_min) [cell_index];

      //complete gas-water cap pressure table
      for (i_table = 0; i_table < n_table; i_table++)
        {
          sat_cell[i_s_w] = scale_not_table (swl, sw_min, swu, sw_max, w_region.Sp[i_table]);
          sat_cell[i_s_g] = 1. - sat_cell[i_s_w];

          process_capillary (cell_index, sat_cell, sat_reg, perm_array, poro, cap_cell, 0);

          pcgw_table[i_table] = cap_cell[i_s_g] - cap_cell[i_s_w];
        }

      interpolate (pcgw_table, w_region.Sp, pcgw, s, std::greater <item_t> ());

      sw = scale_not_table (swl, sw_min, swu, sw_max, s);
      sg = (1.0 - sw > 0.) ? 1.0 - sw : 0.;
    }

  protected:
    void
    process (index_t cell_index, const item_t *sat, index_t sat_reg, item_t *kr, item_t *d_kr) const
    {
      BS_ASSERT (sat);
      BS_ASSERT (kr);
      BS_ASSERT (d_kr);

      if (is_w && n_phases == 2)
        {
          process_water (cell_index, sat, sat_reg, kr[i_w], d_kr[i_w_w], kr[i_o], d_kr[i_o_o]);
        }
      else if (is_g && n_phases == 2)
        {
          process_gas (cell_index, sat, sat_reg, kr[i_g], d_kr[i_g_g], kr[i_o], d_kr[i_o_o]);
        }
      else if (rpo_model == STONE2_MODEL)
        {
          process_3p_stone2 (cell_index,
                             sat, sat_reg,
                             kr[i_w], kr[i_g], kr[i_o],
                             d_kr[i_w_w], d_kr[i_g_g],
                             d_kr[i_o_w], d_kr[i_o_g], d_kr[i_o_o]);
        }
      else
        {
          process_3p_default (cell_index,
                              sat, sat_reg,
                              kr[i_w], kr[i_g], kr[i_o],
                              d_kr[i_w_w], d_kr[i_g_g],
                              d_kr[i_o_w], d_kr[i_o_g], d_kr[i_o_o]);
        }
    }

    void
    process_capillary (index_t cell_index, const item_t *sat, index_t sat_reg, const item_t *perm, item_t poro, item_t *cap, item_t *d_cap) const
    {
      BS_ASSERT (sat);
      BS_ASSERT (perm);
      BS_ASSERT (cap);

      if (is_w)
        {
          BS_ASSERT (water_jfunc);
          process_capillary (cell_index, sat[i_s_w], *water_scale, water_data->get_region (sat_reg), perm, poro, water_jfunc, cap[i_s_w], d_cap ? &d_cap[i_s_w] : 0);
        }
      if (is_g)
        {
          BS_ASSERT (gas_jfunc);
          process_capillary (cell_index, sat[i_s_g], *gas_scale, gas_data->get_region (sat_reg), perm, poro, gas_jfunc, cap[i_s_g], d_cap ? &d_cap[i_s_g] : 0);
        }
    }

    void 
    process_water (index_t cell_index, const item_t *sat, index_t sat_reg, item_t &kr, item_t &d_kr, item_t &kro, item_t &d_kro) const
    {
      BS_ASSERT (is_w);
      BS_ASSERT (i_s_w != -1);
      BS_ASSERT (sat);

      const scal_region_t &region = water_data->get_region (sat_reg);

      if (is_g)
        {
          const scal_region_t &g_region = gas_data->get_region (sat_reg);
          item_t sg_min = g_region.get_phase_sat_min ();
          item_t sgl = gas_scale->get (sl, sg_min) [cell_index];

          region.process_2phase (cell_index, sat[i_w], sat[i_o], *water_scale, sg_min, sgl, kr, d_kr, kro, d_kro, is_scalecrs);
        }
      else
        {
          region.process_2phase (cell_index, sat[i_w], sat[i_o], *water_scale, 0, 0, kr, d_kr, kro, d_kro, is_scalecrs);
        }
    }

    void 
    process_gas (index_t cell_index, const item_t *sat, index_t sat_reg, item_t &kr, item_t &d_kr, item_t &kro, item_t &d_kro) const
    {
      BS_ASSERT (is_g);
      BS_ASSERT (i_s_g != -1);
      BS_ASSERT (sat);

      const scal_region_t &region = gas_data->get_region (sat_reg);

      if (is_w)
        {
          scal_region_t w_region  = water_data->get_region (sat_reg);
          item_t sw_min = w_region.get_phase_sat_min ();
          item_t swl = water_scale->get (sl, sw_min) [cell_index];

          region.process_2phase (cell_index, sat[i_g], sat[i_o], *gas_scale, sw_min, swl, kr, d_kr, kro, d_kro, is_scalecrs);
        }
      else
        {
          region.process_2phase (cell_index, sat[i_g], sat[i_o], *gas_scale, 0, 0, kr, d_kr, kro, d_kro, is_scalecrs);
        }
    }

    void
    process_3p_default (index_t cell_index, const item_t *sat, index_t sat_reg, 
      item_t &krw, item_t &krg, item_t &kro,
      item_t &d_krww, item_t &d_krgg, item_t &d_krow, item_t &d_krog, item_t &d_kroo) const
    {
      BS_ASSERT (sat);

      item_t krow = 0.0, krog = 0.0;
      process_water (cell_index, sat, sat_reg, krw, d_krww, krow, d_krow);
      process_gas   (cell_index, sat, sat_reg, krg, d_krgg, krog, d_krog);

      scal_region_t region  = water_data->get_region (sat_reg);
      item_t swco         = region.get_phase_sat_min ();
      item_t swt          = (sat[i_w] - swco) > 0.0 ? (sat[i_w] - swco) : 0;
      item_t sat_swt      = sat[i_g] + swt;

      if (sat_swt > EPS_DIFF)
        {
          kro = (sat[i_g] * krog + swt * krow) / (sat_swt);
          if (kro < 0.0)
            {
              kro = 0.0;
            }
          else if (kro > 1.0)
            {
              kro = 1.0;
            }
          else
            {
              d_kroo  = (sat[i_g] * d_krog + swt * d_krow) / sat_swt;
              d_krow  = (krow - kro) / (sat_swt);
              d_krog  = (krog - kro) / (sat_swt);
            }
        }
      else
        {
          // TODO:
          kro         = krow;
          d_krow      = d_krow;
          d_kroo      = d_krow;
          d_krog      = d_krog;
        }
    }

    void 
    process_3p_stone2  (index_t cell_index, const item_t *sat, index_t sat_reg,
      item_t &krw, item_t &krg, item_t &kro, item_t &d_krww, item_t &d_krgg, item_t &d_krow, item_t &d_krog, item_t &d_kroo) const
    {
      BS_ASSERT (sat);

      item_t krow = 0.0, krog = 0.0;
      process_water (cell_index, sat, sat_reg, krw, d_krww, krow, d_krow);
      process_gas   (cell_index, sat, sat_reg, krg, d_krgg, krog, d_krow);

      scal_region_t region = water_data->get_region (sat_reg);
      item_t rporw = region.get_krop_max ();
      rporw = water_scale->get (krop, rporw) [cell_index];

      if (fabs (rporw) > EPS_DIFF)
        {
          kro     = rporw * (krow / rporw + krw) * (krog / rporw + krg) - rporw * (krw + krg);
          //d_krow  = 0.0;
          //d_krog  = 0.0;
          //d_kroo  = 0.0;

          // if oil relative permeability produced by 2-nd Stone's model is negative
          // automatically change value of k_ro to zero (see Eclipse technical description)
          if (kro < 0.0)
            {
              kro = 0.0;
            }
          else if (kro > 1.0)
            {
              kro = 1.0;
            }
          else
            {
              d_krow  = rporw * d_krww * (krog / rporw + krg - 1);
              d_krog  = rporw * d_krgg * (krow / rporw + krw - 1);
              d_kroo  = (d_krow * krog + d_krog * krow) / rporw + d_krow * krg + d_krog * krw;
            }
        }
      else
        {
          kro = (krow + krw) * (krog + krg) - (krw + krg);
          if (kro < 0.0)
            {
              kro = 0.0;
            }
          else if (kro > 1.0)
            {
              kro = 1.0;
            }
          else
            {
              d_krow  = d_krww * (krog + krg - 1);
              d_krog  = d_krgg * (krow + krw - 1);
              d_kroo  = (d_krow * krog + d_krog * krow) + d_krow * krg + d_krog * krw;
            }
        }
    }

    void 
    process_capillary (index_t cell_index, item_t sat, const scale_array_holder_t &scale_arrays, const scal_region_t &region,
      const item_t *perm_array, item_t poro, const sp_jfunction_t jfunc,
      item_t &cap, item_t *d_cap) const
    {
      region.process_capillary (cell_index, sat, scale_arrays, cap, d_cap);

      BS_ASSERT (jfunc);
      if (jfunc->valid ())
        {
          item_t perm = jfunc->get_perm (perm_array);
          item_t mult = 0;
          if (perm > EPS_DIFF)
            mult = jfunc->st_phase * pow (poro, jfunc->alpha) / pow (perm, jfunc->beta);
          cap         = cap * mult;
          if (d_cap)
            *d_cap       = *d_cap * mult;
        }
      else  // jfunc not valid, looking for PCW or PCG
        {
          item_t pcp_max    = region.get_pcp_max ();
          item_t pcp        = scale_arrays.get (blue_sky::pcp, pcp_max) [cell_index];

          if (fabs (pcp_max) > EPS_DIFF)
            {
              item_t mult   = pcp / pcp_max;
              cap = cap * mult;
              if (d_cap)
                *d_cap  = *d_cap * mult;
            }
        }  
    }

    void
    process_init (index_t cell_index, item_t cap, const scale_array_holder_t &scale_arrays,
      const scal_region_t &region,
      const item_t *perm_array, item_t poro, 
      const sp_jfunction_t jfunc,
      item_t &sat, item_t &pc_first, item_t &pc_last) const
    {
      BS_ASSERT (jfunc);

      pc_first = region.Pcp.front ();
      pc_last = region.Pcp.back ();

      if (jfunc->valid ())
        {
          item_t perm = jfunc->get_perm (perm_array);
          item_t mult = 0;
          cap = 0;
          if (perm > EPS_DIFF)
            {
              mult = jfunc->st_phase * pow (poro, jfunc->alpha) / pow (perm, jfunc->beta);
              cap = cap / mult;
            }
          pc_first    *= mult;
          pc_last     *= mult;
        }

      region.process_init (cell_index, cap, scale_arrays, sat);
    }

    item_t
    get_water_oil_sat (const item_t *sat) const
    {
      BS_ASSERT (i_s_w != -1);
      BS_ASSERT (i_s_g != -1);

      return is_w && is_g ? (1.0 - sat[i_s_w] - sat[i_s_g]) : (1.0 - sat[i_s_w]);
    }

    item_t
    get_gas_oil_sat (const item_t *sat) const
    {
      BS_ASSERT (i_s_w != -1);
      BS_ASSERT (i_s_g != -1);

      return is_w && is_g ? (1.0 - sat[i_s_w] - sat[i_s_g]) : (1.0 - sat[i_s_g]);
    }

    bool
    is_water () const
    {
      return is_w;
    }

    bool
    is_gas () const
    {
      return is_g;
    }

    bool
    is_oil () const
    {
      return is_o;
    }


  protected:

    sp_scale_array_holder_t   water_scale;
    sp_scale_array_holder_t   gas_scale;
    sp_scal_2p_data_holder_t  water_data;
    sp_scal_2p_data_holder_t  gas_data;

    sp_jfunction_t            water_jfunc;
    sp_jfunction_t            gas_jfunc;

    index_t i_w, i_g, i_o;
    index_t i_w_w, i_w_g, i_w_o;
    index_t i_g_w, i_g_g, i_g_o;
    index_t i_o_w, i_o_g, i_o_o;

    index_t i_s_w, i_s_g;
    bool is_scalecrs;
  };

  //////////////////////////////////////////////////////////////////////////
  scal_3p::scal_3p (bs_type_ctor_param param /* = NULL */)
  : water_data  (BS_KERNEL.create_object (scal_2p_data_holder_t::bs_type ()))
  , gas_data    (BS_KERNEL.create_object (scal_2p_data_holder_t::bs_type ()))
  , water_scale (BS_KERNEL.create_object (scale_array_holder_t::bs_type ()))
  , gas_scale   (BS_KERNEL.create_object (scale_array_holder_t::bs_type ()))
  , impl_ (0)
  {
  }

  scal_3p::scal_3p (const this_t& s)
  : bs_refcounter (s), scal_3p_iface (s)
  {
    bs_throw_exception ("NOT_IMPL_YET");
  }


  void
  scal_3p::process_init (index_t cell_index, const item_t *pressure, index_t sat_reg, const item_t *perm_array, item_t poro,
                                     item_t *sat, item_t *pc_limit) const
  {
    impl_->process_init (cell_index, pressure, sat_reg, perm_array, poro, sat, pc_limit);
  }

  void
  scal_3p::calc_pcp (index_t cell_index, const item_t sat, index_t sat_reg, item_t cap, item_t &pcp) const
  {
    impl_->calc_pcp (cell_index, sat, sat_reg, cap, pcp);
  }

  void
  scal_3p::calc_gas_water_zone (index_t cell_index, index_t sat_reg, const item_t *perm_array, item_t poro, item_t pcgw,
      item_t &sw, item_t &sg) const
  {
    impl_->calc_gas_water_zone (cell_index, sat_reg, perm_array, poro, pcgw, sw, sg);
  }

  void
  scal_3p::init (bool is_w, bool is_g, bool is_o, const phase_d_t &phase_d, const phase_d_t &sat_d, RPO_MODEL_ENUM r, bool is_scalecrs_)
  {
    if (!water_jfunc)
      {
        bs_throw_exception ("Water jfunction not set");
      }

    if (!gas_jfunc)
      {
        bs_throw_exception ("Gas jfunction not set");
      }

    if (r == STONE2_MODEL)
      {
        if (is_w && is_g && is_o)
          impl_ = new scal_3p_impl <true, true, true, STONE2_MODEL> (water_scale, gas_scale, water_data, gas_data, water_jfunc, gas_jfunc, phase_d, sat_d, is_scalecrs_);
        else if (is_w && is_o)
          impl_ = new scal_3p_impl <true, false, true, STONE2_MODEL> (water_scale, gas_scale, water_data, gas_data, water_jfunc, gas_jfunc, phase_d, sat_d, is_scalecrs_);
        else if (is_g && is_o)
          impl_ = new scal_3p_impl <false, true, true, STONE2_MODEL> (water_scale, gas_scale, water_data, gas_data, water_jfunc, gas_jfunc, phase_d, sat_d, is_scalecrs_);
        else if (is_w)
          impl_ = new scal_3p_impl <true, false, false, STONE2_MODEL> (water_scale, gas_scale, water_data, gas_data, water_jfunc, gas_jfunc, phase_d, sat_d, is_scalecrs_);
        else if (is_g)
          impl_ = new scal_3p_impl <false, true, false, STONE2_MODEL> (water_scale, gas_scale, water_data, gas_data, water_jfunc, gas_jfunc, phase_d, sat_d, is_scalecrs_);
        else if (is_o)
          impl_ = new scal_3p_impl <false, false, true, STONE2_MODEL> (water_scale, gas_scale, water_data, gas_data, water_jfunc, gas_jfunc, phase_d, sat_d, is_scalecrs_);
        else
          {
            bs_throw_exception ("Unkown phase value");
          }
      }
    else
      {
        if (is_w && is_g && is_o)
          impl_ = new scal_3p_impl <true, true, true, RPO_DEFAULT_MODEL> (water_scale, gas_scale, water_data, gas_data, water_jfunc, gas_jfunc, phase_d, sat_d, is_scalecrs_);
        else if (is_w && is_o)
          impl_ = new scal_3p_impl <true, false, true, RPO_DEFAULT_MODEL> (water_scale, gas_scale, water_data, gas_data, water_jfunc, gas_jfunc, phase_d, sat_d, is_scalecrs_);
        else if (is_g && is_o)
          impl_ = new scal_3p_impl <false, true, true, RPO_DEFAULT_MODEL> (water_scale, gas_scale, water_data, gas_data, water_jfunc, gas_jfunc, phase_d, sat_d, is_scalecrs_);
        else if (is_w)
          impl_ = new scal_3p_impl <true, false, false, RPO_DEFAULT_MODEL> (water_scale, gas_scale, water_data, gas_data, water_jfunc, gas_jfunc, phase_d, sat_d, is_scalecrs_);
        else if (is_g)
          impl_ = new scal_3p_impl <false, true, false, RPO_DEFAULT_MODEL> (water_scale, gas_scale, water_data, gas_data, water_jfunc, gas_jfunc, phase_d, sat_d, is_scalecrs_);
        else if (is_o)
          impl_ = new scal_3p_impl <false, false, true, RPO_DEFAULT_MODEL> (water_scale, gas_scale, water_data, gas_data, water_jfunc, gas_jfunc, phase_d, sat_d, is_scalecrs_);
        else
          {
            bs_throw_exception ("Unkown phase value");
          }
      }

    water_data->init_regions ();
    gas_data->init_regions ();
  }


  void 
  scal_3p::set_water_jfunction (sp_jfunction_t jfunc)
  {
    water_jfunc = jfunc;
  }

  void 
  scal_3p::set_gas_jfunction (sp_jfunction_t jfunc)
  {
    gas_jfunc = jfunc;
  }

  void
  scal_3p::update_gas_data ()
  {
    if (impl_->is_water () && impl_->is_gas ())
      gas_data->update_so (water_data);
  }

  void
  scal_3p::process (const spv_double &saturation,
    const spv_long &sat_regions,
    const stdv_float &perm,
    const stdv_float &poro,
    data_array_t &data) const
  {
    impl_->process (saturation, sat_regions, perm, poro, data);
  }

  void
  scal_3p::get_relative_perm (index_t cell_index,
    const sp_array_item_t saturation,
    const sp_array_index_t sat_regions,
    sp_array_item_t relative_perm,
    sp_array_item_t s_deriv_relative_perm) const
  {
    impl_->get_relative_perm (cell_index, saturation, sat_regions, relative_perm, s_deriv_relative_perm);
  }

  void
  scal_3p::get_capillary (index_t cell_index,
    const sp_array_item_t saturation,
    const sp_array_index_t sat_regions,
    const sp_array_item_t perm,
    const sp_array_item_t poro,
    sp_array_item_t cap,
    sp_array_item_t s_deriv_cap) const
  {
    impl_->get_capillary (cell_index, saturation, sat_regions, perm, poro, cap, s_deriv_cap);
  }

  scal_3p::~scal_3p ()
  {
    delete impl_;
  }

  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE (scale_array_holder);
  BLUE_SKY_TYPE_STD_COPY (scale_array_holder);
  BLUE_SKY_TYPE_IMPL(scale_array_holder, scale_array_holder_iface, "scale_array_holder", "scale_array_holder calculation class", "scale_array_holder calculation");

  BLUE_SKY_TYPE_STD_CREATE (scal_2p_data_holder);
  BLUE_SKY_TYPE_STD_COPY (scal_2p_data_holder);
  BLUE_SKY_TYPE_IMPL(scal_2p_data_holder, scal_2p_data_holder_iface, "scal_2p_data_holder", "scal_2p_data_holder calculation class", "scal_2p_data_holder calculation");

  BLUE_SKY_TYPE_STD_CREATE (scal_3p);
  BLUE_SKY_TYPE_STD_COPY (scal_3p);
  BLUE_SKY_TYPE_IMPL(scal_3p, scal_3p_iface, "scal_3p", "scal_3p calculation class", "scal_3p calculation");
  //////////////////////////////////////////////////////////////////////////
  bool scal_register_types (const blue_sky::plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, scale_array_holder::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, scal_2p_data_holder::bs_type ()); BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, scal_3p::bs_type ()); BS_ASSERT (res);

    return res;
  }
} // namespace blue_sky
