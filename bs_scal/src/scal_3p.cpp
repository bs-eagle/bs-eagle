/**
 * \file scal_3p.cpp
 * \brief
 * \author Sergey Miryanov
 * \date 19.05.2008
 * */
#include "bs_scal_stdafx.h"

#include "scal_3p.h"
#include "scal_3p_impl.h"
#include "scal_data_source.h"
#include "scal_data_vector.h"
#include "scale_arrays_placement_strategies.h"
#include "scal_data_placement_strategies.h"

#include "scale_array_holder.h"
#include "scal_region_info.h"
#include "scal_region.h"
#include "scal_2p_data_holder.h"

#include "jfunction.h"

#include <sstream>

namespace blue_sky
  {
    BS_TYPE_IMPL_T_EXT_MEM(bs_array, 2, (calc_model_data, bs_array_shared));
    BS_TYPE_IMPL_T_EXT_MEM(bs_array, 2, (calc_model_data, vector_traits));

  //////////////////////////////////////////////////////////////////////////

  scale_array_holder::scale_array_holder (bs_type_ctor_param /*param  = NULL */)
  {
    data_pool = BS_KERNEL.create_object ("float_var_table");
    BS_ASSERT (data_pool);

    data_pool->init ( scale_array_name_total);
    
    data_pool->set_col_name (socr, L"SOCR");
    data_pool->set_col_name (scr, L"SCR");
    data_pool->set_col_name (su, L"SU");
    data_pool->set_col_name (sl, L"SL");
    data_pool->set_col_name (pcp, L"PCP");
    data_pool->set_col_name (krp, L"KRP");
    data_pool->set_col_name (krop, L"KROP");
    data_pool->set_col_name (krpr, L"KRPR");
    data_pool->set_col_name (krorp, L"KRORP");
    BOOST_STATIC_ASSERT (krorp == scale_array_name_total - 1);
  }

  scale_array_holder::scale_array_holder (const scale_array_holder& s)
  : bs_refcounter (s), scale_array_holder_iface (s)
  {
    BS_ASSERT (false && "NOT IMPL YET");
  }

  //////////////////////////////////////////////////////////////////////////
  scal_3p::scal_3p (bs_type_ctor_param /*param  = NULL */)
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
  scal_3p::process_init_2 (const item_t *pressure, index_t sat_reg, item_t perm, item_t poro,
                           item_t *sat, item_t *pc_limit) const
  {
    impl_->process_init_2 (pressure, sat_reg, perm, poro, sat, pc_limit);
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
  scal_3p::calc_gas_water_zone_2 (index_t sat_reg, item_t perm, item_t poro, item_t pcgw,
      item_t &sw, item_t &sg) const
  {
    impl_->calc_gas_water_zone_2 (sat_reg, perm, poro, pcgw, sw, sg);
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
          impl_ = new scal_3p_impl (true, true, true, STONE2_MODEL, water_scale, gas_scale, water_data, gas_data, water_jfunc, gas_jfunc, phase_d, sat_d, is_scalecrs_);
        else if (is_w && is_o)
          impl_ = new scal_3p_impl (true, false, true, STONE2_MODEL, water_scale, gas_scale, water_data, gas_data, water_jfunc, gas_jfunc, phase_d, sat_d, is_scalecrs_);
        else if (is_g && is_o)
          impl_ = new scal_3p_impl (false, true, true, STONE2_MODEL, water_scale, gas_scale, water_data, gas_data, water_jfunc, gas_jfunc, phase_d, sat_d, is_scalecrs_);
        else if (is_w)
          impl_ = new scal_3p_impl (true, false, false, STONE2_MODEL, water_scale, gas_scale, water_data, gas_data, water_jfunc, gas_jfunc, phase_d, sat_d, is_scalecrs_);
        else if (is_g)
          impl_ = new scal_3p_impl (false, true, false, STONE2_MODEL, water_scale, gas_scale, water_data, gas_data, water_jfunc, gas_jfunc, phase_d, sat_d, is_scalecrs_);
        else if (is_o)
          impl_ = new scal_3p_impl (false, false, true, STONE2_MODEL, water_scale, gas_scale, water_data, gas_data, water_jfunc, gas_jfunc, phase_d, sat_d, is_scalecrs_);
        else
          {
            bs_throw_exception ("Unkown phase value");
          }
      }
    else
      {
        if (is_w && is_g && is_o)
          impl_ = new scal_3p_impl (true, true, true, RPO_DEFAULT_MODEL, water_scale, gas_scale, water_data, gas_data, water_jfunc, gas_jfunc, phase_d, sat_d, is_scalecrs_);
        else if (is_w && is_o)
          impl_ = new scal_3p_impl (true, false, true, RPO_DEFAULT_MODEL, water_scale, gas_scale, water_data, gas_data, water_jfunc, gas_jfunc, phase_d, sat_d, is_scalecrs_);
        else if (is_g && is_o)
          impl_ = new scal_3p_impl (false, true, true, RPO_DEFAULT_MODEL, water_scale, gas_scale, water_data, gas_data, water_jfunc, gas_jfunc, phase_d, sat_d, is_scalecrs_);
        else if (is_w)
          impl_ = new scal_3p_impl (true, false, false, RPO_DEFAULT_MODEL, water_scale, gas_scale, water_data, gas_data, water_jfunc, gas_jfunc, phase_d, sat_d, is_scalecrs_);
        else if (is_g)
          impl_ = new scal_3p_impl (false, true, false, RPO_DEFAULT_MODEL, water_scale, gas_scale, water_data, gas_data, water_jfunc, gas_jfunc, phase_d, sat_d, is_scalecrs_);
        else if (is_o)
          impl_ = new scal_3p_impl (false, false, true, RPO_DEFAULT_MODEL, water_scale, gas_scale, water_data, gas_data, water_jfunc, gas_jfunc, phase_d, sat_d, is_scalecrs_);
        else
          {
            bs_throw_exception ("Unkown phase value");
          }
      }
    // TODO: see initialization from tables 
    //water_data->init_regions ();
    //gas_data->init_regions ();
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

  void
  scal_3p::init_scal_calc_data ()
  {
	  sp_jfunction_t wat_jfunc = BS_KERNEL.create_object(jfunction::bs_type());
	  sp_jfunction_t gas_jfunc = BS_KERNEL.create_object(jfunction::bs_type());
	  phase_d_t phase_d;
	  sat_d_t sat_d;
	  t_long n_phases = 0;
      t_long phases = 0;
	  //t_long sat_counter = 0;    
      if (is_water)
      {
        phase_d[0] = n_phases++;
        phases |= 1 << FI_PHASE_WATER;
      }
      else 
        phase_d[0] = -1;

      if (is_gas)
      {
        phase_d[1] = n_phases++;
        phases |= 1 << FI_PHASE_GAS;
      }
      else 
        phase_d[1] = -1;

      if (is_oil)
      {
        phase_d[2] = n_phases++;
        phases |= 1 << FI_PHASE_OIL;
      }
      else 
        phase_d[2] = -1;

      for (size_t i = 0, sat_counter = 0; i < FI_PHASE_TOT; ++i)
      {
        if ((phases & (1 << i)) && ((t_long)sat_counter < (t_long)n_phases ))
          sat_d[i] = sat_counter++;
        else
          sat_d[i] = -1;
      }
	  init_from_scal_ex(phase_d, sat_d, wat_jfunc, gas_jfunc);
  }

  void
  scal_3p::init_from_scal_ex(const phase_d_t &phase_d, const sat_d_t &sat_d,
						  sp_jfunction_t water_jfunc,
						  sp_jfunction_t gas_jfunc)
  {
	  init_scal_data_from_input_tables();

	  set_water_jfunction(water_jfunc);
	  set_gas_jfunction(gas_jfunc);
	  init(is_water, is_gas, is_oil, phase_d, sat_d, RPO_DEFAULT_MODEL);
	  update_gas_data();
  }
  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE (scale_array_holder);
  BLUE_SKY_TYPE_STD_COPY (scale_array_holder);
  BLUE_SKY_TYPE_IMPL(scale_array_holder, scale_array_holder_iface, "scale_array_holder", "scale_array_holder calculation class", "scale_array_holder calculation");

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
  
  
  void
  scal_3p::init_scal_input_table_arrays (const t_long n_scal_regions_, 
                                         bool is_oil_, bool is_gas_, bool is_water_)
    {
      BS_ASSERT (n_scal_regions_ > 0);
      n_scal_regions = n_scal_regions_;
      is_oil = is_oil_;
	    is_gas = is_gas_;
	    is_water = is_water_;

      water_input_table.resize (n_scal_regions);
      gas_input_table.resize (n_scal_regions);
      oil_input_table.resize (n_scal_regions);
      
      for (t_long i = 0; i < n_scal_regions; ++i)
        {
          if (is_water && is_oil)
            {
              water_input_table[i] = BS_KERNEL.create_object ("table");
            }
          if (is_gas && is_oil)
            {
              gas_input_table[i] = BS_KERNEL.create_object ("table");
            }  
#if 0          
          if (is_oil && (is_water || is_gas)) 
            {
              oil_input_table[i] = BS_KERNEL.create_object ("table");
            }
#endif             
        }
    }                                 
  
  std::list <BS_SP (table_iface)>
  scal_3p::get_tables_list (t_long index_scal_region) const
    {
      BS_ASSERT (index_scal_region >= 0 && index_scal_region < n_scal_regions);
      std::list<BS_SP( table_iface)> tables;
      
      if (water_input_table[index_scal_region])
        tables.push_back (water_input_table[index_scal_region]);
      if (gas_input_table[index_scal_region])
        tables.push_back (gas_input_table[index_scal_region]);
      if (oil_input_table[index_scal_region])
        tables.push_back (oil_input_table[index_scal_region]);
      return tables;
    } 
       
       
  std::list <BS_SP (table_iface)>
  scal_3p::get_tables_fluid_all_regions (t_long scal_fluid_type) const
    {
      BS_ASSERT (scal_fluid_type >= FI_PHASE_NULL && scal_fluid_type < FI_PHASE_TOT);
      std::list<BS_SP( table_iface)> tables;

      if (scal_fluid_type == FI_PHASE_WATER)
        {
          for (t_long i = 0; i < n_scal_regions; ++i)
            tables.push_back (water_input_table[i]);
        }
      else if (scal_fluid_type == FI_PHASE_GAS)
        {
          for (t_long i = 0; i < n_scal_regions; ++i)
            tables.push_back (gas_input_table[i]);
        }  
      else if (scal_fluid_type == FI_PHASE_OIL)
        {  
          for (t_long i = 0; i < n_scal_regions; ++i)
            tables.push_back (oil_input_table[i]);
        }     
      
      return tables;
    } 
        
  BS_SP (table_iface) 
  scal_3p::get_table (t_int scal_fluid_type, t_long index_scal_region = 0) const
    {
      BS_ASSERT (index_scal_region >= 0 && index_scal_region < n_scal_regions);
      BS_ASSERT (scal_fluid_type >= FI_PHASE_NULL && scal_fluid_type < FI_PHASE_TOT);
    
      if (scal_fluid_type == FI_PHASE_WATER)
        {
          return water_input_table[index_scal_region]; 
        }
      else if (scal_fluid_type == FI_PHASE_GAS)
        {
          return gas_input_table[index_scal_region]; 
        }  
      //else if (scal_fluid_type == FI_PHASE_OIL)
      else
        {  
          return oil_input_table[index_scal_region];
        }     
    }                                 

  void 
  scal_3p::init_scal_data_from_input_tables ()
    {
      BS_ASSERT (n_scal_regions > 0);
      
	    get_water_data ()->clear_regions();
	    get_gas_data ()->clear_regions();

      get_water_data ()->init_table_array (n_scal_regions);
      get_gas_data ()->init_table_array (n_scal_regions);
      
      for (t_long region_index = 0; region_index < n_scal_regions; ++region_index)
        {
          if (water_input_table[region_index].get ())
            {
              t_long n_cols_water = water_input_table[region_index]->get_n_cols ();
              t_long n_rows_water = water_input_table[region_index]->get_n_rows ();
              if (n_cols_water == SPOF_KEYWORD_COLUMNS)
                {
                  get_water_data ()->add_spof (water_input_table[region_index]->convert_to_array (n_rows_water, n_cols_water), region_index, true);
                }
              else 
                {
                  get_water_data ()->add_spfn (water_input_table[region_index]->convert_to_array (n_rows_water, n_cols_water), region_index, true);
                  if (oil_input_table[region_index].get ())
                    {
                      t_long n_cols_oil = oil_input_table[region_index]->get_n_cols ();
                      t_long n_rows_oil = oil_input_table[region_index]->get_n_rows ();
                      if (n_cols_oil == SOF3_KEYWORD_COLUMNS)
                        {
                          get_water_data ()->add_sof3 (oil_input_table[region_index]->convert_to_array (n_rows_oil, n_cols_oil), region_index, true);
                        }
                      else 
                        {
                          get_water_data ()->add_sof2 (oil_input_table[region_index]->convert_to_array (n_rows_oil, n_cols_oil), region_index, true);
                        }  
                    }  
                }  
            }
          
          if (gas_input_table[region_index].get ())
            {
              t_long n_cols_gas = gas_input_table[region_index]->get_n_cols ();
              t_long n_rows_gas = gas_input_table[region_index]->get_n_rows ();
              if (n_cols_gas == SPOF_KEYWORD_COLUMNS)
                {
                  get_gas_data ()->add_spof (gas_input_table[region_index]->convert_to_array (n_rows_gas, n_cols_gas), region_index, false);
                }
              else 
                {
                  get_gas_data ()->add_spfn (gas_input_table[region_index]->convert_to_array (n_rows_gas, n_cols_gas), region_index, false);
                  if (oil_input_table[region_index].get ())
                    {
                      t_long n_cols_oil = oil_input_table[region_index]->get_n_cols ();
                      t_long n_rows_oil = oil_input_table[region_index]->get_n_rows ();
                    
                      if (oil_input_table[region_index]->get_n_cols () == SOF3_KEYWORD_COLUMNS)
                        {
                          get_gas_data ()->add_sof3 (oil_input_table[region_index]->convert_to_array (n_rows_oil, n_cols_oil), region_index, false);
                        }
                      else 
                        {
                          get_gas_data ()->add_sof2 (oil_input_table[region_index]->convert_to_array (n_rows_oil, n_cols_oil), region_index, false);
                        }  
                    }  
                }  
            }  
        }
      get_water_data ()->init_regions_from_tables ();
	  if (is_gas)
		  get_gas_data ()->init_regions_from_tables ();
        
    }  
} // namespace blue_sky
