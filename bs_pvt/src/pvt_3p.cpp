#include "bs_pvt_stdafx.h"

//#include BS_FORCE_PLUGIN_IMPORT ()
#include "pvt_3p.h"
#include "pvt_base.h"
#include "pvt_dead_oil.h"
#include "pvt_oil.h"
#include "pvt_gas.h"
#include "pvt_water.h"
#include "data_class.h"
#include "pvt_3p_iface.h"
//#include BS_STOP_PLUGIN_IMPORT ()

    
namespace blue_sky
{

  pvt_3p::pvt_3p (bs_type_ctor_param param /* = NULL */)
  { 
  
  }
  
  pvt_3p::pvt_3p (const pvt_3p& s)
  : bs_refcounter (s), pvt_3p_iface (s)
  {
    BS_ASSERT (false && "NOT IMPL YET");
  }

  /**
   * \class pvt_helper
   * \todo  describe
   * */
  struct pvt_helper
    {
      template <typename pvt_array_t, typename pvt_vector_t>
      static void
      set_array (const pvt_array_t &pvt, const pvt_vector_t &v)
      {
        for (size_t i = 0; i < v.size(); ++i)
          {
            const typename pvt_array_t::value_type &pvt__(pvt[i]);

            const typename pvt_vector_t::value_type &p = v[i];
            pvt__->insert_vector(*p.main_data_);
            if (p.has_density_)
              {
                pvt__->set_density (p.density_, p.molar_density_);
              }
          }
      }
    };
  
  /*struct pvt_table_helper
    {
      template <typename pvt_array_t, typename pvt_table_array_t>
      static void 
      set_array (const pvt_array_t &pvt, const pvt_table_array_t &t)
        {
          BS_ASSERT (pvt.size () == t.size ());
          for (size_t i = 0; i < t.size(); ++i)
            {
              const typename pvt_array_t::value_type &pvt__(pvt[i]);
              const typename pvt_table_array_t::value_type &p = t[i];
              
//              pvt__->pvt_input_props->copy (p->get_table ());
            }
        }
    };*/
      
  /*void
  pvt_3p::init_pvt_arrays (const t_long n_pvt_regions_, const sp_idata_t idata_,
                           bool is_oil, bool is_gas, bool is_water, 
                           t_float atm_p, t_float min_p, t_float max_p, t_float n_intervals)
  {
    typedef idata_t::pvt_vector    pvt_vector;
    typedef idata_t::pvt_info      pvt_info;

    n_pvt_regions = n_pvt_regions_;
    BS_ASSERT (n_pvt_regions > 0);
    
    pvt_oil_array.resize (n_pvt_regions);
    pvt_gas_array.resize (n_pvt_regions);
    pvt_water_array.resize (n_pvt_regions);

    for (size_t i = 0; i < n_pvt_regions; i++)
      {
        BS_ASSERT (idata_->pvto.size ());
        std::cout << "pvto (" << i << "): " << idata_->pvto.back ().main_data_->empty () << std::endl;

        if (idata_->pvto.back ().main_data_->empty ())
          {
            pvt_oil_array[i] = BS_KERNEL.create_object (pvt_dead_oil_t::bs_type());
          }
        else
          {
            pvt_oil_array[i] = BS_KERNEL.create_object (pvt_oil_t::bs_type());
          }

        pvt_gas_array[i] = BS_KERNEL.create_object (pvt_gas_t::bs_type());
        pvt_water_array[i] = BS_KERNEL.create_object (pvt_water_t::bs_type());
      }

    BS_ASSERT (idata_->pvto.size ());
    if (idata_->pvto.back ().main_data_->empty ())
      {
        pvt_helper::set_array (pvt_oil_array, idata_->pvtdo);
      }
    else
      {
        pvt_helper::set_array (pvt_oil_array, idata_->pvto);
      }

    pvt_helper::set_array (pvt_gas_array, idata_->pvtg);
    pvt_helper::set_array (pvt_water_array, idata_->pvtw);

    for (size_t i = 0; i < n_pvt_regions; i++)
      {
        if (is_oil)
          pvt_oil_array[i]->build (atm_p, min_p, max_p, n_intervals);

        if (is_gas)
          pvt_gas_array[i]->build (atm_p, min_p, max_p, n_intervals);

        if (is_water)
          pvt_water_array[i]->build (atm_p, min_p, max_p, n_intervals);
      }
  }*/


  /*void
  pvt_3p::init_pvt_arrays (const t_long n_pvt_regions_, 
                           const sp_pvt_dummy_iface_array_t &pvt_oil_data,
                           const sp_pvt_dummy_iface_array_t &pvt_gas_data,
                           const sp_pvt_dummy_iface_array_t &pvt_water_data,
                           bool is_oil, bool is_gas, bool is_water, 
                           t_float atm_p, t_float min_p, t_float max_p, t_float n_intervals)
  {
    typedef idata_t::pvt_vector    pvt_vector;
    typedef idata_t::pvt_info      pvt_info;

    n_pvt_regions = n_pvt_regions_;
    BS_ASSERT (n_pvt_regions > 0);
    if (is_oil)
      BS_ASSERT (pvt_oil_data.size () == n_pvt_regions);
    if (is_gas)
      BS_ASSERT (pvt_gas_data.size () == n_pvt_regions);
    if (is_water)
      BS_ASSERT (pvt_water_data.size () == n_pvt_regions);
    
    pvt_oil_array.resize (n_pvt_regions);
    pvt_gas_array.resize (n_pvt_regions);
    pvt_water_array.resize (n_pvt_regions);

    for (size_t i = 0; i < n_pvt_regions; i++)
      {
        if (is_oil) 
          {
            if (!is_gas)
              {
                pvt_oil_array[i] = BS_KERNEL.create_object (pvt_dead_oil_t::bs_type());
              }
            else
              {
                pvt_oil_array[i] = BS_KERNEL.create_object (pvt_oil_t::bs_type());
              }
          }
        if (is_gas)
          pvt_gas_array[i] = BS_KERNEL.create_object (pvt_gas_t::bs_type());
        if (is_water)  
          pvt_water_array[i] = BS_KERNEL.create_object (pvt_water_t::bs_type());
      }
    
    if (is_oil)
      {
        if (!is_gas)
          {
            pvt_table_helper::set_array (pvt_oil_array, pvt_oil_data);
          }
        else
          {
            pvt_table_helper::set_array (pvt_oil_array, pvt_oil_data);
          }
      }
    
    if (is_gas)  
      pvt_table_helper::set_array (pvt_gas_array, pvt_gas_data);
    if (is_water)
      pvt_table_helper::set_array (pvt_water_array, pvt_water_data);

    for (size_t i = 0; i < n_pvt_regions; i++)
      {
        if (is_oil)
          pvt_oil_array[i]->build (atm_p, min_p, max_p, n_intervals);

        if (is_gas)
          pvt_gas_array[i]->build (atm_p, min_p, max_p, n_intervals);

        if (is_water)
          pvt_water_array[i]->build (atm_p, min_p, max_p, n_intervals);
      }
  }*/
  
  struct init_pvt_arr_helper
  {
	  template <typename pvt_base_elem_t>
	  static void
	  set_pvt_base (pvt_base_elem_t elem, BS_SP(table_iface) table)
	  {
		  int n_cols = table->get_n_cols();
		  int n_rows = table->get_n_rows();
          spv_double data = BS_KERNEL.create_object (v_double::bs_type ());
		  data->resize(n_rows * n_cols);
		  for (size_t row = 0; row < n_rows; row++)
		  {
			  for (size_t col = 0; col < n_cols; col++)
			  {
				  (*data)[row*n_cols+col] = table->get_value(row, col);
			  }
		  }
		  elem->insert_vector(*data);
	  }
  };
  
  void
  pvt_3p::init_pvt_arrays (const t_long n_pvt_regions_, 
                           bool is_oil, bool is_gas, bool is_water)
  {

    this->n_pvt_regions = n_pvt_regions_;

    BS_ASSERT (n_pvt_regions > 0);
    
    pvt_oil_array.resize (n_pvt_regions);
    pvt_gas_array.resize (n_pvt_regions);
    pvt_water_array.resize (n_pvt_regions);
    density_table.resize (n_pvt_regions);

    density = BS_KERNEL.create_object (v_float::bs_type ());
    density->resize (n_pvt_regions * FI_PHASE_TOT);
    
    for (size_t i = 0; i < n_pvt_regions; i++)
      {
        density_table[i] = BS_KERNEL.create_object ("table");
        density_table[i]->init (1, FI_PHASE_TOT);
        if (is_oil) 
          {
            if (!is_gas)
              {
                pvt_oil_array[i] = BS_KERNEL.create_object (pvt_dead_oil_t::bs_type());
              }
            else
              {
                pvt_oil_array[i] = BS_KERNEL.create_object (pvt_oil_t::bs_type());
              }
          }
        if (is_gas)
          pvt_gas_array[i] = BS_KERNEL.create_object (pvt_gas_t::bs_type());
        if (is_water)  
          pvt_water_array[i] = BS_KERNEL.create_object (pvt_water_t::bs_type());
      }
  }
  

  BS_SP (pvt_dead_oil) 
  pvt_3p::get_pvt_oil (const t_long index_pvt_region) const
  {
    BS_ASSERT (index_pvt_region >= 0 && index_pvt_region < n_pvt_regions && pvt_oil_array.size () > 0 && index_pvt_region < pvt_oil_array.size ());
    return pvt_oil_array[index_pvt_region];
  }
  
  BS_SP (pvt_gas)
  pvt_3p::get_pvt_gas (const t_long index_pvt_region) const
  {
    BS_ASSERT (index_pvt_region >= 0 && index_pvt_region < n_pvt_regions && pvt_gas_array.size () > 0 && index_pvt_region < pvt_gas_array.size ());
    return pvt_gas_array[index_pvt_region];
  }  
  
  BS_SP (pvt_water)
  pvt_3p::get_pvt_water (const t_long index_pvt_region) const
  {
    BS_ASSERT (index_pvt_region >= 0 && index_pvt_region < n_pvt_regions && pvt_water_array.size () > 0 && index_pvt_region < pvt_water_array.size ());
    return pvt_water_array[index_pvt_region];
  }

  std::list <BS_SP( table_iface)>
  pvt_3p::get_tables_list (t_long pvt_fluid_type) const
  {
    std::list<BS_SP( table_iface)> tables;
    
    if (pvt_fluid_type == FI_PHASE_OIL)
      {
        for (t_long i = 0; i < n_pvt_regions; ++i)
          tables.push_back (get_pvt_oil (i)->get_pvt_input_table ());
      }
    else if (pvt_fluid_type == FI_PHASE_WATER)
      {
        for (t_long i = 0; i < n_pvt_regions; ++i)
          tables.push_back (get_pvt_water (i)->get_pvt_input_table ());
      }  
    else if (pvt_fluid_type == FI_PHASE_GAS)
      {
        for (t_long i = 0; i < n_pvt_regions; ++i)
          tables.push_back (get_pvt_gas (i)->get_pvt_input_table ());
      }  
    return tables;
  }

  //! return input table for fluid type of for defined pvt region 
  BS_SP (table_iface)
  pvt_3p::get_table (t_long index_pvt_region, t_long pvt_fluid_type) const 
  {
    BS_ASSERT (index_pvt_region >= 0 && index_pvt_region < n_pvt_regions);
    BS_ASSERT (pvt_fluid_type >= FI_PHASE_NULL && pvt_fluid_type < FI_PHASE_TOT);
  
    if (pvt_fluid_type == FI_PHASE_OIL)
      {
        return get_pvt_oil (index_pvt_region)->get_pvt_input_table ();
      }
    else if (pvt_fluid_type == FI_PHASE_WATER)
      {
        return get_pvt_water (index_pvt_region)->get_pvt_input_table (); 
      }  
    else if (pvt_fluid_type == FI_PHASE_GAS)
      {  
        return get_pvt_gas (index_pvt_region)->get_pvt_input_table ();
      }     
	  else return NULL;
  }
  
  //! return input table for density for defined pvt region 
  BS_SP (table_iface)
  pvt_3p::get_density_table (t_long index_pvt_region) const 
  {
    BS_ASSERT (index_pvt_region >= 0 && index_pvt_region < n_pvt_regions);
    return density_table[index_pvt_region];
  }
  
  std::list <BS_SP( table_iface)>
  pvt_3p::get_density_list () const
  {
    std::list<BS_SP( table_iface)> tables;
    
    for (t_long i = 0; i < n_pvt_regions; ++i)
      tables.push_back (get_density_table (i));
    return tables;
  }
  
  
  void 
  pvt_3p::set_density_to_pvt_internal ()
  {
    t_float *density_data = &(*density)[0];
    for (t_long i = 0; i < n_pvt_regions; i++)
      {
         density_table[i]->set_value (0, 0, density_data[i * FI_PHASE_TOT + 0]);
         density_table[i]->set_value (0, 1, density_data[i * FI_PHASE_TOT + 1]);
         density_table[i]->set_value (0, 2, density_data[i * FI_PHASE_TOT + 2]);
      }
  }
  
  void
  pvt_3p::set_density_to_pvt_arrays ()
  {
    for (t_long i = 0; i < n_pvt_regions; i++)
      {
         if (pvt_oil_array[i].get ())
           pvt_oil_array[i]->set_surface_density (density_table[i]->get_value (0, 0));
         if (pvt_water_array[i].get ())
           pvt_water_array[i]->set_surface_density (density_table[i]->get_value (0, 1));     
         if (pvt_gas_array[i].get ())
           pvt_gas_array[i]->set_surface_density (density_table[i]->get_value (0, 2)); 
      }      
  }

	//! build pvt internal tables 
	void 
	pvt_3p::init_pvt_calc_data (t_float atm_p, t_float min_p, t_float max_p, t_float n_intervals)
	{
    // set density 
    set_density_to_pvt_arrays ();
    // build calc pvt data 
    for (t_long i = 0; i < n_pvt_regions; i++)
      {
        if (pvt_oil_array[i].get ())
          pvt_oil_array[i]->build (atm_p, min_p, max_p, n_intervals);

        if (pvt_gas_array[i].get ())
          pvt_gas_array[i]->build (atm_p, min_p, max_p, n_intervals);

        if (pvt_water_array[i].get ())
          pvt_water_array[i]->build (atm_p, min_p, max_p, n_intervals);
      }
	}           
   
  BLUE_SKY_TYPE_STD_CREATE (pvt_3p);
  BLUE_SKY_TYPE_STD_COPY (pvt_3p);
  BLUE_SKY_TYPE_IMPL(pvt_3p, pvt_3p_iface, "pvt_3p", "pvt_3p calculation class", "pvt_3p calculation");
  //////////////////////////////////////////////////////////////////////////

} // namespace blue_sky
