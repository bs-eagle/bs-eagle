/*
 * \file scal_3p.h
 * \brief scal 3p
 * \author Sergey Miryanov
 * \date 19.05.2008
 * */
#ifndef BS_SCAL_3P_H_
#define BS_SCAL_3P_H_

#include "scal_3p_iface.hpp"
#include "jfunction.h"
#include "scal_dummy_iface.h"

#include "bs_serialize_decl.h"

namespace blue_sky
{

  class BS_API_PLUGIN scal_region;

  class BS_API_PLUGIN scale_array_holder;

  class BS_API_PLUGIN scal_2p_data_holder;
  
  //class BS_API_PLUGIN jfunction;

  //////////////////////////////////////////////////////////////////////////
  class BS_API_PLUGIN scal_3p : public scal_3p_iface
    {
    public:
      typedef scal_3p             											this_t;

      typedef t_float                    							item_t;
      typedef t_long                       							index_t;
      typedef v_long                                    index_array_t;
      typedef v_float                                   item_array_t;
      typedef smart_ptr <index_array_t, true>           sp_array_index_t;
      typedef smart_ptr <item_array_t, true>            sp_array_item_t;

      typedef scal_region             									scal_region_t;
      typedef scale_array_holder            						scale_array_holder_t;
      typedef scal_2p_data_holder              					scal_2p_data_holder_t;

      
      typedef smart_ptr <scale_array_holder_t, true>		sp_scale_array_holder_t;
      typedef smart_ptr <scal_2p_data_holder_t, true>		sp_scal_2p_data_holder_t;
      typedef boost::array <index_t, FI_PHASE_TOT>			phase_d_t;
      typedef boost::array <index_t, FI_PHASE_TOT>			sat_d_t;

      typedef calc_model_data                               data_t;
      typedef std::vector <calc_model_data>                 data_array_t;

	  typedef smart_ptr <scal_dummy_iface, true>			        sp_scal_dummy_iface;
      typedef unsigned char															    phase_index_t;
      typedef unsigned char															    sat_index_t;

      typedef scal_3p_iface::sp_scal_input_table_t          sp_scal_input_table_t;
      typedef scal_3p_iface::sp_scal_input_table_array_t    sp_scal_input_table_array_t; 

      struct scal_3p_impl_base;

    public:

      ~scal_3p ();

      BS_SP (scale_array_holder_iface)   
      get_water_scale () const
      {
        return water_scale;
      }
      BS_SP (scale_array_holder_iface)
      get_gas_scale () const
      {
        return gas_scale;
      }

      BS_SP (scal_2p_data_holder_iface)
      get_water_data () const
      {
        return water_data;
      }
      BS_SP (scal_2p_data_holder_iface)
      get_gas_data () const
      {
        return gas_data;
      }

      sp_jfunction_t
      get_water_jfunction () const
      {
        return water_jfunc;
      }
      sp_jfunction_t
      get_gas_jfunction () const
      {
        return gas_jfunc;
      }

      void
      get_relative_perm (index_t cell_index, 
        const sp_array_item_t saturation, 
        const sp_array_index_t sat_regions, 
        sp_array_item_t relative_perm, 
        sp_array_item_t s_deriv_relative_perm) const;

        void
        get_capillary (index_t cell_index, 
          const sp_array_item_t saturation, 
          const sp_array_index_t sat_regions, 
          const sp_array_item_t perm, 
          const sp_array_item_t poro, 
          sp_array_item_t cap, 
          sp_array_item_t s_deriv_cap) const;

      void
      process (const spv_double& saturation, 
        const spv_long &sat_regions,
        const stdv_float &perm,
        const stdv_float &poro,
        data_array_t &data) const;

      void
      process_init (index_t cell_index, 
        const item_t *pressure, 
        index_t sat_reg, 
        const item_t *perm_array, 
        item_t poro,
        item_t *sat, 
        item_t *pc_limit) const;

      void
      process_init_2 (
        const item_t *pressure, 
        index_t sat_reg, 
        item_t perm, 
        item_t poro,
        item_t *sat, 
        item_t *pc_limit) const;

      void
      calc_pcp (index_t cell_index, 
        const item_t sat, 
        index_t sat_reg, 
        item_t cap, 
        item_t &pcp) const;

      void
      calc_gas_water_zone (index_t cell_index, 
        index_t sat_reg, 
        const item_t *perm_array, 
        item_t poro, 
        item_t pcgw,
        item_t &sw, 
        item_t &sg) const;

      void
      calc_gas_water_zone_2 (
        index_t sat_reg, 
        item_t perm, 
        item_t poro, 
        item_t pcgw,
        item_t &sw, 
        item_t &sg) const;

      void
      init (bool is_w, bool is_g, bool is_o, 
        const phase_d_t &phase_d, const phase_d_t &sat_d,
        RPO_MODEL_ENUM rpo_model, bool is_scalecrs_ = false);

      void
      set_water_jfunction (sp_jfunction_t jfunc);

      void
      set_gas_jfunction (sp_jfunction_t jfunc);

      t_long 
      get_n_scal_regions ()
        { return n_scal_regions; }

      void
      update_gas_data ();

	    void
	    init_scal_calc_data ();

	    void
	    init_from_scal_ex(const phase_d_t &phase_d, const sat_d_t &sat_d,
					 sp_jfunction_t water_jfunc,
					 sp_jfunction_t gas_jfunc);
	    
	    //! init of input scal tables 
      void
      init_scal_input_table_arrays (const t_long n_scal_regions_, 
                                    bool is_oil_, bool is_gas_, bool is_water_);
      //! return input scal table for defined region and fluid type
      virtual BS_SP (table_iface)
      get_table (t_int scal_fluid_type, t_long index_scal_region) const;
      
      std::list <BS_SP (table_iface)> 
      get_tables_list (t_long index_scal_region) const;

      std::list <BS_SP (table_iface)>
      get_tables_fluid_all_regions (t_long scal_fluid_type) const;
      
      virtual void 
      init_scal_data_from_input_tables ();
      
    private:

      sp_scal_2p_data_holder_t  water_data;
      sp_scal_2p_data_holder_t  gas_data;

      sp_scale_array_holder_t   water_scale;
      sp_scale_array_holder_t   gas_scale;

      sp_jfunction_t            water_jfunc;
      sp_jfunction_t            gas_jfunc;

      scal_3p_impl_base         *impl_;
      
      // INPUT TABLES DATA
      t_long                       n_scal_regions;
      sp_scal_input_table_array_t  water_input_table;
      sp_scal_input_table_array_t  gas_input_table;
      sp_scal_input_table_array_t  oil_input_table;
      
	  bool is_gas, is_oil, is_water;

      friend class blue_sky::bs_serialize;

    public:

      BLUE_SKY_TYPE_DECL (scal_3p);
    };

  bool scal_register_types (const blue_sky::plugin_descriptor &pd);

} // namespace blue_sky

#endif	// #ifndef BS_SCAL_3P_H_
