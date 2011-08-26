#ifndef SCAL_3P_IFACE_H_ada7b1b4_70c4_11e0_9278_47db4902b249
#define SCAL_3P_IFACE_H_ada7b1b4_70c4_11e0_9278_47db4902b249
/**
 *       \file  scal_3p_iface.hpp
 *      \brief  Interface for scal_3p
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  27.04.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "calc_model_data.h"
#include "table_iface.h"

namespace blue_sky 
{

  class BS_API_PLUGIN jfunction;
  class BS_API_PLUGIN scal_input_table;
  
  enum 
  {
    SOF2_KEYWORD_COLUMNS = 2,
    SOF3_KEYWORD_COLUMNS = 3,
    SPFN_KEYWORD_COLUMNS = 3,
    SPOF_KEYWORD_COLUMNS = 4
  };

  enum scale_array_name
  {
    socr,
    scr,
    su,
    sl,
    pcp,
    krp,
    krop,
    krpr,
    krorp,
    scale_array_name_total
  };

  class BS_API_PLUGIN scale_array_holder_iface : public objbase
  {
  public:
    virtual ~scale_array_holder_iface () {}

    virtual void remove (scale_array_name array) = 0;

    virtual t_double get_sl (t_long cell, t_double value) const = 0;
    virtual t_double get_su (t_long cell, t_double value) const = 0;

    virtual void set (scale_array_name array, std::string const &name, spv_float const &data) = 0;
  };

  class BS_API_PLUGIN scal_2p_data_holder_iface : public objbase
  {
  public:
    virtual ~scal_2p_data_holder_iface () {}

    virtual void
    add_spof (spv_float const &swof, t_long index, bool is_water) = 0;

    virtual void
    add_spfn (spv_float const &swfn, t_long index, bool is_water) = 0;

    virtual void
    add_sof3 (spv_float const &sof3, t_long index, bool is_water) = 0;

    virtual void
    add_sof2 (spv_float const &sof3, t_long index, bool is_water) = 0;
    
    virtual t_float
    get_phase_sat_min (t_long region) const = 0;

    virtual t_float
    get_phase_sat_max (t_long region) const = 0;

    virtual t_float
    get_pcp_max (t_long region) const = 0;
    
    virtual void 
    init_table_array (const t_long n_scal_regions) = 0;
    
    virtual void
    init_regions_from_tables () = 0;
    
	virtual void
	clear_regions() = 0;
  };

  class BS_API_PLUGIN scal_3p_iface : public objbase
  {
  public:

    typedef boost::array <t_long, FI_PHASE_TOT>			phase_d_t;
    typedef boost::array <t_long, FI_PHASE_TOT>			sat_d_t;
    typedef std::vector <calc_model_data>           data_array_t;
    
    typedef smart_ptr <table_iface, true>           sp_scal_input_table_t;
    typedef std::vector <sp_scal_input_table_t>     sp_scal_input_table_array_t;    

	typedef jfunction             					jfunction_t;
    typedef smart_ptr <jfunction_t, true>			sp_jfunction_t;

    virtual ~scal_3p_iface () {}

    virtual BS_SP (scale_array_holder_iface)
    get_water_scale () const = 0;

    virtual BS_SP (scale_array_holder_iface) 
    get_gas_scale () const = 0;

    virtual BS_SP (scal_2p_data_holder_iface)
    get_water_data () const = 0;

    virtual BS_SP (scal_2p_data_holder_iface)
    get_gas_data () const = 0;

    virtual BS_SP (jfunction)
    get_water_jfunction () const = 0;

    virtual BS_SP (jfunction)
    get_gas_jfunction () const = 0;

    virtual void
    process (const spv_double & saturation, 
      const spv_long &          sat_regions,
      const stdv_float &        perm,
      const stdv_float &        poro,
      data_array_t &            data) const = 0;


    virtual void
    process_init (t_long  cell_index, 
      const t_double *    pressure, 
      t_long              sat_reg, 
      const t_double *    perm_array, 
      t_double            poro,
      t_double *          sat, 
      t_double *          pc_limit) const = 0;

    virtual void
    process_init_2 (
      const t_double *    pressure, 
      t_long              sat_reg, 
      t_double            perm, 
      t_double            poro,
      t_double *          sat, 
      t_double *          pc_limit) const = 0;


    virtual void
    calc_pcp (t_long      cell_index, 
      const t_double      sat, 
      t_long              sat_reg, 
      t_double            cap, 
      t_double &          pcp) const = 0;

    virtual void
    calc_gas_water_zone (t_long cell_index, 
      t_long              sat_reg, 
      const t_double *    perm_array, 
      t_double            poro, 
      t_double            pcgw,
      t_double &          sw, 
      t_double &          sg) const = 0;

    virtual void
    calc_gas_water_zone_2 (
      t_long              sat_reg, 
      t_double            perm, 
      t_double            poro, 
      t_double            pcgw,
      t_double &          sw, 
      t_double &          sg) const = 0;


    virtual void
    init (bool            is_w, 
      bool                is_g, 
      bool                is_o, 
      const phase_d_t &   phase_d, 
      const phase_d_t &   sat_d,
      RPO_MODEL_ENUM      rpo_model, 
      bool                is_scalecrs_ = false) = 0;
    
    virtual t_long 
    get_n_scal_regions () = 0;
    
    virtual void
    init_from_scal () = 0;

    virtual void
    init_from_scal_ex(const phase_d_t &phase_d, const sat_d_t &sat_d,
		                  sp_jfunction_t water_jfunc,
                      sp_jfunction_t gas_jfunc) = 0;
  
    virtual BS_SP (table_iface)
    get_table (t_int scal_fluid_type, t_long index_scal_region) const = 0;
    
    virtual std::list <BS_SP (table_iface)> 
    get_tables_list (t_long index_scal_region) const = 0;
    
    virtual std::list <BS_SP (table_iface)>
    get_tables_fluid_all_regions (t_long scal_fluid_type) const = 0;
    
    virtual void
    init_scal_input_table_arrays (const t_long    n_scal_regions_, 
                                bool            is_oil, 
                                bool            is_gas, 
                                bool            is_water) = 0;

    virtual void
    set_water_jfunction (BS_SP (jfunction) jfunc) = 0;

    virtual void
    set_gas_jfunction (BS_SP (jfunction) jfunc) = 0;

    virtual void
    update_gas_data () = 0;

  };

}
#endif //

