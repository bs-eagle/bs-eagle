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

namespace blue_sky 
{

  class BS_API_PLUGIN jfunction;
  
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
    add_spof (spv_float const &swof, bool is_water) = 0;

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
  };

  class BS_API_PLUGIN scal_3p_iface : public objbase
  {
  public:

    typedef boost::array <t_long, FI_PHASE_TOT>			phase_d_t;
    typedef boost::array <t_long, FI_PHASE_TOT>			sat_d_t;
    typedef std::vector <calc_model_data>           data_array_t;

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
    init (bool            is_w, 
      bool                is_g, 
      bool                is_o, 
      const phase_d_t &   phase_d, 
      const phase_d_t &   sat_d,
      RPO_MODEL_ENUM      rpo_model, 
      bool                is_scalecrs_ = false) = 0;

    virtual void
    set_water_jfunction (BS_SP (jfunction) jfunc) = 0;

    virtual void
    set_gas_jfunction (BS_SP (jfunction) jfunc) = 0;

    virtual void
    update_gas_data () = 0;

  };

}
#endif //

