/// @file scal_3p_impl_base.h
/// @brief Interface for scal_3p implementation class
/// @author uentity
/// @version 
/// @date 27.01.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef SCAL_3P_IMPL_BASE_XQZHMFQL
#define SCAL_3P_IMPL_BASE_XQZHMFQL

#include "conf.h"
#include "scal_3p.h"
#include <string>

namespace blue_sky {

  //////////////////////////////////////////////////////////////////////////
  struct scal_3p::scal_3p_impl_base
  {
    typedef t_float                           item_t;
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
    process_init_2 (const item_t *pressure, index_t sat_reg, item_t perm, item_t poro, item_t *sat, item_t *pc_limit) const = 0;

    virtual void
    calc_pcp (index_t cell_index, item_t sat, index_t sat_reg, item_t cap, item_t &pcp) const = 0;

    virtual void
    calc_gas_water_zone (index_t cell_index, index_t sat_reg, const item_t *perm_array, item_t poro, item_t pcgw, item_t &sw, item_t &sg) const = 0;

    virtual void
    calc_gas_water_zone_2 (index_t sat_reg, const item_t perm, const item_t poro, item_t pcgw, item_t &sw, item_t &sg) const = 0;

    virtual bool
    is_water () const = 0;

    virtual bool
    is_gas () const = 0;

    virtual bool
    is_oil () const = 0;

    virtual int
    get_rpo_model () const = 0;

    virtual int
    get_n_phases () const = 0;

    virtual std::string
    dump_state() const = 0;

    virtual void
    restore_state(const std::string&) = 0;
  };

} /* blue_sky */

#endif /* end of include guard: SCAL_3P_IMPL_BASE_XQZHMFQL */

