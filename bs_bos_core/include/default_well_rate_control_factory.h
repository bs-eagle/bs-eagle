/**
 * \file default_well_rate_control_factory.h
 * \brief factory to construct well_rate_control for default well model
 * \author Sergey Miryanov
 * \date 21.05.2009
 * */
#ifndef BS_BOS_CORE_DEFAULT_WELL_RATE_CONTROL_FACTORY_H_
#define BS_BOS_CORE_DEFAULT_WELL_RATE_CONTROL_FACTORY_H_

#include "well_rate_control_interface.h"

namespace blue_sky {
namespace wells {

  template <typename strategy_t>
  class BS_API_PLUGIN default_well_rate_control_factory : public well_rate_control_factory <strategy_t>
  {
  public:

    typedef well_rate_control_factory <strategy_t>    base_t;
    typedef typename base_t::sp_calc_model_t          sp_calc_model_t;
    typedef typename base_t::sp_well_rate_control_t   sp_well_rate_control_t;

    sp_well_rate_control_t
    create_control (rate_control_type control_type, bool is_bhp, bool is_production, const sp_calc_model_t &calc_model);

    BLUE_SKY_TYPE_DECL_T (default_well_rate_control_factory);
  };



} // namespace wells
} // namespace blue_sky

#endif  // #ifndef BS_BOS_CORE_DEFAULT_WELL_RATE_CONTROL_FACTORY_H_

