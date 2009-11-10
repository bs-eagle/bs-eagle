/**
 *
 * */
#ifndef BS_CALC_MODEL_TYPE_HELPER_H_
#define BS_CALC_MODEL_TYPE_HELPER_H_

namespace blue_sky
  {

  template <typename strategy_type>
  struct calc_model_type_helper
    {
      typedef strategy_type                                             strategy_t;
      typedef typename strategy_t::item_t                               item_t;
      typedef typename strategy_t::index_t                              index_t;

      typedef boost::array <item_t, FI_PHASE_TOT>                       well_mobility_t;
      typedef boost::array <item_t, FI_PHASE_TOT>                       p_deriv_well_mobility_t;
      typedef boost::array <item_t, FI_PHASE_TOT * FI_PHASE_TOT>        sat_deriv_well_mobility_t;
      typedef boost::array <item_t, FI_PHASE_TOT>                       xref_deriv_well_mobility_t;
      typedef boost::array <item_t, FI_PHASE_TOT * (FI_PHASE_TOT - 1)>  flow_deriv_well_mobility_t;
      typedef boost::array <item_t, FI_PHASE_TOT>                       phase_fvf_t;
    };


} // namespace blue_sky


#endif  // #ifndef BS_CALC_MODEL_TYPE_HELPER_H_

