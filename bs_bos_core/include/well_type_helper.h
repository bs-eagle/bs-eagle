/**
 *
 * */
#ifndef BS_WELL_TYPE_HELPER_H_
#define BS_WELL_TYPE_HELPER_H_

#include BS_FORCE_PLUGIN_IMPORT ()
#include "constants.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {
  namespace wells
    {

    template <typename strategy_type>
    struct type_helper
      {
        typedef strategy_type                             strategy_t;
        typedef typename strategy_t::item_t               item_t;
        typedef typename strategy_t::index_t							index_t;

        typedef boost::array <index_t, FI_PHASE_TOT>      phase_d_t;
        typedef boost::array <index_t, FI_PHASE_TOT>      sat_d_t;

        typedef item_t                                    item_ww_block_t;
        typedef boost::array <item_t, FI_PHASE_TOT * FI_PHASE_TOT>
        item_rr_block_t;

        typedef boost::array <item_t, FI_PHASE_TOT>       item_wr_block_t;
        typedef boost::array <item_t, FI_PHASE_TOT>       item_rw_block_t;
        typedef boost::array <item_t, FI_PHASE_TOT>       item_rhs_block_t;
        typedef boost::array <item_t, FI_PHASE_TOT>       item_ps_block_t;
        typedef boost::array <item_t, FI_PHASE_TOT>       item_sp_block_t;

        typedef boost::array <item_t, FI_PHASE_TOT>       item_q_rate_t;
        typedef boost::array <item_t, FI_PHASE_TOT>       item_q_rate_inflow_t;
        typedef boost::array <item_t, GAS_RATE_TOTAL>     item_gas_rate_t;

        typedef boost::array <item_t, FI_PHASE_TOT>       invers_fvf_avgerage_t;

        typedef boost::array <item_t, FI_PHASE_TOT>       initial_rate_t;
      };


  } // namespace wells
} // namespace blue_sky

#endif  // #ifndef BS_WELL_TYPE_HELPER_H_

