/**
 *       \file  well_type_helper.h
 *      \brief  Type helper, stores common well types
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  09.09.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Obsolete, should be removed
 * */
#ifndef BS_WELL_TYPE_HELPER_H_
#define BS_WELL_TYPE_HELPER_H_

#include "conf.h"
#include BS_FORCE_PLUGIN_IMPORT ()
#include "constants.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {
  namespace wells
    {

    struct type_helper
      {
        typedef t_double               item_t;
        typedef t_long              index_t;

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

