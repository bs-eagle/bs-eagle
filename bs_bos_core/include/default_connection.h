/**
 * \file default_connection.h
 * \brief default impl of connection
 * \author Sergey Miryanov
 * \date 20.05.2009
 * */
#ifndef BS_BOS_CORE_DEFAULT_CONNECTION_H_
#define BS_BOS_CORE_DEFAULT_CONNECTION_H_

#include "well_connection.h"
#include "array_ext.h"

namespace blue_sky {
namespace wells {

  template <typename strategy_t>
  class BS_API_PLUGIN default_connection : public connection <strategy_t>
  {
  public:

    typedef connection <strategy_t>             base_t;
    typedef typename base_t::item_t             item_t;
    typedef typename base_t::rhs_item_t         rhs_item_t;

  public:

    BLUE_SKY_TYPE_DECL_T (default_connection <strategy_t>);

    void clear_data ();
    array_ext <item_t> get_rw_value ();
    array_ext <item_t> get_wr_value ();
    array_ext <item_t> get_rr_value ();
    array_ext <item_t> get_ps_value ();
    array_ext <rhs_item_t> get_rate_value ();

  public:

    enum 
      {
        rr_value_count = FI_PHASE_TOT * FI_PHASE_TOT + FI_PHASE_TOT,
      };

    boost::array <item_t, FI_PHASE_TOT>                 mobility_value;
    boost::array <rhs_item_t, FI_PHASE_TOT>             rate_value;
    boost::array <item_t, FI_PHASE_TOT>                 ps_value;
    boost::array <item_t, rr_value_count>               rr_value;
    boost::array <item_t, FI_PHASE_TOT>                 rw_value;
    boost::array <item_t, FI_PHASE_TOT>                 wr_value;
  };


} // namespace wells
} // namespace blue_sky

#endif // #ifndef BS_BOS_CORE_DEFAULT_CONNECTION_H_


