/**
 * \file default_connection.h
 * \brief impl of
 * \author Sergey Miryanov
 * \date 20.05.2009
 * */
#include "stdafx.h"
#include "default_connection.h"

namespace blue_sky {
namespace wells {

  template <typename strategy_t>
  default_connection <strategy_t>::default_connection (bs_type_ctor_param param /* = NULL */)
  : base_t (param)
  {
    clear_data ();
  }
  template <typename strategy_t>
  default_connection <strategy_t>::default_connection (const default_connection &c)
        : bs_refcounter (), connection <strategy_t> () 
  {
    *this = c;
  }

  template <typename strategy_t>
  array_ext <typename default_connection <strategy_t>::item_t>
  default_connection <strategy_t>::get_rw_value ()
  {
    return array_ext <item_t> (&rw_value[0], rw_value.size ());
  }
  template <typename strategy_t>
  array_ext <typename default_connection <strategy_t>::item_t>
  default_connection <strategy_t>::get_wr_value ()
  {
    return array_ext <item_t> (&wr_value[0], wr_value.size ());
  }
  template <typename strategy_t>
  array_ext <typename default_connection <strategy_t>::item_t>
  default_connection <strategy_t>::get_rr_value ()
  {
    return array_ext <item_t> (&rr_value[0], rr_value.size ());
  }
  template <typename strategy_t>
  array_ext <typename default_connection <strategy_t>::item_t>
  default_connection <strategy_t>::get_ps_value ()
  {
    return array_ext <item_t> (&ps_value[0], ps_value.size ());
  }
  template <typename strategy_t>
  array_ext <typename default_connection <strategy_t>::rhs_item_t>
  default_connection <strategy_t>::get_rate_value ()
  {
    return array_ext <rhs_item_t> (&rate_value[0], rate_value.size ());
  }

  template <typename strategy_t>
  void
  default_connection <strategy_t>::clear_data ()
  {
    assign (mobility_value, 0);
    assign (rate_value, 0);
    assign (ps_value, 0);
    assign (rr_value, 0);
    assign (rw_value, 0);
    assign (wr_value, 0);

    base_t::clear_data ();
  }

  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE_T_DEF (default_connection, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF (default_connection, (class));
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (default_connection<base_strategy_fi>), 1, (connection<base_strategy_fi>), "default_connection_seq_fi", "default_connection_seq_fi", "default_connection_seq_fi", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (default_connection<base_strategy_di>), 1, (connection<base_strategy_di>), "default_connection_seq_di", "default_connection_seq_di", "default_connection_seq_di", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (default_connection<base_strategy_mixi>), 1, (connection<base_strategy_mixi>), "default_connection_seq_mixi", "default_connection_seq_mixi", "default_connection_seq_mixi", false);
  //////////////////////////////////////////////////////////////////////////

} //namespace wells
} //namespace blue_sky

