/**
 *       \file  well_rate_control_compute_deriv_typedef.h
 *      \brief  Get types from base class (well_rate_control_deriv)
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  24.11.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Obsolete, should be removed
 * */
#ifndef BS_WELLS_WELL_RATE_CONTROL_DERIV_TYPEDEF_H_
#define BS_WELLS_WELL_RATE_CONTROL_DERIV_TYPEDEF_H_

#define GET_COMPUTE_DERIV_BASE_TYPES                                            \
  typedef compute_deriv <mobility_calc_t>                 base_t;               \
  typedef typename base_t::data_t                         data_t;               \
  typedef typename base_t::params_t                       params_t;             \
  typedef typename base_t::item_t                         item_t;               \
  typedef typename base_t::rhs_item_t                     rhs_item_t;           \
  typedef typename base_t::item_q_rate_t                  item_q_rate_t;        \
  typedef typename base_t::item_rhs_block_t               item_rhs_block_t;     \
  typedef typename base_t::item_rr_block_t                item_rr_block_t;      \
  typedef typename base_t::item_rw_block_t                item_rw_block_t;      \
  typedef typename base_t::item_wr_block_t                item_wr_block_t;      \
  typedef typename base_t::item_ps_block_t                item_ps_block_t;      \
  typedef typename base_t::sp_connection_t                sp_connection_t;



#endif // #ifndef BS_WELLS_WELL_RATE_CONTROL_DERIV_TYPEDEF_H_

