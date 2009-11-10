/**
 * \file well_rate_control_deriv_typedef.h
 * \brief
 * \author Sergey Miryanov
 * \date 24.11.2008
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

