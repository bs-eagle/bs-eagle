/**
 *       \file  compute_acc_rates.h
 *      \brief  Calculates total rates
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  22.12.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_MAIN_LOOP_CALC_COMPUTE_ACC_RATES_H_
#define BS_MAIN_LOOP_CALC_COMPUTE_ACC_RATES_H_

#include "apply_wefac.h"

namespace blue_sky
  {

  /**
   * \brief  Calculates total rates for well and 
   *         rate and total rate for reservoir
   * */
  template <bool is_w, bool is_g, bool is_o>
  inline void
  main_loop_calc <is_w, is_g, is_o>::compute_acc_rates ()
  {
    typedef facility_manager::well_const_iterator_t well_iterator_t;
    typedef smart_ptr <well> sp_well_t;

    well_iterator_t wb = reservoir_->get_facility_list ()->wells_begin ();
    well_iterator_t we = reservoir_->get_facility_list ()->wells_end ();

    reservoir_t::rate_data_t &rs_rate          = reservoir_->rate_;
    reservoir_t::rate_data_t &rs_rate_rc       = reservoir_->rate_rc_;
    reservoir_t::rate_data_t &rs_rate_wefac    = reservoir_->rate_wefac_;
    reservoir_t::rate_data_t &rs_rate_rc_wefac = reservoir_->rate_rc_wefac_;
    reservoir_t::rate_data_t &rs_rate_initial  = reservoir_->rate_initial_;
    reservoir_t::rate_data_t &rs_rate_total    = reservoir_->rate_total_;

    rs_rate           = 0;
    rs_rate_rc        = 0;
    rs_rate_wefac     = 0;
    rs_rate_rc_wefac  = 0;
    rs_rate_initial   = 0;
    rs_rate_total     = 0;

    for (; wb != we; ++wb)
      {
        sp_well_t well (wb->second, bs_dynamic_cast ());

        bool is_prod    = well->get_well_controller ()->is_production ();
        item_t wefac    = apply_wefac (1.0, well->exploitation_factor_);
        item_t wefac_dt = apply_wefac (dt_, wefac);

        //if (well->is_shut () || !well->well_state_.is_work)
        //  continue;

        // TODO: OPT:
        if (is_prod)
          {
            well->rate_total_.prod += well->rate ().prod * wefac_dt;
          }
        else
          {
            well->rate_total_.inj += well->rate ().inj * wefac_dt;
          }

        rs_rate           += well->rate ();
        rs_rate_rc        += well->rate_rc_;
        rs_rate_wefac     += well->rate () * wefac;
        rs_rate_rc_wefac  += well->rate_rc_ * wefac;
        rs_rate_initial   += well->rate () * wefac;
        rs_rate_total     += well->rate_total_;
      }
  }


} // namespace blue_sky



#endif  // #ifndef BS_MAIN_LOOP_CALC_COMPUTE_ACC_RATES_H_
