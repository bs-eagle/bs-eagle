/**
 *       \file  well_rate_control_rate_deriv.h
 *      \brief  Calculates derivs for well controled by rate
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  24.11.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Obsolete, should be removed
 * */
#ifndef BS_WELLS_WELL_RATE_CONTROL_RATE_DERIV_H_
#define BS_WELLS_WELL_RATE_CONTROL_RATE_DERIV_H_

#include "well_rate_control_deriv.h"
#include "well_rate_control_compute_deriv_typedef.h"

namespace blue_sky
  {
  namespace wells
    {

    template <typename mobility_calc_t>
    struct compute_rate_deriv_3p : public compute_deriv <mobility_calc_t>
      {
        GET_COMPUTE_DERIV_BASE_TYPES;

        compute_rate_deriv_3p (const mobility_calc_t &mobility_calc)
            : base_t (mobility_calc)
        {
        }

        void
        oil_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            shared_vector <item_t> rw_block = c->get_rw_value ();
            shared_vector <item_t> wr_block = c->get_wr_value ();
            shared_vector <item_t> rr_block = c->get_rr_value ();
            shared_vector <item_t> ww_value = params.ww_value;
            item_t rw                   = 0;

            if (mobility_calc_t::is_oil_injection (params))
              {
                rw = rw_block [p3_oil] = /*this->mult * */base_t::compute_oil_pref_deriv (data, params);
              }
            if (mobility_calc_t::is_o_ctrl (params))
              {
                ww_value[0] += rw/* * this->mult*/;
                wr_block[0] += (rr_block[p3_oil * 3 + 0]/* * this->mult*/);
                wr_block[1] += (rr_block[p3_oil * 3 + 1]/* * this->mult*/);
                wr_block[2] += (rr_block[p3_oil * 3 + 2]/* * this->mult*/);
              }
          }
        void
        water_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            shared_vector <item_t> rw_block = c->get_rw_value ();
            shared_vector <item_t> wr_block = c->get_wr_value ();
            shared_vector <item_t> rr_block = c->get_rr_value ();
            shared_vector <item_t> ww_value = params.ww_value;
            item_t rw                   = 0;

            if (mobility_calc_t::is_water_injection (params))
              {
                rw = rw_block [p3_wat] = /*this->mult * */base_t::compute_water_pref_deriv (data, params);
#ifdef _DEBUG
                //BOSOUT (section::wells, level::debug) << boost::format ("[%s : %d] rw: %.20e") % params.well_->name ().c_str () % params.n_block % rw << bs_end;
#endif
              }
            if (mobility_calc_t::is_w_ctrl (params))
              {
                ww_value[0] += rw/* * this->mult*/;
                wr_block[0] += (rr_block[p3_wat * 3 + 0]/* * this->mult*/);
                wr_block[1] += (rr_block[p3_wat * 3 + 1]/* * this->mult*/);
                wr_block[2] += (rr_block[p3_wat * 3 + 2]/* * this->mult*/);
              }
          }
        void
        gas_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            shared_vector <item_t> rw_block = c->get_rw_value ();
            shared_vector <item_t> wr_block = c->get_wr_value ();
            shared_vector <item_t> rr_block = c->get_rr_value ();
            shared_vector <item_t> ww_value = params.ww_value;
            item_t rw                   = 0;

            if (mobility_calc_t::is_gas_injection (params))
              {
                rw = rw_block [p3_gas] = /*this->mult * */base_t::compute_gas_pref_deriv (data, params);
              }
            if (mobility_calc_t::is_g_ctrl (params))
              {
                ww_value[0] += rw/* * this->mult*/;
                wr_block[0] += (rr_block[p3_gas * 3 + 0]/* * this->mult*/);
                wr_block[1] += (rr_block[p3_gas * 3 + 1]/* * this->mult*/);
                wr_block[2] += (rr_block[p3_gas * 3 + 2]/* * this->mult*/);
              }
          }
        void
        update_wr (const sp_connection_t &c, item_t ww) const
          {
            c->get_wr_value ()[0] *= ww;
            c->get_wr_value ()[1] *= ww;
            c->get_wr_value ()[2] *= ww;
          }
        void
        apply_wefac (const sp_connection_t &c, params_t &params) const
        {
          shared_vector <item_t> rw = c->get_rw_value ();
          item_t wefac = params.well_->exploitation_factor_;

          rw[0] = blue_sky::apply_wefac (rw[0], wefac);
          rw[1] = blue_sky::apply_wefac (rw[1], wefac);
          rw[2] = blue_sky::apply_wefac (rw[2], wefac);
        }

        void
        update_rate (const sp_connection_t &c, params_t &params) const
        {
          item_t ww_bw = (params.bw_value[0] / params.ww_value[0]);

          c->get_rate_value ()[0] -= c->get_rw_value () [0] * ww_bw;
          c->get_rate_value ()[1] -= c->get_rw_value () [1] * ww_bw;
          c->get_rate_value ()[2] -= c->get_rw_value () [2] * ww_bw;
        }
      };

    template <typename mobility_calc_t>
    struct compute_rate_deriv_2p_ow : public compute_deriv <mobility_calc_t>
      {
        GET_COMPUTE_DERIV_BASE_TYPES;

        compute_rate_deriv_2p_ow (const mobility_calc_t &mobility_calc)
            : base_t (mobility_calc)
        {
        }

        void
        oil_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            shared_vector <item_t> rw_block = c->get_rw_value ();
            shared_vector <item_t> wr_block = c->get_wr_value ();
            shared_vector <item_t> rr_block = c->get_rr_value ();
            shared_vector <item_t> ww_value = params.ww_value;
            item_t rw                   = 0;

            if (mobility_calc_t::is_oil_injection (params))
              {
                rw = rw_block [p2ow_oil] = /*this->mult * */base_t::compute_oil_pref_deriv (data, params);
              }
            if (mobility_calc_t::is_o_ctrl (params))
              {
                ww_value[0] += rw/* * this->mult*/;
                wr_block[0] += (rr_block[p2ow_oil * 2 + 0]/* * this->mult*/);
                wr_block[1] += (rr_block[p2ow_oil * 2 + 1]/* * this->mult*/);
              }
          }
        void
        water_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            shared_vector <item_t> rw_block = c->get_rw_value ();
            shared_vector <item_t> wr_block = c->get_wr_value ();
            shared_vector <item_t> rr_block = c->get_rr_value ();
            shared_vector <item_t> ww_value = params.ww_value;
            item_t rw                   = 0;

            if (mobility_calc_t::is_water_injection (params))
              {
                rw = rw_block [p2ow_wat] = /*this->mult * */base_t::compute_water_pref_deriv (data, params);
              }
            if (mobility_calc_t::is_w_ctrl (params))
              {
                ww_value[0] += rw/* * this->mult*/;
                wr_block[0] += (rr_block[p2ow_wat * 2 + 0]/* * this->mult*/);
                wr_block[1] += (rr_block[p2ow_wat * 2 + 1]/* * this->mult*/);
              }
          }
        void
        gas_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
          }
        void
        update_wr (const sp_connection_t &c, item_t ww) const
          {
            c->get_wr_value ()[0] *= ww;
            c->get_wr_value ()[1] *= ww;
          }
        void
        apply_wefac (const sp_connection_t &c, params_t &params) const
        {
          shared_vector <item_t> rw = c->get_rw_value ();
          item_t wefac = params.well_->exploitation_factor_;

          rw[0] = blue_sky::apply_wefac (rw[0], wefac);
          rw[1] = blue_sky::apply_wefac (rw[1], wefac);
        }
        void
        update_rate (const sp_connection_t &c, params_t &params) const
        {
          item_t ww_bw = (params.bw_value[0] / params.ww_value[0]);

          c->get_rate_value ()[0] -= c->get_rw_value () [0] * ww_bw;
          c->get_rate_value ()[1] -= c->get_rw_value () [1] * ww_bw;
        }
      };

    template <typename mobility_calc_t>
    struct compute_rate_deriv_2p_og : public compute_deriv <mobility_calc_t>
      {
        GET_COMPUTE_DERIV_BASE_TYPES;

        compute_rate_deriv_2p_og (const mobility_calc_t &mobility_calc)
            : base_t (mobility_calc)
        {
        }

        void
        oil_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            shared_vector <item_t> rw_block = c->get_rw_value ();
            shared_vector <item_t> wr_block = c->get_wr_value ();
            shared_vector <item_t> rr_block = c->get_rr_value ();
            shared_vector <item_t> ww_value = params.ww_value;
            item_t rw                   = 0;

            if (mobility_calc_t::is_oil_injection (params))
              {
                rw = rw_block [p2og_oil] = /*this->mult * */base_t::compute_oil_pref_deriv (data, params);
              }
            if (mobility_calc_t::is_o_ctrl (params))
              {
                ww_value[0] += rw/* * this->mult*/;
                wr_block[0] += (rr_block[p2og_oil * 2 + 0]/* * this->mult*/);
                wr_block[1] += (rr_block[p2og_oil * 2 + 1]/* * this->mult*/);
              }
          }
        void
        water_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
          }
        void
        gas_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            shared_vector <item_t> rw_block = c->get_rw_value ();
            shared_vector <item_t> wr_block = c->get_wr_value ();
            shared_vector <item_t> rr_block = c->get_rr_value ();
            shared_vector <item_t> ww_value = params.ww_value;
            item_t rw                   = 0;

            if (mobility_calc_t::is_gas_injection (params))
              {
                rw = rw_block [p2og_gas] = /*this->mult * */base_t::compute_gas_pref_deriv (data, params);
              }
            if (mobility_calc_t::is_g_ctrl (params))
              {
                ww_value[0] += rw/* * this->mult*/;
                wr_block[0] += (rr_block[p2og_gas * 2 + 0]/* * this->mult*/);
                wr_block[1] += (rr_block[p2og_gas * 2 + 1]/* * this->mult*/);
              }
          }
        void
        update_wr (const sp_connection_t &c, item_t ww) const
          {
            c->get_wr_value ()[0] *= ww;
            c->get_wr_value ()[1] *= ww;
          }
        void
        apply_wefac (const sp_connection_t &c, params_t &params) const
        {
          shared_vector <item_t> rw = c->get_rw_value ();
          item_t wefac = params.well_->exploitation_factor_;

          rw[0] = blue_sky::apply_wefac (rw[0], wefac);
          rw[1] = blue_sky::apply_wefac (rw[1], wefac);
        }
        void
        update_rate (const sp_connection_t &c, params_t &params) const
        {
          item_t ww_bw = (params.bw_value[0] / params.ww_value[0]);

          c->get_rate_value ()[0] -= c->get_rw_value () [0] * ww_bw;
          c->get_rate_value ()[1] -= c->get_rw_value () [1] * ww_bw;
        }
      };

    template <typename mobility_calc_t>
    struct compute_rate_deriv_1p_o : public compute_deriv <mobility_calc_t>
      {
        GET_COMPUTE_DERIV_BASE_TYPES;

        compute_rate_deriv_1p_o (const mobility_calc_t &mobility_calc)
            : base_t (mobility_calc)
        {
        }

        void
        oil_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            shared_vector <item_t> rw_block = c->get_rw_value ();
            shared_vector <item_t> wr_block = c->get_wr_value ();
            shared_vector <item_t> rr_block = c->get_rr_value ();
            shared_vector <item_t> ww_value = params.ww_value;
            item_t rw                   = 0;

            if (mobility_calc_t::is_oil_injection (params))
              {
                rw = rw_block [0] = /*this->mult * */base_t::compute_oil_pref_deriv (data, params);
              }
            if (mobility_calc_t::is_o_ctrl (params))
              {
                ww_value[0] += rw/* * this->mult*/;
                wr_block[0] += (rr_block[0]/* * this->mult*/);
              }
          }
        void
        water_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
          }
        void
        gas_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
          }
        void
        update_wr (const sp_connection_t &c, item_t ww) const
          {
            c->get_wr_value ()[0] *= ww;
          }
        void
        apply_wefac (const sp_connection_t &c, params_t &params) const
        {
          shared_vector <item_t> rw = c->get_rw_value ();
          item_t wefac = params.well_->exploitation_factor_;

          rw[0] = blue_sky::apply_wefac (rw[0], wefac);
        }
        void
        update_rate (const sp_connection_t &c, params_t &params) const
        {
          item_t ww_bw = (params.bw_value[0] / params.ww_value[0]);

          c->get_rate_value ()[0] -= c->get_rw_value () [0] * ww_bw;
        }
      };
    template <typename mobility_calc_t>
    struct compute_rate_deriv_1p_w : public compute_deriv <mobility_calc_t>
      {
        GET_COMPUTE_DERIV_BASE_TYPES;

        compute_rate_deriv_1p_w (const mobility_calc_t &mobility_calc)
            : base_t (mobility_calc)
        {
        }

        void
        oil_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
          }
        void
        water_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            shared_vector <item_t> rw_block = c->get_rw_value ();
            shared_vector <item_t> wr_block = c->get_wr_value ();
            shared_vector <item_t> rr_block = c->get_rr_value ();
            shared_vector <item_t> ww_value = params.ww_value;
            item_t rw                   = 0;

            if (mobility_calc_t::is_water_injection (params))
              {
                rw = rw_block [0] = /*this->mult * */base_t::compute_water_pref_deriv (data, params);
              }
            if (mobility_calc_t::is_w_ctrl (params))
              {
                ww_value[0] += rw/* * this->mult*/;
                wr_block[0] += (rr_block[0]/* * this->mult*/);
              }
          }
        void
        gas_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
          }
        void
        update_wr (const sp_connection_t &c, item_t ww) const
          {
            c->get_wr_value ()[0] *= ww;
          }
        void
        apply_wefac (const sp_connection_t &c, params_t &params) const
        {
          shared_vector <item_t> rw = c->get_rw_value ();
          item_t wefac = params.well_->exploitation_factor_;

          rw[0] = blue_sky::apply_wefac (rw[0], wefac);
        }
        void
        update_rate (const sp_connection_t &c, params_t &params) const
        {
          item_t ww_bw = (params.bw_value[0] / params.ww_value[0]);

          c->get_rate_value ()[0] -= c->get_rw_value () [0] * ww_bw;
        }
      };
    template <typename mobility_calc_t>
    struct compute_rate_deriv_1p_g : public compute_deriv <mobility_calc_t>
    {
      GET_COMPUTE_DERIV_BASE_TYPES;

      compute_rate_deriv_1p_g (const mobility_calc_t &mobility_calc)
          : base_t (mobility_calc)
      {
      }

      void
      oil_function (const sp_connection_t &c, const data_t &data, params_t &params) const
      {
      }
      void
      water_function (const sp_connection_t &c, const data_t &data, params_t &params) const
      {
      }
      void
      gas_function (const sp_connection_t &c, const data_t &data, params_t &params) const
      {
        shared_vector <item_t> rw_block = c->get_rw_value ();
        shared_vector <item_t> wr_block = c->get_wr_value ();
        shared_vector <item_t> rr_block = c->get_rr_value ();
        shared_vector <item_t> ww_value = params.ww_value;
        item_t rw                   = 0;

        if (mobility_calc_t::is_gas_injection (params))
          {
            rw = rw_block [0] = /*this->mult * */base_t::compute_gas_pref_deriv (data, params);
          }
        if (mobility_calc_t::is_g_ctrl (params))
          {
            ww_value[0] += rw/* * this->mult*/;
            wr_block[0] += (rr_block[0]/* * this->mult*/);
          }
      }
      void
      update_wr (const sp_connection_t &c, item_t ww) const
      {
        c->get_wr_value ()[0] *= ww;
      }
      void
      apply_wefac (const sp_connection_t &c, params_t &params) const
      {
        shared_vector <item_t> rw = c->get_rw_value ();
        item_t wefac = params.well_->exploitation_factor_;

        rw[0] = blue_sky::apply_wefac (rw[0], wefac);
      }
      void
      update_rate (const sp_connection_t &c, params_t &params) const
      {
        item_t ww_bw = (params.bw_value[0] / params.ww_value[0]);

        c->get_rate_value ()[0] -= c->get_rw_value () [0] * ww_bw;
      }
    };



  } // namespace wells
} // namespace blue_sky

#endif  // #ifndef BS_WELLS_WELL_RATE_CONTROL_RATE_DERIV_H_

