/**
 *       \file  well_rate_control_bhp_deriv.h
 *      \brief  Calculates well derivs with control by BHP
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  21.11.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Obsolete, should be removed
 * */
#ifndef BS_WELLS_WELL_RATE_CONTROL_BHP_DERIV_H_
#define BS_WELLS_WELL_RATE_CONTROL_BHP_DERIV_H_

#include "well_rate_control_deriv.h"
#include "well_rate_control_compute_deriv_typedef.h"

namespace blue_sky
  {
  namespace wells
    {

    template <typename mobility_calc_t>
    struct compute_bhp_deriv_3p : public compute_deriv <mobility_calc_t>
      {
        GET_COMPUTE_DERIV_BASE_TYPES;

        compute_bhp_deriv_3p (const mobility_calc_t &mobility_calc)
            : base_t (mobility_calc)
        {
        }

        void
        oil_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            if (mobility_calc_t::is_oil_injection (params))
              {
                shared_vector <item_t> rr_block = c->get_rr_value ();
                shared_vector <item_t> ps_block = c->get_ps_value ();

                ps_block [p3_oil]     = base_t::compute_oil_sw_deriv (data, params);
                rr_block [p3_oil_po]  = base_t::compute_oil_po_deriv (data, params);
                rr_block [p3_oil_so]  = base_t::compute_oil_so_deriv (data, params);
                rr_block [p3_oil_sg]  = base_t::compute_oil_sg_deriv (data, params);
              }
          }
        void
        water_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            if (mobility_calc_t::is_water_injection (params))
              {
                shared_vector <item_t> rr_block = c->get_rr_value ();
                shared_vector <item_t> ps_block = c->get_ps_value ();

                ps_block [p3_wat]     = base_t::compute_water_sw_deriv (data, params);
                rr_block [p3_wat_po]  = base_t::compute_water_po_deriv (data, params);
                rr_block [p3_wat_so]  = base_t::compute_water_so_deriv (data, params);
                rr_block [p3_wat_sg]  = base_t::compute_water_sg_deriv (data, params);
              }
          }
        void
        gas_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            if (mobility_calc_t::is_gas_injection (params))
              {
                shared_vector <item_t> rr_block = c->get_rr_value ();
                shared_vector <item_t> ps_block = c->get_ps_value ();

                ps_block [p3_gas]     = base_t::compute_gas_sw_deriv (data, params);
                rr_block [p3_gas_po]  = base_t::compute_gas_po_deriv (data, params);
                rr_block [p3_gas_so]  = base_t::compute_gas_so_deriv (data, params);
                rr_block [p3_gas_sg]  = base_t::compute_gas_sg_deriv (data, params);
              }
          }

        void
        update_rr (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            shared_vector <item_t> rr = c->get_rr_value ();
            shared_vector <item_t> ps = c->get_ps_value ();
            const rhs_item_t *sp    = &params.jmatrix_->get_sp_diagonal ()[params.n_block * params.n_phases];

            M_MINUS_VV_PROD_3x3 (ps, sp, rr);
          }
        void
        update_rhs_flux (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            shared_vector <item_t> ps         = c->get_ps_value ();
            const rhs_item_t *s_rhs_block   = &params.jmatrix_->get_sec_rhs ()[params.calc_model_->n_sec_vars * params.n_block];
            shared_vector <rhs_item_t> r_rhs_block  = c->get_rate_value ();

            V_MINUS_VS_PROD_3x3 (ps, s_rhs_block, r_rhs_block);
          }

        void
        compute_bw_value (params_t &params) const
        {
          params.bw_value[0] = 0;

          if (mobility_calc_t::is_o_ctrl (params) && mobility_calc_t::is_w_ctrl (params))
            {
              if (params.is_prod)
                params.bw_value[0] += (params.rate.prod.oil + params.rate.prod.water) - this->mult * params.limit_rate.prod.liquid;
              else
                params.bw_value[0] += (params.rate.inj.oil + params.rate.inj.water) - this->mult * params.limit_rate.inj.liquid;
            }
          else
            {
              if (mobility_calc_t::is_o_ctrl (params))
                {
                  if (params.is_prod)
                    params.bw_value[0] += params.rate.prod.oil - this->mult * params.limit_rate.prod.oil;
                  else
                    params.bw_value[0] += params.rate.inj.oil - this->mult * params.limit_rate.inj.oil;
                }

              if (mobility_calc_t::is_w_ctrl (params))
                {
                  if (params.is_prod)
                    params.bw_value[0] += params.rate.prod.water - this->mult * params.limit_rate.prod.water;
                  else
                    params.bw_value[0] += params.rate.inj.water - this->mult * params.limit_rate.inj.water;
                }

              if (mobility_calc_t::is_g_ctrl (params))
                {
                  if (params.is_prod)
                    params.bw_value[0] += params.rate.prod.gas - this->mult * params.limit_rate.prod.gas;
                  else
                    params.bw_value[0] += params.rate.inj.gas - this->mult * params.limit_rate.inj.gas;
                }
            }
        }
        void
        apply_wefac (const sp_connection_t &c, params_t &params) const
        {
          shared_vector <item_t> rr = c->get_rr_value ();
          item_t wefac = params.well_->exploitation_factor_;

          rr[0] = blue_sky::apply_wefac (rr[0], wefac);
          rr[1] = blue_sky::apply_wefac (rr[1], wefac);
          rr[2] = blue_sky::apply_wefac (rr[2], wefac);
          rr[3] = blue_sky::apply_wefac (rr[3], wefac);
          rr[4] = blue_sky::apply_wefac (rr[4], wefac);
          rr[5] = blue_sky::apply_wefac (rr[5], wefac);
          rr[6] = blue_sky::apply_wefac (rr[6], wefac);
          rr[7] = blue_sky::apply_wefac (rr[7], wefac);
          rr[8] = blue_sky::apply_wefac (rr[8], wefac);
        }
      };

    template <typename mobility_calc_t>
    struct compute_bhp_deriv_2p_ow : public compute_deriv <mobility_calc_t>
      {
        GET_COMPUTE_DERIV_BASE_TYPES;

        compute_bhp_deriv_2p_ow (const mobility_calc_t &mobility_calc)
            : base_t (mobility_calc)
        {
        }

        void
        oil_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            if (mobility_calc_t::is_oil_injection (params))
              {
                shared_vector <item_t> rr_block = c->get_rr_value ();
                shared_vector <item_t> ps_block = c->get_ps_value ();

                ps_block [p2ow_oil]     = base_t::compute_oil_sw_deriv (data, params);
                rr_block [p2ow_oil_po]  = base_t::compute_oil_po_deriv (data, params);
                rr_block [p2ow_oil_so]  = base_t::compute_oil_so_deriv (data, params);
              }
          }
        void
        water_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            if (mobility_calc_t::is_water_injection (params))
              {
                shared_vector <item_t> rr_block = c->get_rr_value ();
                shared_vector <item_t> ps_block = c->get_ps_value ();

                ps_block [p2ow_wat]     = base_t::compute_water_sw_deriv (data, params);
                rr_block [p2ow_wat_po]  = base_t::compute_water_po_deriv (data, params);
                rr_block [p2ow_wat_so]  = base_t::compute_water_so_deriv (data, params);
              }
          }
        void
        gas_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
          }

        void
        update_rr (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            shared_vector <item_t> rr = c->get_rr_value ();
            shared_vector <item_t> ps = c->get_ps_value ();
            const rhs_item_t *sp  = &params.jmatrix_->get_sp_diagonal ()[params.n_block * params.n_phases];

            M_MINUS_VV_PROD_2x2 (ps, sp, rr);
          }
        void
        update_rhs_flux (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            shared_vector <item_t> ps               = c->get_ps_value ();
            const rhs_item_t *s_rhs_block       = &params.jmatrix_->get_sec_rhs ()[params.calc_model_->n_sec_vars * params.n_block];
            shared_vector <rhs_item_t> r_rhs_block  = c->get_rate_value ();

            V_MINUS_VS_PROD_2x2 (ps, s_rhs_block, r_rhs_block);
          }
        void
        compute_bw_value (params_t &params) const
          {
            params.bw_value[0] = 0;

            if (mobility_calc_t::is_o_ctrl (params))
              {
                if (params.is_prod)
                  params.bw_value[0] += params.rate.prod.oil - this->mult * params.limit_rate.prod.oil;
                else
                  params.bw_value[0] += params.rate.inj.oil - this->mult * params.limit_rate.inj.oil;
              }

            if (mobility_calc_t::is_w_ctrl (params))
              {
                if (params.is_prod)
                  params.bw_value[0] += params.rate.prod.water - this->mult * params.limit_rate.prod.water;
                else
                  params.bw_value[0] += params.rate.inj.water - this->mult * params.limit_rate.inj.water;
              }
          }
        void
        apply_wefac (const sp_connection_t &c, params_t &params) const
        {
          shared_vector <item_t> rr = c->get_rr_value ();
          item_t wefac = params.well_->exploitation_factor_;

          rr[0] = blue_sky::apply_wefac (rr[0], wefac);
          rr[1] = blue_sky::apply_wefac (rr[1], wefac);
          rr[2] = blue_sky::apply_wefac (rr[2], wefac);
          rr[3] = blue_sky::apply_wefac (rr[3], wefac);
        }
      };

    template <typename mobility_calc_t>
    struct compute_bhp_deriv_2p_og : public compute_deriv <mobility_calc_t>
      {
        GET_COMPUTE_DERIV_BASE_TYPES;

        compute_bhp_deriv_2p_og (const mobility_calc_t &mobility_calc)
            : base_t (mobility_calc)
        {
        }

        void
        oil_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            if (mobility_calc_t::is_oil_injection (params))
              {
                shared_vector <item_t> rr_block = c->get_rr_value ();
                shared_vector <item_t> ps_block = c->get_ps_value ();

                ps_block [p2og_oil]     = base_t::compute_oil_so_deriv (data, params);
                rr_block [p2og_oil_po]  = base_t::compute_oil_po_deriv (data, params);
                rr_block [p2og_oil_sg]  = base_t::compute_oil_sg_deriv (data, params);
              }
          }
        void
        water_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
          }
        void
        gas_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            if (mobility_calc_t::is_gas_injection (params))
              {
                shared_vector <item_t> rr_block = c->get_rr_value ();
                shared_vector <item_t> ps_block = c->get_ps_value ();

                ps_block [p2og_gas]     = base_t::compute_gas_so_deriv (data, params);
                rr_block [p2og_gas_po]  = base_t::compute_gas_po_deriv (data, params);
                rr_block [p2og_gas_sg]  = base_t::compute_gas_sg_deriv (data, params);
              }
          }

        void
        update_rr (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            shared_vector <item_t> rr = c->get_rr_value ();
            shared_vector <item_t> ps = c->get_ps_value ();
            const rhs_item_t *sp    = &params.jmatrix_->get_sp_diagonal ()[params.n_block * params.n_phases];

            M_MINUS_VV_PROD_2x2 (ps, sp, rr);
          }
        void
        update_rhs_flux (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            shared_vector <item_t> ps               = c->get_ps_value ();
            const rhs_item_t *s_rhs_block       = &params.jmatrix_->get_sec_rhs ()[params.calc_model_->n_sec_vars * params.n_block];
            shared_vector <rhs_item_t> r_rhs_block  = c->get_rate_value ();

            V_MINUS_VS_PROD_2x2 (ps, s_rhs_block, r_rhs_block);
          }
        void
        compute_bw_value (params_t &params) const
          {
            params.bw_value[0] = 0;

            if (mobility_calc_t::is_o_ctrl (params))
              {
                if (params.is_prod)
                  params.bw_value[0] += params.rate.prod.oil - this->mult * params.limit_rate.prod.oil;
                else
                  params.bw_value[0] += params.rate.inj.oil - this->mult * params.limit_rate.inj.oil;
              }

            if (mobility_calc_t::is_g_ctrl (params))
              {
                if (params.is_prod)
                  params.bw_value[0] += params.rate.prod.gas - this->mult * params.limit_rate.prod.gas;
                else
                  params.bw_value[0] += params.rate.inj.gas - this->mult * params.limit_rate.inj.gas;
              }
          }
        void
        apply_wefac (const sp_connection_t &c, params_t &params) const
        {
          shared_vector <item_t> rr = c->get_rr_value ();
          item_t wefac = params.well_->exploitation_factor_;

          rr[0] = blue_sky::apply_wefac (rr[0], wefac);
          rr[1] = blue_sky::apply_wefac (rr[1], wefac);
          rr[2] = blue_sky::apply_wefac (rr[2], wefac);
          rr[3] = blue_sky::apply_wefac (rr[3], wefac);
        }
      };

    template <typename mobility_calc_t>
    struct compute_bhp_deriv_1p_o : public compute_deriv <mobility_calc_t>
      {
        GET_COMPUTE_DERIV_BASE_TYPES;

        compute_bhp_deriv_1p_o (const mobility_calc_t &mobility_calc)
            : base_t (mobility_calc)
        {
        }

        void
        oil_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            if (mobility_calc_t::is_oil_injection (params))
              {
                shared_vector <item_t> rr_block = c->get_rr_value ();
                rr_block [0] = base_t::compute_oil_po_deriv (data, params);
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
        update_rr (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
          }
        void
        update_rhs_flux (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            shared_vector <item_t> ps               = c->get_ps_value ();
            const rhs_item_t *s_rhs_block       = &params.jmatrix_->get_sec_rhs ()[params.calc_model_->n_sec_vars * params.n_block];
            shared_vector <rhs_item_t> r_rhs_block  = c->get_rate_value ();

            V_MINUS_VS_PROD_1x1 (ps, s_rhs_block, r_rhs_block);
          }
        void
        compute_bw_value (params_t &params) const
          {
            BS_ASSERT (mobility_calc_t::is_o_ctrl (params));

            if (params.is_prod)
              params.bw_value[0] += params.rate.prod.oil - this->mult * params.limit_rate.prod.oil;
            else
              params.bw_value[0] += params.rate.inj.oil - this->mult * params.limit_rate.inj.oil;
          }
        void
        apply_wefac (const sp_connection_t &c, params_t &params) const
        {
          shared_vector <item_t> rr = c->get_rr_value ();
          item_t wefac = params.well_->exploitation_factor_;

          rr[0] = blue_sky::apply_wefac (rr[0], wefac);
        }
      };
    template <typename mobility_calc_t>
    struct compute_bhp_deriv_1p_w : public compute_deriv <mobility_calc_t>
      {
        GET_COMPUTE_DERIV_BASE_TYPES;

        compute_bhp_deriv_1p_w (const mobility_calc_t &mobility_calc)
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
            if (mobility_calc_t::is_water_injection (params))
              {
                shared_vector <item_t> rr_block = c->get_rr_value ();
                rr_block [0] = base_t::compute_water_po_deriv (data, params);
              }
          }
        void
        gas_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
          }

        void
        update_rr (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
          }
        void
        update_rhs_flux (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            shared_vector <item_t> ps               = c->get_ps_value ();
            const rhs_item_t *s_rhs_block       = &params.jmatrix_->get_sec_rhs ()[params.calc_model_->n_sec_vars * params.n_block];
            shared_vector <rhs_item_t> r_rhs_block  = c->get_rate_value ();

            V_MINUS_VS_PROD_1x1 (ps, s_rhs_block, r_rhs_block);
          }
        void
        compute_bw_value (params_t &params) const
          {
            BS_ASSERT (mobility_calc_t::is_w_ctrl (params));

            if (params.is_prod)
              params.bw_value[0] += params.rate.prod.water - this->mult * params.limit_rate.prod.water;
            else
              params.bw_value[0] += params.rate.inj.water - this->mult * params.limit_rate.inj.water;
          }
        void
        apply_wefac (const sp_connection_t &c, params_t &params) const
        {
          shared_vector <item_t> rr = c->get_rr_value ();
          item_t wefac = params.well_->exploitation_factor_;

          rr[0] = blue_sky::apply_wefac (rr[0], wefac);
        }
      };
    template <typename mobility_calc_t>
    struct compute_bhp_deriv_1p_g : public compute_deriv <mobility_calc_t>
      {
        GET_COMPUTE_DERIV_BASE_TYPES;

        compute_bhp_deriv_1p_g (const mobility_calc_t &mobility_calc)
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
            if (mobility_calc_t::is_gas_injection (params))
              {
                shared_vector <item_t> rr_block = c->get_rr_value ();
                rr_block [0] = base_t::compute_gas_po_deriv (data, params);
              }
          }

        void
        update_rr (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
          }
        void
        update_rhs_flux (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            shared_vector <item_t> ps               = c->get_ps_value ();
            const rhs_item_t *s_rhs_block       = &params.jmatrix_->get_sec_rhs ()[params.calc_model_->n_sec_vars * params.n_block];
            shared_vector <rhs_item_t> r_rhs_block  = c->get_rate_value ();

            V_MINUS_VS_PROD_1x1 (ps, s_rhs_block, r_rhs_block);
          }
        void
        compute_bw_value (params_t &params) const
          {
            BS_ASSERT (mobility_calc_t::is_g_ctrl (params));

            if (params.is_prod)
              params.bw_value[0] += params.rate.prod.gas - this->mult * params.limit_rate.prod.gas;
            else
              params.bw_value[0] += params.rate.inj.gas - this->mult * params.limit_rate.inj.gas;
          }
        void
        apply_wefac (const sp_connection_t &c, params_t &params) const
        {
          shared_vector <item_t> rr = c->get_rr_value ();
          item_t wefac = params.well_->exploitation_factor_;

          rr[0] = blue_sky::apply_wefac (rr[0], wefac);
        }
      };

  } // namespace wells
} // namespace blue_sky


#endif  // #ifndef BS_WELLS_WELL_RATE_CONTROL_BHP_DERIV_H_

