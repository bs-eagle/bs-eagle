/**
 *       \file  well_rate_control_impl.h
 *      \brief  Implementation of well_rate_control
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  21.11.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Obsolete, should be removed
 * */
#ifndef BS_WELLS_WELL_RATE_CONTROL_IMPL_H_
#define BS_WELLS_WELL_RATE_CONTROL_IMPL_H_

#include "well_rate_control_interface.h"
#include "well_rate_compute_potentials.h"
#include "well_rate_call_proxy.h"
#include "well_rate_connection_loop.h"

#include "well_rate_compute_params.h"

namespace blue_sky
  {

  template <typename impl_type_t>
  class BS_API_PLUGIN well_rate_control_impl : public wells::well_rate_control_interface <typename impl_type_t::strategy_t>
    {
    public:

      template <typename impl_t>
      struct impl;

      typedef typename impl_type_t::strategy_t                strategy_t;
      typedef impl <typename impl_type_t::inj_impl_t>         inj_impl_t;
      typedef impl <typename impl_type_t::prod_impl_t>        prod_impl_t;

      typedef well_rate_control_impl  <impl_type_t>           this_t;
      typedef wells::well_rate_control_interface <strategy_t> base_t;

      typedef typename base_t::calc_model_t                   calc_model_t;
      typedef typename calc_model_t::data_t                   data_t;
      typedef compute_params <strategy_t>                     compute_params_t;

      typedef typename base_t::sp_calc_model_t                sp_calc_model_t;
      typedef typename base_t::sp_jmatrix_t                   sp_jmatrix_t;
      typedef typename base_t::sp_well_t                      sp_well_t;
      typedef typename base_t::sp_well_controller_t           sp_well_controller_t;

      typedef typename base_t::sp_connection_t                sp_connection_t;

      typedef typename base_t::item_t                         item_t;
      typedef typename base_t::index_t                        index_t;

      void
      compute_rate (const sp_calc_model_t &calc_model, sp_jmatrix_t &jmatrix, sp_well_t &well, const sp_well_controller_t &well_controller) const
      {
        compute_params_t params (calc_model, jmatrix, well, well_controller);
        connection_loop (params, 
          one_call (&inj_impl_, &inj_impl_t::compute_rate), 
          one_call (&prod_impl_, &prod_impl_t::compute_rate));

        compute_bw_value (params);
        if (params.well_->is_shut () || params.well_->exploitation_factor_ <= 0.0)
          {
            if (params.is_prod)
              {
                apply_wefac_connection_loop (params, one_call (&prod_impl_, &prod_impl_t::apply_wefac));
              }
            else
              {
                apply_wefac_connection_loop (params, one_call (&inj_impl_, &inj_impl_t::apply_wefac));
              }
          }
      }
      void
      compute_derivs (const sp_calc_model_t &calc_model, sp_jmatrix_t &jmatrix, sp_well_t &well, const sp_well_controller_t &well_controller) const
      {
        compute_params_t params (calc_model, jmatrix, well, well_controller);
        connection_loop (params, 
          two_call (&inj_impl_, &inj_impl_t::compute_bhp_derivs, &inj_impl_t::compute_rate_derivs),
          two_call (&prod_impl_, &prod_impl_t::compute_bhp_derivs, &prod_impl_t::compute_rate_derivs));

        if (params.is_prod)
          {
            update_wr_connection_loop (params, one_call (&prod_impl_, &prod_impl_t::update_wr));
            apply_wefac_connection_loop (params, one_call (&prod_impl_, &prod_impl_t::apply_wefac));
            apply_wefac_connection_loop (params, one_call (&prod_impl_, &prod_impl_t::update_rate));
          }
        else
          {
            update_wr_connection_loop (params, one_call (&inj_impl_, &inj_impl_t::update_wr));
            apply_wefac_connection_loop (params, one_call (&inj_impl_, &inj_impl_t::apply_wefac));
            apply_wefac_connection_loop (params, one_call (&inj_impl_, &inj_impl_t::update_rate));
          }
      }

      BLUE_SKY_TYPE_DECL_T (well_rate_control_impl);

    public:

      template <typename impl_t>
      struct impl
      {
        typedef typename impl_t::mobility_t     mobility_t;
        typedef typename impl_t::rate_t         rate_t;
        typedef typename impl_t::bhp_deriv_t    bhp_deriv_t;
        typedef typename impl_t::rate_deriv_t   rate_deriv_t;
        typedef compute_params_t                params_t;

        impl ()
          : rate_ (mobility_)
          , bhp_deriv_ (mobility_)
          , rate_deriv_ (mobility_)
        {

        }

        void
        compute_rate (const sp_connection_t &c, const data_t &data, params_t &params) const
        {
          rate_.oil_function (c, data, params);
          rate_.water_function (c, data, params);
          rate_.gas_function (c, data, params);
        }

        void
        compute_bhp_derivs (const sp_connection_t &c, const data_t &data, params_t &params) const
        {
          bhp_deriv_.oil_function (c, data, params);
          bhp_deriv_.water_function (c, data, params);
          bhp_deriv_.gas_function (c, data, params);

          bhp_deriv_.update_rr (c, data, params);
          bhp_deriv_.update_rhs_flux (c, data, params);
        }
        void
        compute_rate_derivs (const sp_connection_t &c, const data_t &data, params_t &params) const
        {
          rate_deriv_.oil_function (c, data, params);
          rate_deriv_.water_function (c, data, params);
          rate_deriv_.gas_function (c, data, params);
        }

        void
        update_wr (const sp_connection_t &c, item_t ww) const
        {
          rate_deriv_.update_wr (c, ww);
        }

        void
        apply_wefac (const sp_connection_t &c, params_t &params) const
        {
          bhp_deriv_.apply_wefac (c, params);
          rate_deriv_.apply_wefac (c, params);
        }

        void
        update_rate (const sp_connection_t &c, params_t &params) const
        {
          rate_deriv_.update_rate (c, params);
        }

        mobility_t          mobility_;
        rate_t              rate_;
        bhp_deriv_t         bhp_deriv_;
        rate_deriv_t        rate_deriv_;
      };

    private:


      void
      compute_bw_value (compute_params_t &params) const
      {
        params.is_prod 
          ? prod_impl_.bhp_deriv_.compute_bw_value (params)
          : inj_impl_.bhp_deriv_.compute_bw_value (params)
          ;
      }

      inj_impl_t  inj_impl_;
      prod_impl_t prod_impl_;
    };


} // namespace blue_sky


#endif  // #ifndef BS_WELLS_WELL_RATE_CONTROL_IMPL_H

