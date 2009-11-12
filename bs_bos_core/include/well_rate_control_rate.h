/**
 *       \file  well_rate_control_rate.h
 *      \brief  Calculates rates for well
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  24.11.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Obsolete, should be removed
 * */
#ifndef BS_WELLS_WELL_RATE_CONTROL_RATE_H_
#define BS_WELLS_WELL_RATE_CONTROL_RATE_H_

#include "well_rate_control_deriv.h"
#include "well_rate_control_compute_deriv_typedef.h"

namespace blue_sky
  {
  namespace wells
    {


    template <typename mobility_calc_t>
    struct compute_rate_3p : public compute_deriv <mobility_calc_t>
      {
        GET_COMPUTE_DERIV_BASE_TYPES;

        compute_rate_3p (const mobility_calc_t &mobility_calc)
            : base_t (mobility_calc)
        {
        }

        void
        oil_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            if (mobility_calc_t::is_oil_injection (params))
              {
                c->get_rate_value () [p3_oil] = base_t::compute_oil_rate (data, params);
                if (params.is_prod)
                  {
                    c->rate_.prod.oil        = c->get_rate_value () [p3_oil];
                    c->rate_.prod.liquid     = c->get_rate_value () [p3_oil];
                    params.rate.prod.oil    += c->get_rate_value () [p3_oil];
                    params.rate.prod.liquid += c->get_rate_value () [p3_oil];

                    c->rate_rc_.prod.oil                = c->rate_.prod.oil / params.calc_model_->invers_fvf_average[params.phase_d[FI_PHASE_OIL]];
                    c->rate_rc_.prod.liquid             = c->rate_rc_.prod.oil;
                    params.well_->rate_rc_.prod.oil    += c->rate_rc_.prod.oil;
                    params.well_->rate_rc_.prod.liquid += c->rate_rc_.prod.oil;
                  }
                else
                  {
                    c->rate_.inj.oil         = c->get_rate_value () [p3_oil];
                    c->rate_.inj.liquid      = c->get_rate_value () [p3_oil];
                    params.rate.inj.oil     += c->get_rate_value () [p3_oil];
                    params.rate.inj.liquid  += c->get_rate_value () [p3_oil];

                    c->rate_rc_.inj.oil                = c->rate_.inj.oil / params.calc_model_->invers_fvf_average[params.phase_d[FI_PHASE_OIL]];
                    c->rate_rc_.inj.liquid             = c->rate_rc_.inj.oil;
                    params.well_->rate_rc_.inj.oil    += c->rate_rc_.inj.oil;
                    params.well_->rate_rc_.inj.liquid += c->rate_rc_.inj.oil;
                  }
              }
          }
        void
        water_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            if (mobility_calc_t::is_water_injection (params))
              {
                c->get_rate_value () [p3_wat] = base_t::compute_water_rate (data, params);
                if (params.is_prod)
                  {
                    c->rate_.prod.water      = c->get_rate_value () [p3_wat];
                    c->rate_.prod.liquid     = c->get_rate_value () [p3_wat];
                    params.rate.prod.water  += c->get_rate_value () [p3_wat];
                    params.rate.prod.liquid += c->get_rate_value () [p3_wat];

                    c->rate_rc_.prod.water              = c->rate_.prod.water / params.calc_model_->invers_fvf_average[params.phase_d[FI_PHASE_WATER]];
                    c->rate_rc_.prod.liquid             = c->rate_rc_.prod.water;
                    params.well_->rate_rc_.prod.water  += c->rate_rc_.prod.water;
                    params.well_->rate_rc_.prod.liquid += c->rate_rc_.prod.water;
                  }
                else
                  {
                    c->rate_.inj.water       = c->get_rate_value () [p3_wat];
                    c->rate_.inj.liquid      = c->get_rate_value () [p3_wat];
                    params.rate.inj.water   += c->get_rate_value () [p3_wat];
                    params.rate.inj.liquid  += c->get_rate_value () [p3_wat];

                    c->rate_rc_.inj.water              = c->rate_.inj.water / params.calc_model_->invers_fvf_average[params.phase_d[FI_PHASE_WATER]];
                    c->rate_rc_.inj.liquid             = c->rate_rc_.inj.water;
                    params.well_->rate_rc_.inj.water  += c->rate_rc_.inj.water;
                    params.well_->rate_rc_.inj.liquid += c->rate_rc_.inj.water;
                  }
              }
          }
        void
        gas_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            if (mobility_calc_t::is_gas_injection (params))
              {
                c->get_rate_value () [p3_gas] = base_t::compute_gas_rate (data, params);
                if (params.is_prod)
                  {
                    c->rate_.prod.gas          = c->get_rate_value () [p3_gas];
                    c->rate_.free_gas          = base_t::compute_free_gas_rate (data, params);
                    c->rate_.solution_gas      = base_t::compute_solution_gas_rate (data, params);
                    params.rate.prod.gas      += c->get_rate_value () [p3_gas];
                    params.rate.free_gas      += c->rate_.free_gas;
                    params.rate.solution_gas  += c->rate_.solution_gas;
                    params.well_->gor_         = params.rate.solution_gas / params.rate.prod.oil;

                    c->rate_rc_.free_gas       = c->rate_.free_gas / params.calc_model_->invers_fvf_average[params.phase_d[FI_PHASE_GAS]];
                    c->rate_rc_.solution_gas   = c->rate_.solution_gas;
                    c->rate_rc_.prod.gas       = c->rate_.free_gas + c->rate_.solution_gas / INVERS_FVF (data, params.phase_d, FI_PHASE_GAS);
                    params.well_->rate_rc_.prod.gas += c->rate_rc_.prod.gas;
                    params.well_->rate_rc_.free_gas += c->rate_rc_.free_gas;
                  }
                else
                  {
                    c->rate_.inj.gas       = c->get_rate_value () [p3_gas];
                    params.rate.inj.gas   += c->get_rate_value () [p3_gas];

                    c->rate_rc_.inj.gas    = base_t::compute_free_gas_rate (data, params) / params.calc_model_->invers_fvf_average[params.phase_d[FI_PHASE_GAS]];
                    params.well_->rate_rc_.inj.gas += c->rate_rc_.inj.gas;
                  }
              }
          }
      };

    template <typename mobility_calc_t>
    struct compute_rate_2p_ow : public compute_deriv <mobility_calc_t>
      {
        GET_COMPUTE_DERIV_BASE_TYPES;

        compute_rate_2p_ow (const mobility_calc_t &mobility_calc)
            : base_t (mobility_calc)
        {
        }

        void
        oil_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            if (mobility_calc_t::is_oil_injection (params))
              {
                c->get_rate_value () [p2ow_oil] = base_t::compute_oil_rate (data, params);
                if (params.is_prod)
                  {
                    c->rate_.prod.oil        = c->get_rate_value () [p2ow_oil];
                    c->rate_.prod.liquid     = c->get_rate_value () [p2ow_oil];
                    params.rate.prod.oil    += c->get_rate_value () [p2ow_oil];
                    params.rate.prod.liquid += c->get_rate_value () [p2ow_oil];

                    c->rate_rc_.prod.oil                = c->rate_.prod.oil / params.calc_model_->invers_fvf_average[params.phase_d[FI_PHASE_OIL]];
                    c->rate_rc_.prod.liquid             = c->rate_rc_.prod.oil;
                    params.well_->rate_rc_.prod.oil    += c->rate_rc_.prod.oil;
                    params.well_->rate_rc_.prod.liquid += c->rate_rc_.prod.oil;
                  }
                else
                  {
                    c->rate_.inj.oil         = c->get_rate_value () [p2ow_oil];
                    c->rate_.inj.liquid      = c->get_rate_value () [p2ow_oil];
                    params.rate.inj.oil     += c->get_rate_value () [p2ow_oil];
                    params.rate.inj.liquid  += c->get_rate_value () [p2ow_oil];

                    c->rate_rc_.inj.oil                = c->rate_.inj.oil / params.calc_model_->invers_fvf_average[params.phase_d[FI_PHASE_OIL]];
                    c->rate_rc_.inj.liquid             = c->rate_rc_.inj.oil;
                    params.well_->rate_rc_.inj.oil    += c->rate_rc_.inj.oil;
                    params.well_->rate_rc_.inj.liquid += c->rate_rc_.inj.oil;
                  }
              }
          }
        void
        water_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            if (mobility_calc_t::is_water_injection (params))
              {
                c->get_rate_value () [p2ow_wat] = base_t::compute_water_rate (data, params);
                if (params.is_prod)
                  {
                    c->rate_.prod.water      = c->get_rate_value () [p2ow_wat];
                    c->rate_.prod.liquid     = c->get_rate_value () [p2ow_wat];
                    params.rate.prod.water  += c->get_rate_value () [p2ow_wat];
                    params.rate.prod.liquid += c->get_rate_value () [p2ow_wat];

                    c->rate_rc_.prod.water              = c->rate_.prod.water / params.calc_model_->invers_fvf_average[params.phase_d[FI_PHASE_WATER]];
                    c->rate_rc_.prod.liquid             = c->rate_rc_.prod.water;
                    params.well_->rate_rc_.prod.water  += c->rate_rc_.prod.water;
                    params.well_->rate_rc_.prod.liquid += c->rate_rc_.prod.water;
                  }
                else
                  {
                    c->rate_.inj.water       = c->get_rate_value () [p2ow_wat];
                    c->rate_.inj.liquid      = c->get_rate_value () [p2ow_wat];
                    params.rate.inj.water   += c->get_rate_value () [p2ow_wat];
                    params.rate.inj.liquid  += c->get_rate_value () [p2ow_wat];

                    c->rate_rc_.inj.water              = c->rate_.inj.water / params.calc_model_->invers_fvf_average[params.phase_d[FI_PHASE_WATER]];
                    c->rate_rc_.inj.liquid             = c->rate_rc_.inj.water;
                    params.well_->rate_rc_.inj.water  += c->rate_rc_.inj.water;
                    params.well_->rate_rc_.inj.liquid += c->rate_rc_.inj.water;
                  }
              }
          }
        void
        gas_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
          }
      };

    template <typename mobility_calc_t>
    struct compute_rate_2p_og : public compute_deriv <mobility_calc_t>
      {
        GET_COMPUTE_DERIV_BASE_TYPES;

        compute_rate_2p_og (const mobility_calc_t &mobility_calc)
            : base_t (mobility_calc)
        {
        }

        void
        oil_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            if (mobility_calc_t::is_oil_injection (params))
              {
                c->get_rate_value () [p2og_oil] = base_t::compute_oil_rate (data, params);
                if (params.is_prod)
                  {
                    c->rate_.prod.oil        = c->get_rate_value () [p2og_oil];
                    c->rate_.prod.liquid     = c->get_rate_value () [p2og_oil];
                    params.rate.prod.oil    += c->get_rate_value () [p2og_oil];
                    params.rate.prod.liquid += c->get_rate_value () [p2og_oil];

                    c->rate_rc_.prod.oil                = c->rate_.prod.oil / params.calc_model_->invers_fvf_average[params.phase_d[FI_PHASE_OIL]];
                    c->rate_rc_.prod.liquid             = c->rate_rc_.prod.oil;
                    params.well_->rate_rc_.prod.oil    += c->rate_rc_.prod.oil;
                    params.well_->rate_rc_.prod.liquid += c->rate_rc_.prod.oil;
                  }
                else
                  {
                    c->rate_.inj.oil         = c->get_rate_value () [p2og_oil];
                    c->rate_.inj.liquid      = c->get_rate_value () [p2og_oil];
                    params.rate.inj.oil     += c->get_rate_value () [p2og_oil];
                    params.rate.inj.liquid  += c->get_rate_value () [p2og_oil];

                    c->rate_rc_.inj.oil                = c->rate_.inj.oil / params.calc_model_->invers_fvf_average[params.phase_d[FI_PHASE_OIL]];
                    c->rate_rc_.inj.liquid             = c->rate_rc_.inj.oil;
                    params.well_->rate_rc_.inj.oil    += c->rate_rc_.inj.oil;
                    params.well_->rate_rc_.inj.liquid += c->rate_rc_.inj.oil;
                  }
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
                c->get_rate_value () [p2og_gas] = base_t::compute_gas_rate (data, params);
                c->rate_.free_gas               = base_t::compute_free_gas_rate (data, params) / INVERS_FVF (data, params.phase_d, FI_PHASE_GAS);
                if (params.is_prod)
                  {
                    c->rate_.prod.gas = c->get_rate_value () [p2og_gas];
                    params.rate.prod.gas += c->get_rate_value () [p2og_gas];

                    c->rate_.free_gas          = base_t::compute_free_gas_rate (data, params);
                    c->rate_.solution_gas      = base_t::compute_solution_gas_rate (data, params);

                    params.rate.free_gas      += c->rate_.free_gas;
                    params.rate.solution_gas  += c->rate_.solution_gas;
                    params.well_->gor_         = params.rate.solution_gas / params.rate.prod.oil;

                    c->rate_rc_.free_gas       = c->rate_.free_gas / params.calc_model_->invers_fvf_average[params.phase_d[FI_PHASE_GAS]];
                    c->rate_rc_.solution_gas   = c->rate_.solution_gas;
                    c->rate_rc_.prod.gas       = c->rate_.free_gas + c->rate_.solution_gas / INVERS_FVF (data, params.phase_d, FI_PHASE_GAS);
                    params.well_->rate_rc_.prod.gas += c->rate_rc_.prod.gas;
                    params.well_->rate_rc_.free_gas += c->rate_rc_.free_gas;
                  }
                else
                  {
                    c->rate_.inj.gas = c->get_rate_value () [p2og_gas];
                    params.rate.inj.gas += c->get_rate_value () [p2og_gas];

                    c->rate_rc_.inj.gas    = base_t::compute_free_gas_rate (data, params) / params.calc_model_->invers_fvf_average[params.phase_d[FI_PHASE_GAS]];
                    params.well_->rate_rc_.inj.gas += c->rate_rc_.inj.gas;
                  }
              }
          }
      };

    template <typename mobility_calc_t>
    struct compute_rate_1p_o : public compute_deriv <mobility_calc_t>
      {
        GET_COMPUTE_DERIV_BASE_TYPES;

        compute_rate_1p_o (const mobility_calc_t &mobility_calc)
            : base_t (mobility_calc)
        {
        }

        void
        oil_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
            if (mobility_calc_t::is_oil_injection (params))
              {
                c->get_rate_value () [0] = base_t::compute_oil_rate (data, params);
                if (params.is_prod)
                  {
                    c->rate_.prod.oil        = c->get_rate_value () [0];
                    c->rate_.prod.liquid     = c->get_rate_value () [0];
                    params.rate.prod.oil    += c->get_rate_value () [0];
                    params.rate.prod.liquid += c->get_rate_value () [0];

                    c->rate_rc_.prod.oil                = c->rate_.prod.oil / params.calc_model_->invers_fvf_average[params.phase_d[FI_PHASE_OIL]];
                    c->rate_rc_.prod.liquid             = c->rate_rc_.prod.oil;
                    params.well_->rate_rc_.prod.oil    += c->rate_rc_.prod.oil;
                    params.well_->rate_rc_.prod.liquid += c->rate_rc_.prod.oil;
                  }
                else
                  {
                    c->rate_.inj.oil         = c->get_rate_value () [0];
                    c->rate_.inj.liquid      = c->get_rate_value () [0];
                    params.rate.inj.oil     += c->get_rate_value () [0];
                    params.rate.inj.liquid  += c->get_rate_value () [0];

                    c->rate_rc_.inj.oil                = c->rate_.inj.oil / params.calc_model_->invers_fvf_average[params.phase_d[FI_PHASE_OIL]];
                    c->rate_rc_.inj.liquid             = c->rate_rc_.inj.oil;
                    params.well_->rate_rc_.inj.oil    += c->rate_rc_.inj.oil;
                    params.well_->rate_rc_.inj.liquid += c->rate_rc_.inj.oil;
                  }
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
      };
    template <typename mobility_calc_t>
    struct compute_rate_1p_w : public compute_deriv <mobility_calc_t>
      {
        GET_COMPUTE_DERIV_BASE_TYPES;

        compute_rate_1p_w (const mobility_calc_t &mobility_calc)
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
                c->get_rate_value () [0] = base_t::compute_water_rate (data, params);
                if (params.is_prod)
                  {
                    c->rate_.prod.water      = c->get_rate_value () [0];
                    c->rate_.prod.liquid     = c->get_rate_value () [0];
                    params.rate.prod.water  += c->get_rate_value () [0];
                    params.rate.prod.liquid += c->get_rate_value () [0];

                    c->rate_rc_.prod.water              = c->rate_.prod.water / params.calc_model_->invers_fvf_average[params.phase_d[FI_PHASE_WATER]];
                    c->rate_rc_.prod.liquid             = c->rate_rc_.prod.water;
                    params.well_->rate_rc_.prod.water  += c->rate_rc_.prod.water;
                    params.well_->rate_rc_.prod.liquid += c->rate_rc_.prod.water;
                  }
                else
                  {
                    c->rate_.inj.water       = c->get_rate_value () [0];
                    c->rate_.inj.liquid      = c->get_rate_value () [0];
                    params.rate.inj.water   += c->get_rate_value () [0];
                    params.rate.inj.liquid  += c->get_rate_value () [0];

                    c->rate_rc_.inj.water              = c->rate_.inj.water / params.calc_model_->invers_fvf_average[params.phase_d[FI_PHASE_WATER]];
                    c->rate_rc_.inj.liquid             = c->rate_rc_.inj.water;
                    params.well_->rate_rc_.inj.water  += c->rate_rc_.inj.water;
                    params.well_->rate_rc_.inj.liquid += c->rate_rc_.inj.water;
                  }
              }
          }
        void
        gas_function (const sp_connection_t &c, const data_t &data, params_t &params) const
          {
          }
      };
    template <typename mobility_calc_t>
    struct compute_rate_1p_g : public compute_deriv <mobility_calc_t>
      {
        GET_COMPUTE_DERIV_BASE_TYPES;

        compute_rate_1p_g (const mobility_calc_t &mobility_calc)
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
                c->get_rate_value () [0] = base_t::compute_gas_rate (data, params);
                c->rate_.free_gas        = base_t::compute_free_gas_rate (data, params) / INVERS_FVF (data, params.phase_d, FI_PHASE_GAS);
                if (params.is_prod)
                  {
                    c->rate_.prod.gas          = c->get_rate_value () [0];
                    c->rate_.free_gas          = base_t::compute_free_gas_rate (data, params);
                    c->rate_.solution_gas      = base_t::compute_solution_gas_rate (data, params);
                    params.rate.prod.gas      += c->get_rate_value () [0];
                    params.rate.free_gas      += c->rate_.free_gas;
                    params.rate.solution_gas  += c->rate_.solution_gas;
                    params.well_->gor_         = params.rate.solution_gas / params.rate.prod.oil;

                    c->rate_rc_.free_gas       = c->rate_.free_gas / params.calc_model_->invers_fvf_average[params.phase_d[FI_PHASE_GAS]];
                    c->rate_rc_.solution_gas   = c->rate_.solution_gas;
                    c->rate_rc_.prod.gas       = c->rate_.free_gas + c->rate_.solution_gas / INVERS_FVF (data, params.phase_d, FI_PHASE_GAS);
                    params.well_->rate_rc_.prod.gas += c->rate_rc_.prod.gas;
                    params.well_->rate_rc_.free_gas += c->rate_rc_.free_gas;
                  }
                else
                  {
                    c->rate_.inj.gas       = c->get_rate_value () [0];
                    params.rate.inj.gas   += c->get_rate_value () [0];

                    c->rate_rc_.inj.gas    = base_t::compute_free_gas_rate (data, params) / params.calc_model_->invers_fvf_average[params.phase_d[FI_PHASE_GAS]];
                    params.well_->rate_rc_.inj.gas += c->rate_rc_.inj.gas;
                  }
              }
          }
      };

  } // namespace wells
} // namespace blue_sky

#endif // #ifndef BS_WELLS_WELL_RATE_CONTROL_RATE_H_

