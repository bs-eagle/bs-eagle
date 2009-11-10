/** 
 * \file well_rate_compute_params.h
 * \brief params for compute well derivs
 * \author Sergey Miryanov
 * \date 14.01.2009
 * */
#ifndef BS_WELLS_WELL_RATE_COMPUTE_PARAMS_H_
#define BS_WELLS_WELL_RATE_COMPUTE_PARAMS_H_

namespace blue_sky {

  template <typename strategy_t>
  struct compute_params
  {
    typedef typename strategy_t::index_t      index_t;
    typedef typename strategy_t::item_t       item_t;
    typedef typename strategy_t::item_array_t item_array_t;
    typedef typename strategy_t::index_array_t  index_array_t;

    typedef calc_model <strategy_t>             calc_model_t;
    typedef well <strategy_t>                   well_t;
    typedef wells::well_controller <strategy_t>        well_controller_t;
    typedef jacobian_matrix <strategy_t>        jmatrix_t;

    typedef typename calc_model_t::data_t       data_t;
    typedef typename calc_model_t::data_array_t data_array_t;
    typedef typename calc_model_t::main_var_array_t main_var_array_t;
    typedef typename calc_model_t::phase_d_t    phase_d_t;

    typedef typename well_t::rate_data_t        rate_data_t;

    typedef smart_ptr <calc_model_t>            sp_calc_model_t;
    typedef smart_ptr <well_t>                  sp_well_t;
    typedef smart_ptr <const well_controller_t> sp_well_controller_t;
    typedef smart_ptr <jmatrix_t>               sp_jmatrix_t;

    typedef typename well_t::sp_connection_t    sp_connection_t;

    struct inj_params
    {
      inj_params (wells::injection_type injection)
        : is_o_ctrl (injection == wells::injection_oil),
        is_w_ctrl (injection == wells::injection_water),
        is_g_ctrl (injection == wells::injection_gas),
        is_oil_injection (injection == wells::injection_oil),
        is_water_injection (injection == wells::injection_water),
        is_gas_injection (injection == wells::injection_gas)
      {
      }

      bool    is_o_ctrl;
      bool    is_w_ctrl;
      bool    is_g_ctrl;
      bool    is_oil_injection;
      bool    is_water_injection;
      bool    is_gas_injection;

      item_t  krow_tetaow;
      item_t  krp_tetap;
      item_t  krp_dtetap_dpp;
      item_t  dkrow_dsw_tetaow_krw_dtetaw_dsw;
      item_t  dkrog_dsg_tetaog_krg_dtetag_dsg;
    };

    struct prod_params
    {
      prod_params (wells::rate_control_type control)
        : is_o_ctrl (control == wells::oil_rate_control || control == wells::liquid_rate_control),
        is_w_ctrl (control == wells::water_rate_control || control == wells::liquid_rate_control),
        is_g_ctrl (control == wells::gas_rate_control)
      {
      }

      bool is_o_ctrl;
      bool is_w_ctrl;
      bool is_g_ctrl;
    };

    compute_params (const sp_calc_model_t       &calc_model,
                    sp_jmatrix_t                &jmatrix,
                    sp_well_t                   &well,
                    const sp_well_controller_t  &well_controller)
      : calc_model_ (calc_model),
      jmatrix_ (jmatrix),
      well_ (well),
      well_controller_ (well_controller),
      main_vars (calc_model_->main_variable),
      pressure (calc_model_->pressure),
      gas_oil_ratio (calc_model_->gas_oil_ratio),
      gravity (calc_model_->internal_constants.gravity_constant),
      n_block (-1),
      perf_bhp (0),
      diff_depth (0),
      depth (0), 
      rho (0),
      main_var (FI_NULL),
      data_array (calc_model_->data),
      phase_d (calc_model_->phase_d),
      n_phases (calc_model_->n_phases),
      is_o (calc_model_->is_oil ()),
      is_w (calc_model_->is_water ()),
      is_g (calc_model_->is_gas ()),
      bw_value (well_->get_bw_value ()),
      ww_value (well_->get_ww_value ()),
      rate (well_->rate_),
      limit_rate (well_controller_->rate ()),
      inj_params_ (inj_params (well_controller_->injection ())),
      prod_params_ (prod_params (well_controller_->get_control_type ())),
      is_prod (well_controller_->is_production ()),
      Po (0),
      Pw (0),
      Pg (0),
      gw (0)
    {
    }

    void
    compute_perf_vars (const data_t &data, inj_params &params)
    {
      params.krow_tetaow
        = (is_o ? RELATIVE_PERM (data, phase_d, FI_PHASE_OIL)   * INVERS_VISCOSITY (data, phase_d, FI_PHASE_OIL) : 0)
        + (is_w ? RELATIVE_PERM (data, phase_d, FI_PHASE_WATER) * INVERS_VISCOSITY (data, phase_d, FI_PHASE_WATER) : 0);

      params.krp_tetap
        = params.krow_tetaow
        + (is_g ? RELATIVE_PERM (data, phase_d, FI_PHASE_GAS) * INVERS_VISCOSITY (data, phase_d, FI_PHASE_GAS) : 0);

      params.krp_dtetap_dpp
        = (is_o ? RELATIVE_PERM (data, phase_d, FI_PHASE_OIL)   * P_DERIV_INVERS_VISCOSITY (data, phase_d, FI_PHASE_OIL)    : 0)
        + (is_w ? RELATIVE_PERM (data, phase_d, FI_PHASE_WATER) * P_DERIV_INVERS_VISCOSITY (data, phase_d, FI_PHASE_WATER)  : 0)
        + (is_g ? RELATIVE_PERM (data, phase_d, FI_PHASE_GAS)   * P_DERIV_INVERS_VISCOSITY (data, phase_d, FI_PHASE_GAS)    : 0);

      params.dkrow_dsw_tetaow_krw_dtetaw_dsw
        = (is_o && is_w ? S_DERIV_RELATIVE_PERM (data, phase_d, n_phases, FI_PHASE_OIL, FI_PHASE_WATER)   * INVERS_VISCOSITY (data, phase_d, FI_PHASE_OIL) : 0)
        + (is_w         ? S_DERIV_RELATIVE_PERM (data, phase_d, n_phases, FI_PHASE_WATER, FI_PHASE_WATER) * INVERS_VISCOSITY (data, phase_d, FI_PHASE_WATER) : 0)
        + (is_w         ? RELATIVE_PERM (data, phase_d, FI_PHASE_WATER) * S_DERIV_INVERS_VISCOSITY (data, phase_d, FI_PHASE_WATER) * CAP_PRESSURE (data, phase_d, FI_PHASE_WATER) : 0)
        - (is_w         ? P_DERIV_INVERS_FVF (data, phase_d, FI_PHASE_WATER) * CAP_PRESSURE (data, phase_d, FI_PHASE_WATER) : 0)
        ;

      params.dkrog_dsg_tetaog_krg_dtetag_dsg
        = (is_o && is_g ? S_DERIV_RELATIVE_PERM (data, phase_d, n_phases, FI_PHASE_OIL, FI_PHASE_GAS) * INVERS_VISCOSITY (data, phase_d, FI_PHASE_OIL) : 0)
        + (is_g         ? S_DERIV_RELATIVE_PERM (data, phase_d, n_phases, FI_PHASE_GAS, FI_PHASE_GAS) * INVERS_VISCOSITY (data, phase_d, FI_PHASE_GAS) : 0)
        + (is_g         ? RELATIVE_PERM (data, phase_d, FI_PHASE_GAS) * S_DERIV_INVERS_VISCOSITY (data, phase_d, FI_PHASE_GAS) * CAP_PRESSURE (data, phase_d, FI_PHASE_GAS) : 0)
        ;
    }

public:
    const sp_calc_model_t       &calc_model_;
    sp_jmatrix_t                &jmatrix_;
    sp_well_t                   &well_;
    const sp_well_controller_t  &well_controller_;
    const main_var_array_t    &main_vars;
    const item_array_t        &pressure;
    const item_array_t        &gas_oil_ratio;
    item_t                    gravity;
    index_t                   n_block;
    item_t                    perf_bhp;
    item_t                    diff_depth;
    item_t                    depth;
    item_t                    rho;
    main_var_type             main_var;
    const data_array_t        &data_array;
    const phase_d_t           &phase_d;
    index_t                   n_phases;
    bool                      is_o;
    bool                      is_w;
    bool                      is_g;
    array_ext <item_t>        bw_value;
    array_ext <item_t>        ww_value;
    rate_data_t               &rate;
    const rate_data_t         &limit_rate;
    inj_params                inj_params_;
    prod_params               prod_params_;
    bool                      is_prod;
    
    item_t                    Po;
    item_t                    Pw;
    item_t                    Pg;
    item_t                    gw;
  };


} // namespace blue_sky


#endif  // #ifndef BS_WELLS_WELL_RATE_COMPUTE_PARAMS_H_

