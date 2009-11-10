/**
 * \file well_rate_control_inj_mobility.h
 * \brief mobility computation for injection wells
 * \author Sergey Miryanov
 * \date 21.11.2008
 * */
#ifndef BS_WELLS_WELL_RATE_INJ_MOBILITY_H_
#define BS_WELLS_WELL_RATE_INJ_MOBILITY_H_

#include "calc_model.h"
#include "calc_model_data_accessors.h"
#include "well_rate_compute_params.h"

namespace blue_sky
  {
  namespace wells
    {

    template <typename strategy_type>
    struct mobility_calc_inj
      {
public:
        typedef strategy_type                               strategy_t;
        typedef typename strategy_t::item_t                 item_t;
        typedef typename strategy_t::index_t                index_t;
        typedef typename strategy_t::item_array_t           item_array_t;

        typedef calc_model <strategy_t>                     calc_model_t;
        typedef typename calc_model_t::data_t               data_t;
        typedef jacobian_matrix <strategy_t>                jmatrix_t;
        typedef well <strategy_t>                           well_t;
        typedef well_controller <strategy_t>                well_controller_t;
        typedef compute_params <strategy_t>                 params_t;

        typedef typename calc_model_t::data_array_t         data_array_t;
        typedef typename calc_model_t::phase_d_t            phase_d_t;
        typedef typename calc_model_t::main_var_array_t     main_var_array_t;

        typedef typename wells::type_helper <strategy_t>    type_helper_t;
        typedef typename type_helper_t::item_rhs_block_t    item_rhs_block_t;

        typedef smart_ptr <calc_model_t, true>              sp_calc_model_t;
        typedef smart_ptr <jmatrix_t, true>                 sp_jmatrix_t;
        typedef smart_ptr <well_t, true>                    sp_well_t;
        typedef smart_ptr <well_controller_t, true>         sp_well_controller_t;

        enum { mult = 1, };
        enum { is_inj = 1, };

public:

        BS_FORCE_INLINE item_t
        get_mult () const
        {
          return 1.0;
        }

        static BS_FORCE_INLINE bool
        is_oil_injection (const params_t &params) 
        {
          return params.inj_params_.is_oil_injection;
        }
        static BS_FORCE_INLINE bool
        is_water_injection (const params_t &params) 
        {
          return params.inj_params_.is_water_injection;
        }
        static BS_FORCE_INLINE bool
        is_gas_injection (const params_t &params) 
        {
          return params.inj_params_.is_gas_injection;
        }
        static BS_FORCE_INLINE bool
        is_o_ctrl (const params_t &params) 
        {
          return params.inj_params_.is_o_ctrl;
        }
        static BS_FORCE_INLINE bool
        is_w_ctrl (const params_t &params) 
        {
          return params.inj_params_.is_w_ctrl;
        }
        static BS_FORCE_INLINE bool
        is_g_ctrl (const params_t &params) 
        {
          return params.inj_params_.is_g_ctrl;
        }

        BS_FORCE_INLINE item_t
        get_oil_mobility (const data_t &data, const params_t &params) const
          {
            const item_t &krp_tetap      = params.inj_params_.krp_tetap;
            const item_t &bo             = INVERS_FVF (data, params.phase_d, FI_PHASE_OIL);
            item_t mo                    = params.main_var == FI_SG_VAR ? bo * krp_tetap : bo * params.inj_params_.krow_tetaow;
            return mo;
          }
        BS_FORCE_INLINE item_t
        get_water_mobility (const data_t &data, const params_t &params) const
          {
            item_t krp_tetap      = params.inj_params_.krp_tetap;
            item_t bw             = INVERS_FVF (data, params.phase_d, FI_PHASE_WATER);
            item_t mw             = params.main_var == FI_SG_VAR ? bw * krp_tetap : bw * params.inj_params_.krow_tetaow;
            return mw;
          }
        BS_FORCE_INLINE item_t
        get_gas_mobility (const data_t &data, const params_t &params) const
          {
            const item_t &krp_tetap   = params.inj_params_.krp_tetap;
            //const item_t &kow_tetaow  = params.inj_params_.krow_tetaow;
            const item_t &bg          = INVERS_FVF (data, params.phase_d, FI_PHASE_GAS);
            item_t mg                 = params.main_var == FI_SG_VAR ? bg * krp_tetap : bg * params.inj_params_.krow_tetaow;
            return mg;
          }

        BS_FORCE_INLINE item_t
        get_mo_po_deriv (const data_t &data, const params_t &params) const
          {
            const item_t &dbo_dpo        = P_DERIV_INVERS_FVF (data, params.phase_d, FI_PHASE_OIL);
            const item_t &krp_tetap      = params.inj_params_.krp_tetap;
            const item_t &krp_dtetap_dpp = params.inj_params_.krp_dtetap_dpp;
            const item_t &bo             = INVERS_FVF (data, params.phase_d, FI_PHASE_OIL);

            return dbo_dpo * krp_tetap + bo * krp_dtetap_dpp;
          }

        BS_FORCE_INLINE item_t
        get_mo_sw_deriv (const data_t &data, const params_t &params) const
          {
            const item_t &bo = INVERS_FVF (data, params.phase_d, FI_PHASE_OIL);
            const item_t &x  = params.inj_params_.dkrow_dsw_tetaow_krw_dtetaw_dsw;

            return bo * x;
          }

        BS_FORCE_INLINE item_t
        get_mo_so_deriv (const data_t &data, const params_t &params) const
          {
            const item_t &dbo_dso        = 0;
            const item_t &krp_tetap      = params.inj_params_.krp_tetap;
            const item_t &bo             = INVERS_FVF (data, params.phase_d, FI_PHASE_OIL);
            const item_t &dkro_dso       = S_DERIV_RELATIVE_PERM (data, params.phase_d, params.n_phases, FI_PHASE_OIL, FI_PHASE_OIL);

            return dbo_dso * krp_tetap + bo * dkro_dso;
          }

        BS_FORCE_INLINE item_t
        get_mo_sg_deriv (const data_t &data, const params_t &params) const
          {
            if (params.main_var == FI_SG_VAR)
              {
                return get_mo_sg_deriv_sg_var (data, params);
              }
            else
              {
                return get_mo_sg_deriv_ro_var (data, params);
              }
          }

        BS_FORCE_INLINE item_t
        get_mo_sg_deriv_sg_var (const data_t &data, const params_t &params) const
        {
          item_t bo         = INVERS_FVF (data, params.phase_d, FI_PHASE_OIL);
          item_t x          = params.inj_params_.dkrog_dsg_tetaog_krg_dtetag_dsg;

          return bo * x;
        }
        BS_FORCE_INLINE item_t
        get_mo_sg_deriv_ro_var (const data_t &data, const params_t &params) const
        {
          item_t kro        = RELATIVE_PERM (data, params.phase_d, FI_PHASE_OIL);
          item_t bo         = INVERS_FVF (data, params.phase_d, FI_PHASE_OIL);
          item_t dtetao_dro = GOR_DERIV_INVERS_VISCOSITY (data);
          item_t tetao      = INVERS_VISCOSITY (data, params.phase_d, FI_PHASE_OIL);
          item_t dbo_dro    = GOR_DERIV_INVERS_FVF (data);

          return kro * (bo * dtetao_dro + tetao * dbo_dro);
        }

        BS_FORCE_INLINE item_t
        get_mw_po_deriv (const data_t &data, const params_t &params) const
          {
            return (RELATIVE_PERM_O * P_DERIV_INVERS_VISCOSITY_O +
                    RELATIVE_PERM_W * P_DERIV_INVERS_VISCOSITY_W) -
                   (INVERS_FVF_W * params.inj_params_.krp_tetap * P_DERIV_INVERS_FVF_W);

            //const item_t &dbw_dpw        = P_DERIV_INVERS_FVF (data, params.phase_d, FI_PHASE_WATER);
            //const item_t &krp_tetap      = params.krp_tetap;
            //const item_t &krp_dtetap_dpp = params.krp_dtetap_dpp;
            //const item_t &bw             = INVERS_FVF (data, params.phase_d, FI_PHASE_WATER);

            //return (dbw_dpw * krp_tetap + bw * krp_dtetap_dpp);
          }
        BS_FORCE_INLINE item_t
        get_mw_sw_deriv (const data_t &data, const params_t &params) const
          {
            return (S_DERIV_RELATIVE_PERM_WW * INVERS_VISCOSITY_W +
                    S_DERIV_RELATIVE_PERM_OW * INVERS_VISCOSITY_O +
                    RELATIVE_PERM_W * S_DERIV_INVERS_VISCOSITY_W) -
                   (INVERS_FVF_W * params.inj_params_.krp_tetap * P_DERIV_INVERS_FVF_W);

            //const item_t &dbw_dsw   = S_DERIV_INVERS_FVF (data, params.phase_d, FI_PHASE_WATER);
            //const item_t &krp_tetap = params.krp_tetap;
            //const item_t &bw        = INVERS_FVF (data, params.phase_d, FI_PHASE_WATER);
            //const item_t &x         = params.dkrow_dsw_tetaow_krw_dtetaw_dsw;

            //return (dbw_dsw * krp_tetap + bw * x);
          }
        BS_FORCE_INLINE item_t
        get_mw_so_deriv (const data_t &data, const params_t &params) const
          {
            return S_DERIV_RELATIVE_PERM_OO * INVERS_VISCOSITY_O;

            //const item_t &bo        = INVERS_FVF (data, params.phase_d, FI_PHASE_OIL);
            //const item_t &dkrw_dso  = S_DERIV_RELATIVE_PERM (data, params.phase_d, params.n_phases, FI_PHASE_OIL, FI_PHASE_OIL);
            //return bo * dkrw_dso;
          }
        BS_FORCE_INLINE item_t
        get_mw_sg_deriv (const data_t &data, const params_t &params) const
          {
            const item_t &bw = INVERS_FVF (data, params.phase_d, FI_PHASE_WATER);

            if (params.main_var == FI_SG_VAR)
              {
                item_t x = params.inj_params_.dkrog_dsg_tetaog_krg_dtetag_dsg;
                return bw * x;
              }
            else
              {
                item_t kro        = RELATIVE_PERM (data, params.phase_d, FI_PHASE_OIL);
                item_t dtetao_dro = GOR_DERIV_INVERS_VISCOSITY (data);

                return bw * kro * dtetao_dro;
              }
          }

        BS_FORCE_INLINE item_t
        get_mg_po_deriv (const data_t &data, const params_t &params) const
          {
            const item_t &dbg_dpg        = P_DERIV_INVERS_FVF (data, params.phase_d, FI_PHASE_GAS);
            const item_t &krp_tetap      = params.krp_tetap;
            const item_t &krp_dtetap_dpp = params.krp_dtetap_dpp;
            const item_t &bg             = INVERS_FVF (data, params.phase_d, FI_PHASE_GAS);

            return dbg_dpg * krp_tetap + bg * krp_dtetap_dpp;
          }
        BS_FORCE_INLINE item_t
        get_gor_po_deriv (const data_t &data, const params_t &params) const
          {
            return 0;
          }
        BS_FORCE_INLINE item_t
        get_gor (const data_t &data, const params_t &params) const
          {
            return 0;
          }

        BS_FORCE_INLINE item_t
        get_mg_sw_deriv (const data_t &data, const params_t &params) const
          {
            const item_t &bg = INVERS_FVF (data, params.phase_d, FI_PHASE_GAS);
            const item_t &x  = params.inj_params_.dkrow_dsw_tetaow_krw_dtetaw_dsw;

            return bg * x;
          }
        BS_FORCE_INLINE item_t
        get_mg_so_deriv (const data_t &data, const params_t &params) const
          {
            const item_t &bo        = INVERS_FVF (data, params.phase_d, FI_PHASE_OIL);
            const item_t &dkrg_dso  = S_DERIV_RELATIVE_PERM (data, params.phase_d, params.n_phases, FI_PHASE_GAS, FI_PHASE_OIL);

            return bo * dkrg_dso;
          }
        BS_FORCE_INLINE item_t
        get_gas_sw_deriv (const data_t &data, const params_t &params) const
        {
          return get_mg_sw_deriv (data, params) * params.Pg;
        }
        BS_FORCE_INLINE item_t 
        get_gas_so_deriv (const data_t &data, const params_t &params) const
        {
          return get_mg_so_deriv (data, params) * params.Pg;
        }
        BS_FORCE_INLINE item_t
        get_gas_po_deriv (const data_t &data, const params_t &params) const
          {
            item_t dbg_dpg        = P_DERIV_INVERS_FVF (data, params.phase_d, FI_PHASE_GAS);
            item_t krp_tetap      = params.inj_params_.krp_tetap;
            item_t krp_dtetap_dpp = params.inj_params_.krp_dtetap_dpp;
            item_t bg             = INVERS_FVF (data, params.phase_d, FI_PHASE_GAS);
            item_t mg             = params.main_var == FI_SG_VAR ? bg * krp_tetap : bg * params.inj_params_.krow_tetaow;
            item_t Pg             = params.Pg;
            item_t gw             = params.gw;

            item_t po_deriv       = -gw * ((dbg_dpg * krp_tetap + bg * krp_dtetap_dpp) * Pg - mg);
            return po_deriv;
          }
        BS_FORCE_INLINE item_t
        get_gas_sg_deriv (const data_t &data, const params_t &params) const
          {
            item_t gw             = params.gw;
            item_t bg             = INVERS_FVF (data, params.phase_d, FI_PHASE_GAS);
            item_t Pg             = params.Pg;

            if (params.main_var == FI_SG_VAR)
              {
                item_t dbg_dsg    = S_DERIV_INVERS_FVF (data, params.phase_d, FI_PHASE_GAS);
                item_t krp_tetap  = params.inj_params_.krp_tetap;
                item_t x          = params.inj_params_.dkrog_dsg_tetaog_krg_dtetag_dsg;
                item_t mg         = bg * krp_tetap;
                item_t dpcgo_dsg  = S_DERIV_CAP_PRESSURE (data, params.phase_d, FI_PHASE_GAS);

                item_t sg_deriv   = -gw * ((dbg_dsg * krp_tetap + bg * x) * Pg - mg * dpcgo_dsg);
                return sg_deriv;
              }
            else
              {
                item_t kro        = RELATIVE_PERM (data, params.phase_d, FI_PHASE_OIL);
                item_t dtetao_dro = GOR_DERIV_INVERS_VISCOSITY (data);

                item_t ro_deriv   = -gw * bg * kro * dtetao_dro * Pg;
                return ro_deriv;
              }
          }
        BS_FORCE_INLINE item_t
        get_gas_pref_deriv (const data_t &data, const params_t &params) const
          {
            const item_t &gw  = params.gw;
            const item_t &mg  = get_gas_mobility (data, params);

            item_t pref_deriv = -gw * mg;
            return pref_deriv;
          }
        BS_FORCE_INLINE item_t
        get_gas_rate (const data_t &data, const params_t &params) const
          {
            const item_t &Pg          = params.Pg;
            const item_t &gw          = params.gw;
            const item_t &mg          = get_gas_mobility (data, params);

            //item_t Rso = GAS_OIL_RATIO (data);
#ifdef _DEBUG
            //BOSOUT (section::wells, level::debug) << boost::format ("[%s : %d] qg (%.20e) gw: %.20e mg: %.20e pg: %.20e gor: %.20e") % params.well_->name ().c_str () % params.n_block % (gw * mg * Pg) % gw % mg % Pg % Rso << bs_end;
#endif

            item_t rate               = gw * mg * Pg;
            return rate;
          }

        BS_FORCE_INLINE item_t
        get_free_gas_rate (const data_t &data, const params_t &params) const
          {
            const item_t &Pg          = params.Pg;
            const item_t &gw          = params.gw;
            const item_t &mg          = get_gas_mobility (data, params);

            item_t rate               = gw * mg * Pg;
            return rate;
          }

        BS_FORCE_INLINE item_t
        get_solution_gas_rate (const data_t &data, const params_t &params) const
        {
          return 0.0;
        }
      };


  } // namespace wells
} // namespace blue_sky

#endif // #ifndef BS_WELLS_WELL_RATE_INJ_MOBILITY_H_

