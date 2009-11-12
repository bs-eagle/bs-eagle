/**
 *       \file  well_rate_control_prod_mobility.h
 *      \brief  Calculates mobility for production wells
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  21.11.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Obsolete, should be removed
 * */
#ifndef BS_WELLS_WELL_RATE_PROD_MOBILITY_H_
#define BS_WELLS_WELL_RATE_PROD_MOBILITY_H_

#include "calc_model.h"
#include "calc_model_data_accessors.h"
#include "well_rate_compute_params.h"

namespace blue_sky
  {
  namespace wells
    {

    template <typename strategy_type>
    struct mobility_calc_prod
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

        enum { mult = -1, };
        enum { is_inj = 0, };

public:

        BS_FORCE_INLINE item_t
        get_mult () const
          {
            return item_t (mult);
          }

        static BS_FORCE_INLINE bool
        is_oil_injection (const params_t &params) 
        {
          return true;
        }
        static BS_FORCE_INLINE bool
        is_water_injection (const params_t &params) 
        {
          return true;
        }
        static BS_FORCE_INLINE bool
        is_gas_injection (const params_t &params) 
        {
          return true;
        }
        static BS_FORCE_INLINE bool
        is_o_ctrl (const params_t &params) 
        {
          return params.prod_params_.is_o_ctrl;
        }
        static BS_FORCE_INLINE bool
        is_w_ctrl (const params_t &params) 
        {
          return params.prod_params_.is_w_ctrl;
        }
        static BS_FORCE_INLINE bool
        is_g_ctrl (const params_t &params) 
        {
          return params.prod_params_.is_g_ctrl;
        }

        BS_FORCE_INLINE item_t
        get_oil_mobility (const data_t &data, const params_t &params) const
        {
          return MOBILITY (data, params.phase_d, FI_PHASE_OIL);
        }
        BS_FORCE_INLINE item_t
        get_water_mobility (const data_t &data, const params_t &params) const
        {
          return MOBILITY (data, params.phase_d, FI_PHASE_WATER);
        }
        BS_FORCE_INLINE item_t
        get_gas_mobility (const data_t &data, const params_t &params) const
        {
          return MOBILITY (data, params.phase_d, FI_PHASE_GAS);
        }

        BS_FORCE_INLINE item_t
        get_mo_po_deriv (const data_t &data, const params_t &params) const
        {
          return P_DERIV_MOBILITY_O;

          //const item_t &kro       = RELATIVE_PERM (data, params.phase_d, FI_PHASE_OIL);
          //const item_t &depso_dpo = P_DERIV_INVERS_VISC_FVF (data, params.phase_d, FI_PHASE_OIL);

          //return kro * depso_dpo;
        }

        BS_FORCE_INLINE item_t
        get_mo_sw_deriv (const data_t &data, const params_t &params) const
        {
          return S_DERIV_MOBILITY_OW;

          //const item_t &dkro_dsw   = S_DERIV_RELATIVE_PERM (data, params.phase_d, params.n_phases, FI_PHASE_OIL, FI_PHASE_WATER);
          //const item_t &epso       = INVERS_VISC_FVF (data, params.phase_d, FI_PHASE_OIL);

          //return dkro_dsw * epso;
        }

        BS_FORCE_INLINE item_t
        get_mo_so_deriv (const data_t &data, const params_t &params) const
        {
          return S_DERIV_MOBILITY_OO;
        }

        BS_FORCE_INLINE item_t
        get_mo_sg_deriv (const data_t &data, const params_t &params) const
        {
          //return S_DERIV_MOBILITY_OG;

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
          return S_DERIV_MOBILITY_OG;

          item_t dkro_dsg   = S_DERIV_RELATIVE_PERM (data, params.phase_d, params.n_phases, FI_PHASE_OIL, FI_PHASE_GAS);
          item_t epso       = INVERS_VISC_FVF (data, params.phase_d, FI_PHASE_OIL);

          return dkro_dsg * epso;
        }
        BS_FORCE_INLINE item_t
        get_mo_sg_deriv_ro_var (const data_t &data, const params_t &params) const
        {
          item_t kro        = RELATIVE_PERM (data, params.phase_d, FI_PHASE_OIL);
          item_t depso_dro  = GOR_DERIV_INVERS_VISC_FVF (data);

          return kro * depso_dro;
        }

        BS_FORCE_INLINE item_t
        get_mw_po_deriv (const data_t &data, const params_t &params) const
        {
          return P_DERIV_MOBILITY_W;

          //item_t krw        = RELATIVE_PERM (data, params.phase_d, FI_PHASE_WATER);
          //item_t depsw_dpo  = P_DERIV_INVERS_VISC_FVF (data, params.phase_d, FI_PHASE_WATER);

          //return krw * depsw_dpo;
        }
        BS_FORCE_INLINE item_t
        get_mw_sw_deriv (const data_t &data, const params_t &params) const
        {
          return S_DERIV_MOBILITY_WW;

          //item_t dkrw_dsw   = S_DERIV_RELATIVE_PERM (data, params.phase_d, params.n_phases, FI_PHASE_WATER, FI_PHASE_WATER);
          //item_t epsw       = INVERS_VISC_FVF (data, params.phase_d, FI_PHASE_WATER);
          //item_t krw        = RELATIVE_PERM (data, params.phase_d, FI_PHASE_WATER);
          //item_t dbw_dsw    = S_DERIV_INVERS_FVF (data, params.phase_d, FI_PHASE_WATER);
          //item_t tetaw      = INVERS_VISCOSITY (data, params.phase_d, FI_PHASE_WATER);
          //item_t bw         = INVERS_FVF (data, params.phase_d, FI_PHASE_WATER);
          //item_t dtetaw_dsw = S_DERIV_INVERS_VISCOSITY (data, params.phase_d, FI_PHASE_WATER);

          //item_t sw_deriv   = (dkrw_dsw * epsw + krw * (dbw_dsw * tetaw + bw * dtetaw_dsw));
          //return sw_deriv;
        }
        BS_FORCE_INLINE item_t
        get_mw_so_deriv (const data_t &data, const params_t &params) const
        {
          return S_DERIV_MOBILITY_WO;
        }
        BS_FORCE_INLINE item_t
        get_mw_sg_deriv (const data_t &data, const params_t &params) const
        {
          return S_DERIV_MOBILITY_WG;
        }

        //BS_FORCE_INLINE item_t
        //get_mg_sw_deriv (const data_t &data, const params_t &params) const
        //{
        //  const item_t &Rso      = GAS_OIL_RATIO (data);
        //  const item_t &dmo_dsw  = S_DERIV_MOBILITY_OW;

        //  return Rso * dmo_dsw;
        //}
        //BS_FORCE_INLINE item_t
        //get_mg_so_deriv (const data_t &data, const params_t &params) const
        //{
        //  return S_DERIV_MOBILITY_GO + S_DERIV_MOBILITY_OO * GAS_OIL_RATIO (data);
        //}

        BS_FORCE_INLINE item_t
        get_gas_sw_deriv (const data_t &data, const params_t &params) const
        {
          return S_DERIV_MOBILITY_GW * params.Pg + S_DERIV_MOBILITY_OW * GAS_OIL_RATIO (data) * params.Po;
        }

        BS_FORCE_INLINE item_t
        get_gas_so_deriv (const data_t &data, const params_t &params) const
        {
          return S_DERIV_MOBILITY_GO * params.Pg + S_DERIV_MOBILITY_OO * GAS_OIL_RATIO (data) * params.Po;
        }

        BS_FORCE_INLINE item_t
        get_gas_po_deriv (const data_t &data, const params_t &params) const
        {
          item_t gw             = params.gw;
          item_t mo             = MOBILITY (data, params.phase_d, FI_PHASE_OIL);
          item_t Po             = params.Po;
          item_t dmo_dpo        = P_DERIV_MOBILITY_O;
          item_t Rso            = GAS_OIL_RATIO (data);
          item_t ro_deriv       = Rso * (dmo_dpo * Po - mo);

          item_t po_deriv       = 0;
          if (params.main_var == FI_SG_VAR)
            {
              item_t dmg_dpo    = P_DERIV_MOBILITY_G;
              item_t Pg         = params.Pg;
              item_t mg         = MOBILITY (data, params.phase_d, FI_PHASE_GAS);
              item_t dRso_dpo   = P_DERIV_GAS_OIL_RATIO (data);

              po_deriv = -gw * (dmg_dpo * Pg - mg + dRso_dpo * mo * Po + ro_deriv);
            }
          else
            {
              po_deriv = -gw * ro_deriv;
            }
          return po_deriv;
        }
        BS_FORCE_INLINE item_t
        get_gas_sg_deriv (const data_t &data, const params_t &params) const
        {
          if (params.main_var == FI_SG_VAR)
            {
              item_t gw         = params.gw;
              item_t Pg         = params.Pg;
              item_t Po         = params.Po;

              //item_t dkrg_dsg   = S_DERIV_RELATIVE_PERM_GG;
              //item_t dkro_dsg   = S_DERIV_RELATIVE_PERM_OG;

              //item_t epsg       = INVERS_VISC_FVF_G;
              //item_t epso       = INVERS_VISC_FVF_O;

              item_t dpcgo_dsg  = S_DERIV_CAP_PRESSURE_G;
              item_t Rso        = GAS_OIL_RATIO (data);
              item_t mg         = MOBILITY_G;
              item_t dmg_dsg    = S_DERIV_MOBILITY_GG;
              item_t dmo_dsg    = S_DERIV_MOBILITY_OG;

              return -gw * (dmg_dsg * Pg - mg * dpcgo_dsg + Rso * dmo_dsg * Po);
            }
          else
            {
              item_t gw         = params.gw;
              item_t Po         = params.Po;
              item_t Ro         = GAS_OIL_RATIO (data);
              item_t mo         = MOBILITY_O;
              item_t dmo_dsg    = S_DERIV_MOBILITY_OG;

              item_t ro_deriv   = -gw * (mo + Ro * dmo_dsg) * Po;
              return ro_deriv;
            }
        }
        BS_FORCE_INLINE item_t
        get_gas_pref_deriv (const data_t &data, const params_t &params) const
        {
          if (params.main_var == FI_SG_VAR)
            {
              item_t gw     = params.gw;
              item_t mg     = MOBILITY (data, params.phase_d, FI_PHASE_GAS);
              item_t mo     = MOBILITY (data, params.phase_d, FI_PHASE_OIL);
              item_t Rso    = GAS_OIL_RATIO (data);

              item_t pref_deriv = -gw * (mg + Rso * mo);
              return pref_deriv;
            }
          else
            {
              item_t gw     = params.gw;
              item_t Ro     = GAS_OIL_RATIO (data);
              item_t mo     = MOBILITY (data, params.phase_d, FI_PHASE_OIL);

              item_t pref_deriv = -gw * Ro * mo;
              return pref_deriv;
            }
        }
        BS_FORCE_INLINE item_t
        get_gas_rate (const data_t &data, const params_t &params) const
        {
          const item_t &gw  = params.gw;
          const item_t &Pg  = params.Pg;
          const item_t &Po  = params.Po;
          const item_t &mg  = MOBILITY (data, params.phase_d, FI_PHASE_GAS);
          const item_t &mo  = MOBILITY (data, params.phase_d, FI_PHASE_OIL);
          const item_t &Rso = GAS_OIL_RATIO (data);

#ifdef _DEBUG
          BOSOUT (section::wells, level::debug) << boost::format ("[%s : %d] qg (%.20e) gw: %.20e mg: %.20e pg: %.20e gor: %.20e") % params.well_->name ().c_str () % params.n_block % (gw * (mg * Pg + Rso * mo * Po)) % gw % mg % Pg % Rso << bs_end;
#endif

          return gw * (mg * Pg + Rso * mo * Po);
        }

        BS_FORCE_INLINE item_t
        get_free_gas_rate (const data_t &data, const params_t &params) const
        {
          const item_t &gw  = params.gw;
          const item_t &Pg  = params.Pg;
          const item_t &mg  = MOBILITY (data, params.phase_d, FI_PHASE_GAS);

          return gw * mg * Pg;
        }

        BS_FORCE_INLINE item_t
        get_solution_gas_rate (const data_t &data, const params_t &params) const
        {
          const item_t &gw  = params.gw;
          const item_t &Po  = params.Po;
          const item_t &mo  = MOBILITY (data, params.phase_d, FI_PHASE_OIL);
          const item_t &Rso = GAS_OIL_RATIO (data);

          return gw * Rso * mo * Po;
        }
      };

  } // namespace wells
} // namespace blue_sky

#endif  // #ifndef BS_WELLS_WELL_RATE_PROD_MOBILITY_H_

