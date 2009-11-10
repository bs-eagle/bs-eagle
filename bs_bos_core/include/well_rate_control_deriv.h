/**
 * \file well_rate_control_deriv.h
 * \brief
 * \author Sergey Miryanov
 * \date 21.11.2008
 * */
#ifndef BS_WELLS_WELL_RATE_CONTROL_DERIV_H_
#define BS_WELLS_WELL_RATE_CONTROL_DERIV_H_

#include "apply_wefac.h"

namespace blue_sky
  {
  namespace wells
    {

    template <typename mobility_calc_t>
    struct compute_deriv
      {
        typedef typename mobility_calc_t::strategy_t          strategy_t;
        typedef mobility_calc_t                               mobility_t;
        typedef typename strategy_t::item_t                   item_t;
        typedef typename strategy_t::rhs_item_t               rhs_item_t;
        typedef calc_model <strategy_t>                       calc_model_t;
        typedef well <strategy_t>                             well_t;
        typedef typename calc_model_t::data_t                 data_t;
        typedef typename mobility_calc_t::params_t            params_t;
        typedef typename wells::type_helper <strategy_t>      type_helper_t;
        typedef typename type_helper_t::item_rr_block_t       item_rr_block_t;
        typedef typename type_helper_t::item_rw_block_t       item_rw_block_t;
        typedef typename type_helper_t::item_wr_block_t       item_wr_block_t;
        typedef typename type_helper_t::item_q_rate_t         item_q_rate_t;
        typedef typename type_helper_t::item_rhs_block_t      item_rhs_block_t;
        typedef typename type_helper_t::item_ps_block_t       item_ps_block_t;

        typedef typename well_t::sp_connection_t              sp_connection_t;

protected:

        compute_deriv (const mobility_calc_t &mobility_calc)
            : mobility_calc_ (mobility_calc)
            , mult (mobility_calc_t::mult)
        {
        }

        BS_FORCE_INLINE item_t
        compute_oil_rate (const data_t &data, params_t &params) const
          {
            const item_t &mo  = mobility_calc_.get_oil_mobility (data, params);
            const item_t &Po  = params.Po;
            const item_t &gw  = params.gw;

#ifdef _DEBUG
            BOSOUT (section::wells, level::debug) << boost::format ("[%s : %d] qo (%.20e) gw: %.20e mo: %.20e po: %.20e") % params.well_->name ().c_str () % params.n_block % (gw * mo * Po) % gw % mo % Po << bs_end;
#endif

            item_t rate       = gw * mo * Po;
            return rate;
          }
        BS_FORCE_INLINE item_t
        compute_oil_po_deriv (const data_t &data, params_t &params) const
          {
            const item_t &dmo_dpo = mobility_calc_.get_mo_po_deriv (data, params);
            const item_t &mo      = mobility_calc_.get_oil_mobility (data, params);
            item_t gw             = params.gw;
            item_t Po             = params.Po;

#ifdef _DEBUG
            //BOSOUT (section::wells, level::debug) << boost::format ("[%s : %d] oil_po: %.20e") % params.well_->name ().c_str () % params.n_block % (-gw * (dmo_dpo * Po - mo)) << bs_end;
#endif

            item_t po_deriv       = -gw * (dmo_dpo * Po - mo);
            return po_deriv;
          }

        BS_FORCE_INLINE item_t
        compute_oil_sw_deriv (const data_t &data, params_t &params) const
          {
            const item_t &dmo_dsw = mobility_calc_.get_mo_sw_deriv (data, params);
            item_t gw         = params.gw;
            item_t Po         = params.Po;

#ifdef _DEBUG
            //BOSOUT (section::wells, level::debug) << boost::format ("[%s : %d] oil_sw: %.20e") % params.well_->name ().c_str () % params.n_block % (-gw * dmo_dsw * Po) << bs_end;
#endif

            item_t sw_deriv   = -gw * dmo_dsw * Po;
            return sw_deriv;
          }
        BS_FORCE_INLINE item_t
        compute_oil_so_deriv (const data_t &data, params_t &params) const
          {
            const item_t &dmo_dso = mobility_calc_.get_mo_so_deriv (data, params);
            item_t gw             = params.gw;
            item_t Po             = params.Po;

#ifdef _DEBUG
            //BOSOUT (section::wells, level::debug) << boost::format ("[%s : %d] oil_so: %.20e") % params.well_->name ().c_str () % params.n_block % (-gw * dmo_dso * Po) << bs_end;
#endif

            item_t so_deriv       = -gw * dmo_dso * Po;
            return so_deriv;
          }
        BS_FORCE_INLINE item_t
        compute_oil_sg_deriv (const data_t &data, params_t &params) const
          {
            const item_t &dmo_dsg = mobility_calc_.get_mo_sg_deriv (data, params);
            item_t gw             = params.gw;
            item_t Po             = params.Po;

#ifdef _DEBUG
            //BOSOUT (section::wells, level::debug) << boost::format ("[%s : %d] oil_sg: %.20e") % params.well_->name ().c_str () % params.n_block % (-gw * dmo_dsg * Po) << bs_end;
#endif

            item_t sg_deriv       = -gw * dmo_dsg * Po;
            return sg_deriv;
          }
        BS_FORCE_INLINE item_t
        compute_oil_pref_deriv (const data_t &data, params_t &params) const
          {
            const item_t &mo  = mobility_calc_.get_oil_mobility (data, params);
            const item_t &gw  = params.gw;

#ifdef _DEBUG
            //BOSOUT (section::wells, level::debug) << boost::format ("[%s : %d] oil_pref: %.20e") % params.well_->name ().c_str () % params.n_block % (-gw * mo) << bs_end;
#endif

            item_t pref_deriv = -gw * mo;
            return pref_deriv;
          }

        BS_FORCE_INLINE item_t
        compute_water_rate (const data_t &data, params_t &params) const
          {
            const item_t &mw  = mobility_calc_.get_water_mobility (data, params);
            const item_t &Pw  = params.Pw;
            const item_t &gw  = params.gw;

#ifdef _DEBUG
            BOSOUT (section::wells, level::debug) << boost::format ("[%s : %d] qw (%.20e) gw: %.20e mw: %.20e pw: %.20e") % params.well_->name ().c_str () % params.n_block % (gw * mw * Pw) % gw % mw % Pw << bs_end;
#endif

            item_t rate       = gw * mw * Pw;
            return rate;
          }
        BS_FORCE_INLINE item_t
        compute_water_po_deriv (const data_t &data, params_t &params) const
          {
            const item_t &dmw_dpo = mobility_calc_.get_mw_po_deriv (data, params);
            const item_t &mw      = mobility_calc_.get_water_mobility (data, params);
            item_t Pw             = params.Pw;
            item_t gw             = params.gw;

#ifdef _DEBUG
            //BOSOUT (section::wells, level::debug) << boost::format ("[%s : %d] wat_po: %.20e") % params.well_->name ().c_str () % params.n_block % (-gw * (dmw_dpo * Pw - mw)) << bs_end;
#endif

            item_t po_deriv       = -gw * (dmw_dpo * Pw - mw);
            return po_deriv;
          }
        BS_FORCE_INLINE item_t
        compute_water_sw_deriv (const data_t &data, params_t &params) const
          {
            const item_t &dmw_dsw = mobility_calc_.get_mw_sw_deriv (data, params);
            const item_t &mw      = mobility_calc_.get_water_mobility (data, params);
            item_t gw             = params.gw;
            item_t Pw             = params.Pw;
            item_t dpcwo_dsw      = S_DERIV_CAP_PRESSURE (data, params.phase_d, FI_PHASE_WATER);

#ifdef _DEBUG
            //BOSOUT (section::wells, level::debug) << boost::format ("[%s : %d] wat_sw: %.20e") % params.well_->name ().c_str () % params.n_block % (-gw * (dmw_dsw * Pw - mw * dpcwo_dsw)) << bs_end;
#endif

            item_t sw_deriv       = -gw * (dmw_dsw * Pw - mw * dpcwo_dsw);
            return sw_deriv;
          }
        BS_FORCE_INLINE item_t
        compute_water_so_deriv (const data_t &data, params_t &params) const
          {
            const item_t &dmw_dso = mobility_calc_.get_mw_so_deriv (data, params);
            item_t gw             = params.gw;
            item_t Pw             = params.Pw;

#ifdef _DEBUG
            //BOSOUT (section::wells, leve::debug) << boost::format ("[%s : %d] wat_so: %.20e") % params.well_->name ().c_str () % params.n_block % (-gw * dmw_dso * Pw) << bs_end;
#endif

            item_t so_deriv       = -gw * dmw_dso * Pw;
            return so_deriv;
          }
        BS_FORCE_INLINE item_t
        compute_water_sg_deriv (const data_t &data, params_t &params) const
          {
            const item_t &dmw_dsg = mobility_calc_.get_mw_sg_deriv (data, params);
            item_t gw             = params.gw;
            item_t Pw             = params.Pw;

#ifdef _DEBUG
            //BOSOUT (section::wells, level::debug) << boost::format ("[%s : %d] wat_sg: %.20e") % params.well_->name ().c_str () % params.n_block % (-gw * dmw_dsg * Pw) << bs_end;
#endif

            item_t sg_deriv       = -gw * dmw_dsg * Pw;
            return sg_deriv;
          }
        BS_FORCE_INLINE item_t
        compute_water_pref_deriv (const data_t &data, params_t &params) const
          {
            const item_t &mw  = mobility_calc_.get_water_mobility (data, params);
            const item_t &gw  = params.gw;

#ifdef _DEBUG
            //BOSOUT (section::wells, level::debug) << boost::format ("[%s : %d] wat_pref: %.20e") % params.well_->name ().c_str () % params.n_block % (-gw * mw) << bs_end;
#endif

            item_t pref_deriv = -gw * mw;
            return pref_deriv;
          }

        BS_FORCE_INLINE item_t
        compute_gas_rate (const data_t &data, params_t &params) const
          {
            return mobility_calc_.get_gas_rate (data, params);
          }

        BS_FORCE_INLINE item_t
        compute_free_gas_rate (const data_t &data, params_t &params) const
        {
          return mobility_calc_.get_free_gas_rate (data, params);
        }

        BS_FORCE_INLINE item_t
        compute_solution_gas_rate (const data_t &data, params_t &params) const
        {
          return mobility_calc_.get_solution_gas_rate (data, params);
        }

        BS_FORCE_INLINE item_t
        compute_gas_po_deriv (const data_t &data, params_t &params) const
          {
#ifdef _DEBUG
            //BOSOUT (section::wells, level::debug) << boost::format ("[%s : %d] gas_po: %.20e") % params.well_->name ().c_str () % params.n_block % mobility_calc_.get_gas_po_deriv (data, params) << bs_end;
#endif

            return mobility_calc_.get_gas_po_deriv (data, params);
          }
        BS_FORCE_INLINE item_t
        compute_gas_sw_deriv (const data_t &data, params_t &params) const
          {
            const item_t &dmg_dsw = mobility_calc_.get_gas_sw_deriv (data, params);//get_mg_sw_deriv (data, params);

            item_t gw         = params.gw;
            //item_t Pg         = params.Pg;

#ifdef _DEBUG
            //BOSOUT (section::wells, level::debug) << boost::format ("[%s : %d] gas_sw: %.20e") % params.well_->name ().c_str () % params.n_block % (-gw * dmg_dsw * Pg) << bs_end;
#endif

            item_t sw_deriv   = -gw * dmg_dsw;// * Pg;
            return sw_deriv;
          }
        BS_FORCE_INLINE item_t
        compute_gas_so_deriv (const data_t &data, params_t &params) const
          {
            const item_t &dmg_dso = mobility_calc_.get_gas_so_deriv (data, params);//get_mg_so_deriv (data, params);
            item_t gw             = params.gw;
            //item_t Pg             = params.Pg;

#ifdef _DEBUG
            //BOSOUT (section::wells, level::debug) << boost::format ("[%s : %d] gas_so: %.20e") % params.well_->name ().c_str () % params.n_block % (-gw * dmg_dso * Pg) << bs_end;
#endif

            item_t so_deriv       = -gw * dmg_dso;// * Pg;
            return so_deriv;
          }
        BS_FORCE_INLINE item_t
        compute_gas_sg_deriv (const data_t &data, params_t &params) const
          {
#ifdef _DEBUG
            //BOSOUT (section::wells, level::debug) << boost::format ("[%s : %d] gas_sg: %.20e") % params.well_->name ().c_str () % params.n_block % mobility_calc_.get_gas_sg_deriv (data, params) << bs_end;
#endif

            return mobility_calc_.get_gas_sg_deriv (data, params);
          }
        BS_FORCE_INLINE item_t
        compute_gas_pref_deriv (const data_t &data, params_t &params) const
          {
#ifdef _DEBUG
            //BOSOUT (section::wells, level::debug) << boost::format ("[%s : %d] gas_pref: %.20e") % params.well_->name ().c_str () % params.n_block % mobility_calc_.get_gas_pref_deriv (data, params) << bs_end;
#endif

            return mobility_calc_.get_gas_pref_deriv (data, params);
          }

        const mobility_calc_t &mobility_calc_;
        item_t                mult;
      };

    template <typename mobility_calc_t>
    struct dummy_deriv
      {
        typedef typename mobility_calc_t::strategy_t  strategy_t;
        typedef mobility_calc_t                       mobility_t;
        typedef typename strategy_t::item_t           item_t;
        typedef calc_model <strategy_t>               calc_model_t;
        typedef typename calc_model_t::data_t         data_t;
        typedef typename mobility_calc_t::params_t    params_t;
        typedef well <strategy_t>                     well_t;

        typedef typename well_t::sp_connection_t      sp_connection_t;

public:

        dummy_deriv (const mobility_calc_t &)
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
        }
        void
        update_wr (const sp_connection_t &c, item_t ww) const
        {
        }
        void
        update_rr (const sp_connection_t &c, const data_t &data, params_t &params) const
        {
        }
        void
        update_rhs_flux (const sp_connection_t &c, const data_t &data, params_t &params) const
        {
        }
        void
        compute_bw_value (params_t &params) const
        {
        }
        void
        apply_wefac (const sp_connection_t &c, params_t &params) const
        {
        }
        void
        update_rate (const sp_connection_t &c, params_t &params) const
        {
        }
      };

  } // namespace wells
} // namespace blue_sky

#endif  // #ifndef BS_WELLS_WELL_RATE_CONTROL_DERIV_H_

