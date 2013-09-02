/**
 *       \file  fip_results_storage.h
 *      \brief  Storage for fip results
 *     \author  Ilshat Sayfullin
 *       \date  10.02.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Describe
 * */

#ifndef __BS_FIP_RESULTS_STORAGE_H
#define __BS_FIP_RESULTS_STORAGE_H

#include "calc_model.h"
#include "well_type_helper.h"

namespace blue_sky
  {

  enum fip_d_params
  {
    FIP_D_PARAM_COIP = 0,
    FIP_D_PARAM_CMOIP,
    FIP_D_PARAM_CWIP,
    FIP_D_PARAM_CGIP,
    FIP_D_PARAM_HPROD,
    FIP_D_PARAM_PAV,
    FIP_D_PARAM_CPVOL,

    FIP_D_PARAM_OIL_IN,
    FIP_D_PARAM_WATER_IN,
    FIP_D_PARAM_WATER_OUT,

    FIP_D_PARAM_SWAV,
    FIP_D_PARAM_SOAV,
    FIP_D_PARAM_SGAV,
    FIP_D_PARAM_GOR,

    FIP_D_PARAM_COUNT
  };

  class BS_API_PLUGIN fip_results
    {
      // methods
    public:
      // default constructor
      fip_results ();
      // default destructor
      ~fip_results ();

      // variables
    public:
      typedef std::vector<double> dates_type;
      typedef std::vector<float> d_params_internal_type;
      typedef std::vector<d_params_internal_type> d_params_type;

      dates_type dates;                           //!< dates array
      d_params_type d_params;                     //!< set of float parameter arrays
    };


  class BS_API_PLUGIN fip_results_storage : public objbase
    {
      BLUE_SKY_TYPE_DECL (fip_results_storage);
      // methods
    public:
      ~fip_results_storage ();
      // get summary fip data for all time steps for region
      int get_summary_fip_data_for_region ();

      // variables
    public:
      typedef std::map <int, fip_results> fip_type;
      typedef std::pair <int, fip_results> fip_pair_type;
      fip_type fip;

    };


  struct save_fip_data
    {
      typedef t_long                                  index_t;

      typedef calc_model::data_array_t                data_array_t;
      typedef calc_model::item_t                      item_t;
      typedef wells::type_helper::item_rr_block_t     item_rr_block_t;
      typedef wells::type_helper::item_rw_block_t     item_rw_block_t;
      typedef wells::type_helper::item_wr_block_t     item_wr_block_t;
      typedef wells::connection                       connection_t;
      typedef well                                    well_t;
      typedef facility_manager::well_iterator_t       well_iterator_t;

      typedef smart_ptr <calc_model, true>            sp_calc_model_t;
      typedef smart_ptr <well_t, true>                sp_well_t;
      typedef smart_ptr <connection_t, true>          sp_connection_t;
      typedef smart_ptr<fip_results_storage, true>    sp_fip_results_storage;

//! copy current fip data to storage
      void
      copy_fip_data_to_storage (sp_calc_model_t &calc_model, item_t dt)
      {
        /*
        const sp_fip_results_storage &f_res (calc_model->fip_res);
        int i;
        fip_data *fd;
        int fip_region = rsv_status->fip_region;
        float fip_params[FIP_D_PARAM_COUNT] = {0};

        CH_PTR (rsv_status);

        for (i = 0; i < fip_region; i++)
          {
            fip_results &fr = f_res.fip[i];
            fd = &(rsv_status->current_fip[i]);

            fr.dates.push_back (time + time_step_length);

            fip_params[FIP_D_PARAM_COIP] += (float)fd->COIP;
            fip_params[FIP_D_PARAM_CMOIP] += (float)fd->CMOIP;
            fip_params[FIP_D_PARAM_CWIP] += (float)fd->CWIP;
            fip_params[FIP_D_PARAM_CGIP] += (float)fd->CGIP;
            fip_params[FIP_D_PARAM_HPROD] += (float)fd->HPROD;
            fip_params[FIP_D_PARAM_CPVOL] += (float)fd->CPVOL;
        #ifdef _MPI
            fip_params[FIP_D_PARAM_OIL_IN] += (float)fd->mpi_oil_in;
            fip_params[FIP_D_PARAM_WATER_IN] += (float)fd->mpi_water_in;
            fip_params[FIP_D_PARAM_WATER_OUT] += (float)fd->mpi_water_out;
        #else //_MPI
            fip_params[FIP_D_PARAM_OIL_IN] += (float)fd->oil_in;
            fip_params[FIP_D_PARAM_WATER_IN] += (float)fd->water_in;
            fip_params[FIP_D_PARAM_WATER_OUT] += (float)fd->water_out;
        #endif //_MPI
            fip_params[FIP_D_PARAM_SWAV] += (float)fd->SWAV;
            fip_params[FIP_D_PARAM_SOAV] += (float)fd->SOAV;
            fip_params[FIP_D_PARAM_SGAV] += (float)fd->SGAV;

            fr.d_params[FIP_D_PARAM_COIP].push_back ((float)fd->COIP);
            fr.d_params[FIP_D_PARAM_CMOIP].push_back ((float)fd->CMOIP);
            fr.d_params[FIP_D_PARAM_CWIP].push_back ((float)fd->CWIP);
            fr.d_params[FIP_D_PARAM_CGIP].push_back ((float)fd->CGIP);
            fr.d_params[FIP_D_PARAM_HPROD].push_back ((float)fd->HPROD);
            fr.d_params[FIP_D_PARAM_CPVOL].push_back ((float)fd->CPVOL);
            fr.d_params[FIP_D_PARAM_PAV].push_back ((float)fd->PAV);
        #ifdef _MPI
            fr.d_params[FIP_D_PARAM_OIL_IN].push_back ((float)fd->mpi_oil_in);
            fr.d_params[FIP_D_PARAM_WATER_IN].push_back ((float)fd->mpi_water_in);
            fr.d_params[FIP_D_PARAM_WATER_OUT].push_back ((float)fd->mpi_water_out);
        #else //_MPI
            fr.d_params[FIP_D_PARAM_OIL_IN].push_back ((float)fd->oil_in);
            fr.d_params[FIP_D_PARAM_WATER_IN].push_back ((float)fd->water_in);
            fr.d_params[FIP_D_PARAM_WATER_OUT].push_back ((float)fd->water_out);
        #endif //_MPI
            fr.d_params[FIP_D_PARAM_SWAV].push_back ((float)fd->SWAV);
            fr.d_params[FIP_D_PARAM_SOAV].push_back ((float)fd->SOAV);
            fr.d_params[FIP_D_PARAM_SGAV].push_back ((float)fd->SGAV);
            fr.d_params[FIP_D_PARAM_GOR].push_back ((float)fd->GOR);
          }

        fip_results &fr = f_res.fip[fip_region];

        fr.dates.push_back (time + time_step_length);
        fr.d_params[FIP_D_PARAM_COIP].push_back (fip_params[FIP_D_PARAM_COIP]);
        fr.d_params[FIP_D_PARAM_CMOIP].push_back (fip_params[FIP_D_PARAM_CMOIP]);
        fr.d_params[FIP_D_PARAM_CWIP].push_back (fip_params[FIP_D_PARAM_CWIP]);
        fr.d_params[FIP_D_PARAM_CGIP].push_back (fip_params[FIP_D_PARAM_CGIP]);
        fr.d_params[FIP_D_PARAM_HPROD].push_back (fip_params[FIP_D_PARAM_HPROD]);
        fr.d_params[FIP_D_PARAM_CPVOL].push_back (fip_params[FIP_D_PARAM_CPVOL]);
        fr.d_params[FIP_D_PARAM_PAV].push_back (float (rsv_status->pav));
        fr.d_params[FIP_D_PARAM_OIL_IN].push_back (fip_params[FIP_D_PARAM_OIL_IN]);
        fr.d_params[FIP_D_PARAM_WATER_IN].push_back (fip_params[FIP_D_PARAM_WATER_IN]);
        fr.d_params[FIP_D_PARAM_WATER_OUT].push_back (fip_params[FIP_D_PARAM_WATER_OUT]);
        fr.d_params[FIP_D_PARAM_SWAV].push_back (fip_params[FIP_D_PARAM_SWAV]);
        fr.d_params[FIP_D_PARAM_SOAV].push_back (fip_params[FIP_D_PARAM_SOAV]);
        fr.d_params[FIP_D_PARAM_SGAV].push_back (fip_params[FIP_D_PARAM_SGAV]);
        fr.d_params[FIP_D_PARAM_GOR].push_back (float (rsv_status->gor));
        */
      }
    };

} //blue_sky
#endif //__BS_FIP_RESULTS_STORAGE_H
