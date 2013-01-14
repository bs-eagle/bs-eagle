/**
 *       \file  explicit_model.cpp
 *      \brief  EXPLICIT Initialization model
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  03.05.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "stdafx.h"
#include "explicit_model.hpp"
#include "data_class.h"
#include "rs_mesh_iface.h"
#include "calc_model.h"

namespace blue_sky 
{
  void
  init_pressure (BS_SP (calc_model) model, BS_SP (idata) data, BS_SP (rs_mesh_iface) mesh)
  {
    BS_ASSERT (model);
    BS_ASSERT (data);
    BS_ASSERT (mesh);
    
    // FIXME: don't check size, check init mode!
    if (data->prvd.size ())
      {
        t_long n_cells = mesh->get_n_active_elements();
        if (!data->contains_i_array ("EQLNUM"))
          {
            spv_int eqlnum = data->create_i_array ("EQLNUM", &i_pool_sizes[ARRAY_POOL_TOTAL * EQLNUM], i_pool_default_values[EQLNUM]);
            eqlnum->assign (1);
          }
        spv_int eqlnum_ = data->get_i_array ("EQLNUM");
        t_int const *eqlnum = eqlnum_->data ();
        t_float const *depths = mesh->get_depths ()->data ();

        t_double *pressure = model->pressure->data ();
        if (static_cast <t_long> (model->pressure->size ()) != n_cells)
          {
            bs_throw_exception (boost::format ("Pressure array size mismatch: %d == %d")
                                % model->pressure->size () % n_cells);
          }

        //interpolate by depth
        t_long const *int_to_ext = mesh->get_int_to_ext ()->data ();
        for (t_long i = 0; i < n_cells; ++i)
          {
            t_long i_eql = eqlnum[int_to_ext[i]] - 1;
            pressure[i] = data->prvd[i_eql].interpolate_linear (depths[i]);
          }
      }
    else if (data->contains_fp_array ("PRESSURE"))
      {
        convert_arrays (mesh->get_n_active_elements (), mesh->get_int_to_ext (), *(model->pressure), data->get_fp_array ("PRESSURE"));
      }
    else
      {
        bs_throw_exception ("Initial pressure has not been specified");
      }
  }

  void
  init_saturation (BS_SP (calc_model) model, BS_SP (idata) data, BS_SP (rs_mesh_iface) mesh)
  {
    BS_ASSERT (model);
    BS_ASSERT (data);
    BS_ASSERT (mesh);
    BS_ASSERT (model->n_phases > 1) (model->n_phases);

    const t_long d_w = model->phase_d[FI_PHASE_WATER];
    const t_long d_o = model->phase_d[FI_PHASE_OIL];
    const t_long d_g = model->phase_d[FI_PHASE_GAS];
    int fix_soil_bug = model->ts_params->get_bool(fi_params::FIX_SOIL_BUG);

    model->saturation_3p->init (mesh->get_n_active_elements() * model->n_phases, 0);
    t_double *saturation_3p = model->saturation_3p->data ();
    t_long const n_phases = model->n_phases;

    t_long const *original_element_num = mesh->get_int_to_ext()->data ();
    t_long n = mesh->get_n_active_elements();

    // check number of unknowns
    int un_counter = 0
      + (model->is_water () && !data->contains_fp_array ("SWAT"))
      + (model->is_oil ()   && !data->contains_fp_array ("SOIL"))
      + (model->is_gas ()   && !data->contains_fp_array ("SGAS"))
      ;
    if (un_counter > 1)
      {
        bs_throw_exception ("Not enought phases saturations specified");
      }

    spv_float soil_ = data->get_fp_array ("SOIL");
    spv_float swat_ = data->get_fp_array ("SWAT");
    spv_float sgas_ = data->get_fp_array ("SGAS");

    t_float const *soil = soil_->data ();
    t_float const *swat = swat_->data ();
    t_float const *sgas = sgas_->data ();

    // FIXME: work with ->size ()!

    // 2ph water oil
    if (n_phases == 2 && model->is_water () && model->is_oil ())
      {
        if (soil_->size () && swat_->size ())
          {
            BOSWARN (section::init_data, level::warning)
            << "Oil saturation will ignored in 2 phase water-oil system." << bs_end;
          }
        else if (!swat_->size () && !soil_->size ())
          {
            bs_throw_exception ("Water or oil saturation has not been specified");
          }

        for (t_long i = 0; i < n; ++i)
          {
            t_long blk_i = original_element_num[i];
            if (swat_->size ())
              {
                if (swat[blk_i] > (t_double)1.0 + EPS_DIFF || swat[blk_i] < -EPS_DIFF)
                  {
                    BOSWARN (section::init_data, level::warning) << "Water saturation is out of range" << bs_end;
                    //saturation[i] = (t_double)1.0;
                    saturation_3p[i * n_phases + d_w] = (t_double)1.0;
                    saturation_3p[i * n_phases + d_o] = (t_double)0;
                  }
                else
                  {
                    //saturation[i] = swat[blk_i];
                    saturation_3p[i * n_phases + d_w] = swat[blk_i];
                    saturation_3p[i * n_phases + d_o] = (t_double)1.0 - swat[blk_i];
                  }
              }
            else if (soil_->size ())
              {
                if (soil[blk_i] > (t_double)1.0 + EPS_DIFF || soil[blk_i] < -EPS_DIFF)
                  {
                    BOSWARN (section::init_data, level::warning) << "Oil saturation is out of range" << bs_end;
                    //saturation[i] = (t_double)1.0;
                    saturation_3p[i * n_phases + d_w] = (t_double)0;
                    saturation_3p[i * n_phases + d_o] = (t_double)1.0;
                  }
                else
                  {
                    //saturation[i] = (t_double)1. - swat[blk_i];
                    saturation_3p[i * n_phases + d_w] = (t_double)1.0 - soil[blk_i];
                    saturation_3p[i * n_phases + d_o] = soil[blk_i];
                  }
              }
          }
      }
    // 2ph water gas system
    else if (n_phases == 2 && model->is_water () && model->is_gas ())
      {
        if (swat_->size () && sgas_->size ())
          {
            BOSWARN (section::init_data, level::warning)
            << "Gas saturation will be ignored in 2 phase water-gas system." << bs_end;
          }
        if (!swat_->size () && !sgas_->size ())
          {
            bs_throw_exception ("Water or gas saturation has not been specified");
          }

        for (t_long i = 0; i < n; ++i)
          {
            t_long blk_i = original_element_num[i];
            if (swat_->size ())
              {
                if (swat[blk_i] > (t_double)1.0 + EPS_DIFF || swat[blk_i] < -EPS_DIFF)
                  {
                    BOSWARN (section::init_data, level::warning)
                    << "Water saturation is out of range." << bs_end;
                    //saturation[i] = (t_double)1.0;
                    saturation_3p[i * n_phases + d_w] = (t_double)1.0;
                    saturation_3p[i * n_phases + d_g] = (t_double)0;
                  }
                else
                  {
                    //saturation[i] = swat[blk_i];
                    saturation_3p[i * n_phases + d_w] = swat[blk_i];
                    saturation_3p[i * n_phases + d_g] = (t_double)1.0 - swat[blk_i];
                  }
              }
            else if (sgas_->size ())
              {
                if (sgas[blk_i] > (t_double)1.0 + EPS_DIFF || sgas[blk_i] < -EPS_DIFF)
                  {
                    BOSWARN (section::init_data, level::warning)
                    << "Gas saturation is out of range." << bs_end;
                    //saturation[i] = (t_double)1.0;
                    saturation_3p[i * n_phases + d_w] = (t_double)0;
                    saturation_3p[i * n_phases + d_g] = (t_double)1.0;
                  }
                else
                  {
                    //saturation[i] = (t_double)1. - sgas[blk_i];
                    saturation_3p[i * n_phases + d_g] = sgas[blk_i];
                    saturation_3p[i * n_phases + d_w] = (t_double)1.0 - sgas[blk_i];
                  }
              }
          }
      }
    // 2ph oil gas system
    else if (n_phases == 2 && model->is_oil () && model->is_gas ())
      {
        if (soil_->size () && sgas_->size ())
          {
            BOSWARN (section::init_data, level::warning)
            << "Oil saturation will be ignored in 2 phase gas-oil system." << bs_end;
          }
        else if (!soil_->size () && !sgas_->size ())
          {
            bs_throw_exception ("Gas or oil saturation has not been specified");
          }

        for (t_long i = 0; i < n; ++i)
          {
            t_long blk_i = original_element_num[i];
            if (soil_->size ())
              {
                if (soil[blk_i] > (t_double) + EPS_DIFF || soil[blk_i] < -EPS_DIFF)
                  {
                    BOSWARN (section::init_data, level::warning)
                    << "Oil saturation is out of range." << bs_end;
                    //saturation[i] = (t_double)0.0;
                    saturation_3p[i * n_phases + d_g] = (t_double)0;
                    saturation_3p[i * n_phases + d_o] = (t_double)1.0;
                  }
                else
                  {
                    //saturation[i] = (t_double)1.0 - soil[blk_i];
                    saturation_3p[i * n_phases + d_g] = (t_double)1.0 - soil[blk_i];
                    saturation_3p[i * n_phases + d_o] = soil[blk_i];
                  }
              }
            else if (sgas_->size ())
              {
                if (sgas[blk_i] > (t_double)1.0 + EPS_DIFF || sgas[blk_i] < -EPS_DIFF)
                  {
                    BOSWARN (section::init_data, level::warning)
                    << "Gas saturation is out of range." << bs_end;
                    //saturation[i] = (t_double)1.0;
                    saturation_3p[i * n_phases + d_g] = (t_double)1.0;
                    saturation_3p[i * n_phases + d_o] = (t_double)0;
                  }
                else
                  {
                    //saturation[i] = sgas[blk_i];
                    saturation_3p[i * n_phases + d_o] = (t_double)1.0 - sgas[blk_i];
                    saturation_3p[i * n_phases + d_g] = sgas[blk_i];
                  }
              }
          }
      }
    // 3ph water oil gas system
    else if (n_phases == 3)
      {
        if (swat_->size () && soil_->size () && sgas_->size ())
          {
            BOSWARN (section::init_data, level::warning)
            << "Oil saturation will be ignored in 3 phase water-gas-oil system." << bs_end;
          }

        if (swat_->size () && sgas_->size ())
          {
            for (t_long i = 0; i < n; ++i)
              {
                t_long blk_i = original_element_num[i];
                t_double sw = swat[blk_i];
                t_double sg = sgas[blk_i];
                if ((sw + sg) < -EPS_DIFF || (sw + sg) > 1 + EPS_DIFF)
                  {
                    bs_throw_exception ("Gas saturation plus water saturation is out of range");
                  }
                if (fix_soil_bug && fabs (sw - (t_double)1.0) < (t_double)0.01)
                  {
                    sw = static_cast <t_double> (0.99);
                  }

                saturation_3p[i * n_phases + d_w] = sw;
                saturation_3p[i * n_phases + d_g] = sg;
                saturation_3p[i * n_phases + d_o] = (t_double)1.0 - sw - sg;
              }
          }

        else if (swat_->size () && soil_->size ())
          {
            for (t_long i = 0; i < n; ++i)
              {
                t_long blk_i = original_element_num[i];
                t_double sw = swat[blk_i];
                t_double sg = (t_double)1.0 - sw - soil[blk_i];
                if (sg < -EPS_DIFF || sg > 1 + EPS_DIFF)
                  {
                    bs_throw_exception ("Oil saturation plus water saturation is out of range");
                  }
                if (fix_soil_bug && fabs(sw - (t_double)1.0) < (t_double)0.001)
                  {
                    sw = 0.999f;
                  }

                saturation_3p[i * n_phases + d_w] = sw;
                saturation_3p[i * n_phases + d_g] = sg;
                saturation_3p[i * n_phases + d_o] = (t_double)1.0 - sw - sg;
              }
          }

        else if (sgas_->size () && soil_->size ())
          {
            for (t_long i = 0; i < n; ++i)
              {
                t_long blk_i = original_element_num[i];
                t_double sg = sgas[blk_i];
                t_double sw = (t_double)1.0 - sg - soil[blk_i];
                if (sw < -EPS_DIFF || sw > 1 + EPS_DIFF)
                  {
                    bs_throw_exception ("Oil saturation plus gas saturation is out of range");
                  }
                if (fix_soil_bug && fabs(sw - (t_double)1.0) < (t_double)0.01)
                  {
                    sw = static_cast <t_double> (0.99);
                  }

                saturation_3p[i * n_phases + d_w] = sw;
                saturation_3p[i * n_phases + d_g] = sg;
                saturation_3p[i * n_phases + d_o] = (t_double)1.0 - sw - sg;
              }
          }
        else if (swat_->size ())
          {
            for (t_long i = 0; i < n; ++i)
              {
                t_long blk_i = original_element_num[i];
                t_double sw = swat[blk_i];
                t_double sg = (t_double)0.0;
                if (sw > 1 + EPS_DIFF || sw < -EPS_DIFF)
                  {
                    bs_throw_exception ("Water saturation is out of range");
                  }
                if (fix_soil_bug && fabs(sw - (t_double)1.0) < (t_double)0.01)
                  {
                    sw = static_cast <t_double> (0.99);
                  }

                saturation_3p[i * n_phases + d_w] = sw;
                saturation_3p[i * n_phases + d_g] = sg;
                saturation_3p[i * n_phases + d_o] = (t_double)1.0 - sw - sg;
              }
          }

        else if (soil_->size ())
          {
            for (t_long i = 0; i < n; ++i)
              {
                t_long blk_i = original_element_num[i];
                t_double sw = (t_double)1.0 - soil[blk_i];
                t_double sg = (t_double)0.0;
                if (soil[blk_i] > 1 + EPS_DIFF || soil[blk_i] < -EPS_DIFF)
                  {
                    bs_throw_exception ("Oil saturation is out of range");
                  }
                if (fix_soil_bug && fabs(sw - (t_double)1.0) < (t_double)0.01)
                  {
                    sw = static_cast <t_double> (0.99);
                  }

                saturation_3p[i * n_phases + d_w] = sw;
                saturation_3p[i * n_phases + d_g] = sg;
                saturation_3p[i * n_phases + d_o] = (t_double)1.0 - sw - sg;
              }
          }

        else if (sgas_->size ())
          {
            for (t_long i = 0; i < n; ++i)
              {
                t_long blk_i = original_element_num[i];
                t_double sw = (t_double)1.0 - sgas[blk_i];
                t_double sg = sgas[blk_i];
                if (sg > 1 + EPS_DIFF || sg < -EPS_DIFF)
                  {
                    bs_throw_exception ("Gas saturation is out of range");
                  }
                if (fix_soil_bug && fabs(sw - (t_double)1.0) < (t_double)0.01)
                  {
                    sw = static_cast <t_double> (0.99);
                  }

                saturation_3p[i * n_phases + d_w] = sw;
                saturation_3p[i * n_phases + d_g] = sg;
                saturation_3p[i * n_phases + d_o] = (t_double)1.0 - sw - sg;
              }
          }
        else
          {
            bs_throw_exception ("None of saturation specified");
          }
      }
  }

  void
  init_rs (BS_SP (calc_model) model, BS_SP (idata) data, BS_SP (rs_mesh_iface) mesh)
  {
    BS_ASSERT (model);
    BS_ASSERT (data);
    BS_ASSERT (mesh);
    BS_ASSERT (model->is_oil () && model->is_gas ()) (model->is_oil ()) (model->is_gas ());

    spv_float init_pbub_ = data->get_fp_array ("PBUB");
    spv_float init_rs_   = data->get_fp_array ("RS");
    if (!init_pbub_->size () && !init_rs_->size ())
      {
        bs_throw_exception ("Should be specified init_pbub or init_rs");
      }

    t_float const *init_pbub = init_pbub_->data ();
    t_float const *init_rs = init_rs_->data ();

    t_long const *int_to_ext = mesh->get_int_to_ext ()->data ();
    t_double *saturation_3p = model->saturation_3p->data ();
    t_double *gas_oil_ratio = model->gas_oil_ratio->data ();
    t_double const *pressure = model->pressure->data ();
    t_long const n_phases = model->n_phases;

    // main loop through all cells
    t_double cell_pbub = 0;
    for (t_long i = 0, cell_count = mesh->get_n_active_elements(); i < cell_count; ++i)
      {
        // calculate index of gas saturation and oil PVT prop
        t_long i_g = FI_PH_IND (i, model->phase_d[FI_PHASE_GAS], n_phases);
        t_long i_o = FI_PH_IND (i, model->phase_d[FI_PHASE_OIL], n_phases);
        // get cell PVT region index
        t_long reg = (*model->pvt_regions)[i];

        // get PVTO for cell
        BS_SP (pvt_base) pvto = model->pvt_oil_array[reg];
        if (!pvto)
          continue;

        // get cell index in original arrays
        t_long cell_i = int_to_ext[i];

        if (saturation_3p[i_g] < EPS_DIFF && saturation_3p[i_o] < EPS_DIFF)
          {
            gas_oil_ratio[i] = 0;
            model->main_variable[i] = FI_MOMG_VAR;
          }
        // if S_g > 0 --- main var is S_g and system include 3 phases
        // P_bub = P and Rs = Rs (P_bub)
        else if (saturation_3p[i_g] > EPS_DIFF)
          {
            if (init_pbub_->size () && init_pbub[cell_i] < pressure[i])
              cell_pbub = init_pbub[cell_i];
            else
              cell_pbub = pressure[i];

            gas_oil_ratio[i] = pvto->interpolate_and_fix (cell_pbub);

            // set main variable to GAS saturation
            model->main_variable[i] = FI_SG_VAR;
          }
        else
          {
            // if S_g == 0 --- main variable is Rs,
            // if user specify initial bubble point (PBUB)
            if (init_pbub_->size ())
              {
                // initial bubble point can not be greater than cell pressure
                if (init_pbub[cell_i] < pressure[i])
                  cell_pbub = init_pbub[cell_i];
                else
                  cell_pbub = pressure[i];

                gas_oil_ratio[i] = pvto->interpolate_and_fix (cell_pbub);
              }
            // if user specify initial RS
            else if (init_rs_->size ())
              {
                gas_oil_ratio[i] = init_rs[cell_i];
              }
            else
              {
                BS_ASSERT (false && "init_pbub and init_rs is not specified");
              }

            // set main variable to RS
            model->main_variable[i] = FI_RO_VAR;
          }
      }
  }

  explicit_model::explicit_model (bs_type_ctor_param)
  {
  }

  explicit_model::explicit_model (explicit_model const &src)
  : bs_refcounter (src), init_model_iface (src)
  {
  }

  void
  explicit_model::init (BS_SP (calc_model) model, BS_SP (idata) data, BS_SP (rs_mesh_iface) mesh)
  {
    init_pressure (model, data, mesh);
    if (model->n_phases > 1)
      {
        init_saturation (model, data, mesh);
      }
    if (model->is_oil () && model->is_gas ())
      {
        init_rs (model, data, mesh);
      }
  }

  BLUE_SKY_TYPE_STD_CREATE (explicit_model);
  BLUE_SKY_TYPE_STD_COPY (explicit_model);
  BLUE_SKY_TYPE_IMPL (explicit_model, init_model_iface, "explicit_init_model", "explicit_init_model", "EXPLICIT initialization model");
}

