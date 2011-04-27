/**
* @file calc_equil.cpp
* @brief calculating capillary equilibration
* @author Kamaltinova A.
* @date 2008-11-26
*/
#include "stdafx.h"

#include "calc_model.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "data_class.h"
#include "scal_3p.h"
#include "scale_array_holder.h"
#include "rock_grid.h"
#include "scal_region_info.h"
#include "scal_region.h"
#include "scal_2p_data_holder.h"
#include "pvt_dead_oil.h"
#include "pvt_water.h"
#include "pvt_gas.h"
#include "plane_orientation.h"
#include "jfunction.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {

  //! calc pressure on current depth
  int
  calc_model::equil_calc_pressure (item_t prev_press, item_t cur_d, item_t h, index_t phase, index_t i_pvt,
      double rs_type, item_t depth_goc, item_t rs_dat,
      val_vs_depth *rsvd, val_vs_depth *pbvd,
      /*output*/ item_t &p, item_t *rs)
  {
    item_t p0, p1, rho, pb;
    item_t invers_fvf = 0;
    item_t rs_max;
    index_t is_g = FI_CHK_GAS (phases);

    if (phase == FI_PHASE_OIL && is_g)
      {
        if (!rs)
          return -1;

        *rs = rs_dat;
      }

    p0 = 0.;
    p1 = prev_press;
    while (fabs (p1 - p0) > EPS_DIFF)
      {
        p0 = p1;

        if (phase == FI_PHASE_WATER)  //water phase
          {
            pvt_water_array[i_pvt]->calc (0.5 * (prev_press + p0), &invers_fvf, 0, 0, 0, 0, 0);
            rho = pvt_water_array[i_pvt]->get_surface_density () * invers_fvf;
          }
        else if (phase == FI_PHASE_GAS)  //gas phase
          {
            pvt_gas_array[i_pvt]->calc (0.5 * (prev_press + p0), &invers_fvf, 0, 0, 0, 0, 0);
            rho = pvt_gas_array[i_pvt]->get_surface_density () * invers_fvf;
          }
        else  //oil phase
          {
            pvt_oil_array[i_pvt]->calc (is_g, !is_g || (cur_d - 0.5 * h < depth_goc) ? FI_SG_VAR : FI_RO_VAR,
                                        0.5 * (prev_press + p0), rs ? *rs : 0, &invers_fvf,
                                        0, 0, 0, 0, 0, 0, 0, 0,
                                        &rs_max, 0);

            rho = pvt_oil_array[i_pvt]->get_surface_density () * invers_fvf;

            //calc Rs if gas phase is present
            if (is_g)
              {
                //------------------------------- calc Rs ------------------------------------------
                // this is undersaturated oil case under the gas-oil contact
                if (cur_d > depth_goc)
                  {
                    if (rs_type > 0)
                      {
                        if (rsvd)
                          *rs = rsvd->interpolate_linear (cur_d);
                        else if (pbvd)
                          {
                            // bubble-point pressure
                            pb = pbvd->interpolate_linear (cur_d);
                            // calculate Rs which corresponds to Pbub for undersaturated case
                            *rs = pvt_oil_array[i_pvt]->interpolate_and_fix (pb);
                          }
                        else
                          {
                            BOSERR (section::init_data, level::warning)
                            << "Error: you should specify either RSVD or PBVD keyword to use option 7 in keyword EQUIL." << bs_end;
                            return -1;
                          }
                        // Rs calculated by RSVD or PBVD keywords should be less than saturated RS value at the local pressure
                        if (*rs > rs_max)
                          *rs = rs_max;
                      }
                    // the dissolved gas concentration in undersaturated oil equal to the saturated Rs value at the gas-oil contact.
                    else
                      *rs = rs_dat;
                  }
                // this is saturated oil case
                else
                  // use Rs value for saturated oil at the local pressure above the gas-oil contact
                  *rs = rs_max;
                //----------------------------------------------------------------------------------

                rho += *rs * pvt_gas_array[i_pvt]->get_surface_density () * invers_fvf;
              }
          }

        p1 = prev_press + rho * internal_constants.gravity_constant * h;
      }

    p = (p0 + p1) * 0.5;

    return 0;
  }

  // initialize initial conditions
  int
  calc_model::calc_equil (const sp_idata_t &input_data, const sp_mesh_iface_t &mesh)
  {
    typedef std::vector <float> vector_t;

    const sp_idata_t &l_idata(input_data);
    const sp_mesh_iface_t &l_mesh(mesh);

    const index_t n_depth = 100, n_layer = 10;
    index_t i_depth;
    index_t i_eql, i_pvt, i_sat = 0, n_eql;
    index_t i_cell, n_cells;
    stdv_long eqlnum;
    stdv_double min_depth, max_depth;
    item_t d_depth, h;
    item_t depth_dat, press_dat, rs_dat, depth_woc, press_woc, depth_goc, press_goc;
    index_t il, iu, im;
    item_t prev_d, prev_p, cur_d;
    index_t is_w, is_g, is_o, d_w, d_g, d_o, ds_w, ds_g;
    index_t main_phase;
    stdv_double depth, press, rs;      // array of depthes, phases pressure, gas oil ratios
    item_t p_oil, p_water, p_gas;
    item_t p_p[3], s_p[2], pc_limit[2], dcoef;

    const stdv_double &perm = rock_grid_prop->permeability;
    const stdv_double &poro = rock_grid_prop->porosity_p_ref;
    idata::vval_vs_depth &rsvd = l_idata->get_rsvd ();
    idata::vval_vs_depth &pbvd = l_idata->get_pbvd ();
    const spv_int &equil_regions = l_idata->equil_regions;
    const index_t N_equil = 2;

    BS_ASSERT (equil_regions->size ());
    BS_ASSERT (l_idata->contains_i_array ("EQLNUM"));

    is_w = FI_CHK_WATER (phases);
    is_g = FI_CHK_GAS (phases);
    is_o = FI_CHK_OIL (phases);
    d_w = phase_d[FI_PHASE_WATER];
    d_o = phase_d[FI_PHASE_OIL];
    d_g = phase_d[FI_PHASE_GAS];
    ds_w = sat_d[FI_PHASE_WATER];
    ds_g = sat_d[FI_PHASE_GAS];

    index_t i_original_cell, i_layer;
    item_t depth_top, depth_bottom, depth_center;
    item_t depth_step, layer_center;
    item_t s_water, s_gas;
    stdv_double press_w, press_g, press_o;
    spv_float swatinit;
    stdv_double pcw;
    bool is_swatinit = false;

    //get num of eql regions
    n_eql = equil_regions->size (); // FIXME: sergey.miryanov: get number of EQUIL regs
    //get num of elements
    n_cells = l_mesh->get_n_active_elements();

    //get eqlnum array
    eqlnum.resize(n_cells);
    //todo: if zero array
    convert_arrays (mesh->get_n_active_elements (), mesh->get_int_to_ext (), eqlnum, l_idata->get_i_array ("EQLNUM"));

    min_depth.resize(n_eql, 0);
    max_depth.resize(n_eql, 0);
    depth.resize(n_eql * n_depth, 0);
    press.resize(n_eql * n_phases * n_depth, 0);

    if (is_o && is_g)
      rs.resize(n_eql * n_depth, 0);

    //set min and max to datum depth
    for (i_eql = 0; i_eql < n_eql; i_eql++)
      min_depth[i_eql] = max_depth[i_eql] = (*l_idata->equil)[EQUIL_TOTAL * i_eql + EQUIL_DAT_DEPTH];

    const spv_float &cell_depths = l_mesh->get_depths ();

    //get min and max depth in each of eql regions
    for (i_cell = 0; i_cell < n_cells; i_cell++)
      {
        item_t top_depth = 0., bottom_depth = 0.;
        i_eql = eqlnum[i_cell] - 1;
                
        top_depth = (*cell_depths)[i_cell] - 0.5 * l_mesh->get_element_dim3_size (i_cell);
        bottom_depth = (*cell_depths)[i_cell] + 0.5 * l_mesh->get_element_dim3_size (i_cell);

        if (top_depth < min_depth[i_eql])
          min_depth[i_eql] = top_depth;
        if (bottom_depth > max_depth[i_eql])
          max_depth[i_eql] = bottom_depth;
      }

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////// CALCULATE Pp AND Rs BY DEPTH /////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    //loop through all eql regions
    for (i_eql = 0; i_eql < n_eql; i_eql++)
      {
        i_pvt = (*equil_regions)[N_equil * i_eql + 1];    // pvt region for current equil region

        depth_dat = (*l_idata->equil)[EQUIL_TOTAL * i_eql + EQUIL_DAT_DEPTH];
        depth_woc = (*l_idata->equil)[EQUIL_TOTAL * i_eql + EQUIL_WOC_DEPTH];
        depth_goc = (*l_idata->equil)[EQUIL_TOTAL * i_eql + EQUIL_GOC_DEPTH];
        press_dat = (*l_idata->equil)[EQUIL_TOTAL * i_eql + EQUIL_DAT_PRESS];

        //set main_phase
        main_phase = FI_PHASE_OIL;
        if (n_phases > 1 || is_o)
          main_phase = FI_PHASE_OIL;
        else if (is_w)
          main_phase = FI_PHASE_WATER;
        else if (is_g)
          main_phase = FI_PHASE_GAS;
        if (n_phases > 1)
          {
            if (is_w && (depth_dat > depth_woc))
              main_phase = FI_PHASE_WATER;
            else if (is_g && (depth_dat < depth_goc))
              main_phase = FI_PHASE_GAS;
          }

        //set depth step
        d_depth = (max_depth[i_eql] - min_depth[i_eql]) / (item_t)(n_depth - 1);

        //set dep array values
        for (i_depth = 0; i_depth < n_depth; i_depth++)
          depth[i_depth + i_eql * n_depth] = min_depth[i_eql] + i_depth * d_depth;

        // calc Rs at datum pressure from PVTO table
        rs_dat = 0.;
        if (is_o && is_g)
          rs_dat = pvt_oil_array[i_pvt]->interpolate_and_fix (press_dat);

        // search Ddat in dep[] array
        // Ddat can't exceed [Dmin; Dmax]
        BINARY_SEARCH (&depth[i_eql * n_depth], depth_dat, n_depth, iu, im, il);
        if (il == 0)
          iu = 1;
        else
          {
            iu = il;
            il--;
          }

        //---- calc pressure table of main phase ---------
        for (index_t up_down = -1; up_down < 2; up_down += 2)
          {
            //up_down = -1 - go up from datum depth
            //up_down = 1 - go down from datum depth
            if (up_down < 0)
              i_depth = il;
            else
              i_depth = iu;

            prev_d = depth_dat;
            prev_p = press_dat;

            while (i_depth >= 0 && i_depth < n_depth)
              {
                cur_d = depth[i_depth + i_eql * n_depth];
                h = cur_d - prev_d;

                if (equil_calc_pressure (prev_p, cur_d, h, main_phase, i_pvt,
                                         is_g ? (*l_idata->equil)[EQUIL_TOTAL * i_eql + EQUIL_RS_TYPE] : 0,
                                         is_g ? depth_goc : 0, is_g ? rs_dat : 0,
                                         (is_g && rsvd.size ()) ? &rsvd[i_eql] : 0, (is_g && pbvd.size ()) ? &pbvd[i_eql] : 0,
                                         press[i_depth + (i_eql * n_phases + phase_d[main_phase]) * n_depth],
                                         rs.size () ? &rs[i_depth + i_eql * n_depth] : 0))
                  return -1;

                prev_p = press[i_depth + (i_eql * n_phases + phase_d[main_phase]) * n_depth];
                prev_d = cur_d;

                i_depth += up_down;
              }
          }

        if (n_phases > 1)
          {
            item_t p_temp, d_temp, rs_temp;

            if (main_phase != FI_PHASE_OIL)
              {
                //calc oil pressure
                if (main_phase == FI_PHASE_WATER)
                  d_temp = depth_woc;
                else
                  d_temp = depth_goc;

                //search d_temp in depth array
                //Dwoc/Dgoc can exceed [Dmin; Dmax]
                if (d_temp < min_depth[i_eql] || d_temp > max_depth[i_eql])
                  {
                    //calc water/gas pressure at WOC/GOC
                    if (d_temp < min_depth[i_eql])
                      {
                        //get first values from tables
                        prev_d = min_depth[i_eql];
                        prev_p = press[0 + (i_eql * n_phases + phase_d[main_phase]) * n_depth];

                        il = -1;
                        iu = 0;
                      }
                    else
                      {
                        //get last values from tables
                        prev_d = max_depth[i_eql];
                        prev_p = press[n_depth - 1 + (i_eql * n_phases + phase_d[main_phase]) * n_depth];

                        il = n_depth - 1;
                        iu = n_depth;
                      }

                    h = d_temp - prev_d;

                    if (equil_calc_pressure (prev_p, d_temp, h, main_phase, i_pvt,
                                             0, 0, 0, 0, 0,
                                             p_temp))
                      return -1;

                  }
                else
                  {
                    // search Dwoc/Dgoc in dep[] array and interpolate pressure of main phase
                    BINARY_SEARCH (&depth[i_eql * n_depth], d_temp, n_depth, iu, im, il);
                    INTERNAL_INTERPOLATION (&depth[i_eql * n_depth], &press[(i_eql * n_phases + phase_d[main_phase]) * n_depth],
                                            d_temp, p_temp, n_depth);

                    if (il == 0)
                      iu = 1;
                    else
                      {
                        iu = il;
                        il--;
                      }
                  }

                //oil pressure at WOC/GOC
                if (main_phase == FI_PHASE_WATER)
                  p_temp += (*l_idata->equil)[EQUIL_TOTAL * i_eql + EQUIL_WOC_PRESS];
                else
                  p_temp -= (*l_idata->equil)[EQUIL_TOTAL * i_eql + EQUIL_GOC_PRESS];

                //---- calc pressure table of oil -----
                for (index_t up_down = -1; up_down < 2; up_down += 2)
                  {
                    //up_down = -1 - go up from datum depth
                    //up_down = 1 - go down from datum depth
                    if (up_down < 0)
                      i_depth = il;
                    else
                      i_depth = iu;

                    prev_d = d_temp;
                    prev_p = p_temp;

                    while (i_depth >= 0 && i_depth < n_depth)
                      {
                        cur_d = depth[i_depth + i_eql * n_depth];
                        h = cur_d - prev_d;

                        if (equil_calc_pressure (prev_p, cur_d, h, FI_PHASE_OIL, i_pvt,
                                                 is_g ? (*l_idata->equil)[EQUIL_TOTAL * i_eql + EQUIL_RS_TYPE] : 0,
                                                 is_g ? depth_goc : 0, is_g ? rs_dat : 0,
                                                 (is_g && rsvd.size ()) ? &rsvd[i_eql] : 0, (is_g && pbvd.size ()) ? &pbvd[i_eql] : 0,
                                                 press[i_depth + (i_eql * n_phases + phase_d[FI_PHASE_OIL]) * n_depth],
                                                 rs.size () ? &rs[i_depth + i_eql * n_depth] : 0))
                          return -1;

                        prev_p = press[i_depth + (i_eql * n_phases + phase_d[FI_PHASE_OIL]) * n_depth];
                        prev_d = cur_d;

                        i_depth += up_down;
                      }
                  }
              }

            if (is_w && main_phase != FI_PHASE_WATER)
              {
                //calc water pressure

                //search depth_woc in depth array
                //Dwoc can exceed [Dmin; Dmax]
                if (depth_woc < min_depth[i_eql] || depth_woc > max_depth[i_eql])
                  {
                    //calc Po at WOC
                    if (depth_woc < min_depth[i_eql])
                      {
                        //get first values from tables
                        prev_d = min_depth[i_eql];
                        prev_p = press[0 + (i_eql * n_phases + phase_d[FI_PHASE_OIL]) * n_depth];

                        il = -1;
                        iu = 0;
                      }
                    else
                      {
                        //get last values from tables
                        prev_d = max_depth[i_eql];
                        prev_p = press[n_depth - 1 + (i_eql * n_phases + phase_d[FI_PHASE_OIL]) * n_depth];

                        il = n_depth - 1;
                        iu = n_depth;
                      }

                    h = depth_woc - prev_d;

                    //calc oil pressure at WOC
                    if (equil_calc_pressure (prev_p, depth_woc, h, FI_PHASE_OIL, i_pvt,
                                             is_g ? (*l_idata->equil)[EQUIL_TOTAL * i_eql + EQUIL_RS_TYPE] : 0,
                                             is_g ? depth_goc : 0, is_g ? rs_dat : 0,
                                             (is_g && rsvd.size ()) ? &rsvd[i_eql] : 0, (is_g && pbvd.size ()) ? &pbvd[i_eql] : 0,
                                             press_woc, &rs_temp))
                      return -1;

                  }
                else
                  {
                    // search Dwoc in dep[] array and interpolate oil pressure
                    BINARY_SEARCH (&depth[i_eql * n_depth], depth_woc, n_depth, iu, im, il);
                    INTERNAL_INTERPOLATION (&depth[i_eql * n_depth], &press[(i_eql * n_phases + phase_d[FI_PHASE_OIL]) * n_depth],
                                            depth_woc, press_woc, n_depth);

                    if (il == 0)
                      iu = 1;
                    else
                      {
                        iu = il;
                        il--;
                      }
                  }

                //water pressure at WOC
                press_woc -= (*l_idata->equil)[EQUIL_TOTAL * i_eql + EQUIL_WOC_PRESS];

                //---- calc pressure table of water -----
                for (index_t up_down = -1; up_down < 2; up_down += 2)
                  {
                    //up_down = -1 - go up from datum depth
                    //up_down = 1 - go down from datum depth
                    if (up_down < 0)
                      i_depth = il;
                    else
                      i_depth = iu;

                    prev_d = depth_woc;
                    prev_p = press_woc;

                    while (i_depth >= 0 && i_depth < n_depth)
                      {
                        cur_d = depth[i_depth + i_eql * n_depth];
                        h = cur_d - prev_d;

                        if (equil_calc_pressure (prev_p, cur_d, h, FI_PHASE_WATER, i_pvt,
                                                 0, 0, 0, 0, 0,
                                                 press[i_depth + (i_eql * n_phases + phase_d[FI_PHASE_WATER]) * n_depth]))
                          return -1;

                        prev_p = press[i_depth + (i_eql * n_phases + phase_d[FI_PHASE_WATER]) * n_depth];
                        prev_d = cur_d;

                        i_depth += up_down;
                      }
                  }
              }

            if (is_g && main_phase != FI_PHASE_GAS)
              {
                //calc gas pressure

                //search depth_goc in depth array
                //Dgoc can exceed [Dmin; Dmax]
                if (depth_goc < min_depth[i_eql] || depth_goc > max_depth[i_eql])
                  {
                    //calc Po at WOC
                    if (depth_goc < min_depth[i_eql])
                      {
                        //get first values from tables
                        prev_d = min_depth[i_eql];
                        prev_p = press[0 + (i_eql * n_phases + phase_d[FI_PHASE_OIL]) * n_depth];

                        il = -1;
                        iu = 0;
                      }
                    else
                      {
                        //get last values from tables
                        prev_d = max_depth[i_eql];
                        prev_p = press[n_depth - 1 + (i_eql * n_phases + phase_d[FI_PHASE_OIL]) * n_depth];

                        il = n_depth - 1;
                        iu = n_depth;
                      }

                    h = depth_goc - prev_d;

                    //calc oil pressure at GOC
                    if (equil_calc_pressure (prev_p, depth_goc, h, FI_PHASE_OIL, i_pvt,
                                             is_g ? (*l_idata->equil)[EQUIL_TOTAL * i_eql + EQUIL_RS_TYPE] : 0,
                                             is_g ? depth_goc : 0, is_g ? rs_dat : 0,
                                             (is_g && rsvd.size ()) ? &rsvd[i_eql] : 0, (is_g && pbvd.size ()) ? &pbvd[i_eql] : 0,
                                             press_goc, &rs_temp))
                      return -1;

                  }
                else
                  {
                    // search Dgoc in dep[] array and interpolate pressure of oil phase
                    BINARY_SEARCH (&depth[i_eql * n_depth], depth_goc, n_depth, iu, im, il);
                    INTERNAL_INTERPOLATION (&depth[i_eql * n_depth], &press[(i_eql * n_phases + phase_d[FI_PHASE_OIL]) * n_depth],
                                            depth_goc, press_goc, n_depth);

                    if (il == 0)
                      iu = 1;
                    else
                      {
                        iu = il;
                        il--;
                      }
                  }

                //gas pressure at GOC
                press_goc += (*l_idata->equil)[EQUIL_TOTAL * i_eql + EQUIL_GOC_PRESS];

                //---- calc pressure table of gas -----
                for (index_t up_down = -1; up_down < 2; up_down += 2)
                  {
                    //up_down = -1 - go up from datum depth
                    //up_down = 1 - go down from datum depth
                    if (up_down < 0)
                      i_depth = il;
                    else
                      i_depth = iu;

                    prev_d = depth_goc;
                    prev_p = press_goc;

                    while (i_depth >= 0 && i_depth < n_depth)
                      {
                        cur_d = depth[i_depth + i_eql * n_depth];
                        h = cur_d - prev_d;

                        if (equil_calc_pressure (prev_p, cur_d, h, FI_PHASE_GAS, i_pvt,
                                                 0, 0, 0, 0, 0,
                                                 press[i_depth + (i_eql * n_phases + phase_d[FI_PHASE_GAS]) * n_depth]))
                          return -1;

                        prev_p = press[i_depth + (i_eql * n_phases + phase_d[FI_PHASE_GAS]) * n_depth];
                        prev_d = cur_d;

                        i_depth += up_down;
                      }
                  }
              }
          }
      }  //end of loop for equilibration regions

#if 0
    tools::save_seq_vector ("equil_depth_bs.txt").save (depth);
    tools::save_seq_vector ("equil_pressure_bs.txt").save (press);
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////// CALCULATE Pp, Sp AND Rs IN CELLS /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
    //check SWATINIT, JFUNC and PCW
    if (n_phases > 1 && is_w)
      {
        if (l_idata->contains_fp_array ("SWATINIT"))
          {
            is_swatinit = true;
            swatinit = l_idata->get_fp_array ("SWATINIT");

            scal_prop->get_water_jfunction ()->set_valid (false);
            scal_prop->get_water_scale ()->remove_pcp ();
            pcw.resize (n_cells, 0);
          }
        else if (scal_prop->get_water_jfunction ()->valid ())
          scal_prop->get_water_scale ()->remove_pcp ();
      }

    //------------------------------- loop through active cells ----------------------
    for (i_cell = 0; i_cell < n_cells; i_cell++)
      {
        p_oil = 0.;
        p_water = 0.;
        p_gas = 0.;
        if (n_phases > 1 && is_w)
          pc_limit[ds_w] = 0.;
        if (n_phases > 1 && is_g)
          pc_limit[ds_g] = 0.;

        i_original_cell = l_mesh->get_element_int_to_ext (i_cell);
        i_eql = eqlnum[i_cell] - 1;
        if (n_phases > 1)
          i_sat = (*sat_regions)[i_cell];

        item_t gor = 0;
        item_t *rs_reg = 0;
        if (is_o && is_g)
          {
            (*gas_oil_ratio)[i_cell] = 0;
            rs_reg = &rs[i_eql * n_depth];
          }

        if (is_o)
          {
            press_o.assign (&press[(i_eql * n_phases + phase_d[FI_PHASE_OIL]) * n_depth],
                            &press[(i_eql * n_phases + phase_d[FI_PHASE_OIL]) * n_depth] + n_depth);
          }
        if (is_w)
          {
            press_w.assign (&press[(i_eql * n_phases + phase_d[FI_PHASE_WATER]) * n_depth],
                            &press[(i_eql * n_phases + phase_d[FI_PHASE_WATER]) * n_depth] + n_depth);
          }
        if (is_g)
          {
            press_g.assign (&press[(i_eql * n_phases + phase_d[FI_PHASE_GAS]) * n_depth],
                            &press[(i_eql * n_phases + phase_d[FI_PHASE_GAS]) * n_depth] + n_depth);
          }

        depth_center = (*cell_depths)[i_cell];
        depth_top = depth_center - 0.5 * l_mesh->get_element_dim3_size (i_cell);
        depth_bottom = depth_center + 0.5 * l_mesh->get_element_dim3_size (i_cell);

        //calc average values
        depth_step = (depth_bottom - depth_top) / (item_t)n_layer;
        layer_center = depth_top - 0.5 * depth_step;

        s_water = 0.;
        s_gas = 0.;

        for (i_layer = 0; i_layer < n_layer; i_layer++)
          {
            layer_center += depth_step;

            //interpolate oil and phases pressure in layer
            BINARY_SEARCH (&depth[i_eql * n_depth], layer_center, n_depth, iu, im, il);
            if (il == 0)
              {
                if (is_o)
                  {
                    p_p[d_o] = press_o[0];
                    if (is_g)
                      {
                        gor = rs_reg[0];
                      }
                  }
                if (is_w)
                  p_p[d_w] = press_w[0];
                if (is_g)
                  p_p[d_g] = press_g[0];
              }
            else if (il == n_depth)
              {
                if (is_o)
                  {
                   
                    p_p[d_o] = press_o[n_depth - 1];
                    if (is_g)
                      {
                        gor = rs_reg[n_depth - 1];
                      }
                  }
                if (is_w)
                  p_p[d_w] = press_w[n_depth - 1];
                if (is_g)
                  p_p[d_g] = press_g[n_depth - 1];
              }
            else
              {
                iu = il;
                il--;

                dcoef = (depth[il + i_eql * n_depth] - layer_center) / (depth[il + i_eql * n_depth] - depth[iu + i_eql * n_depth]);
                if (is_o)
                  {
                    p_p[d_o] = press_o[il] + dcoef * (press_o[iu] - press_o[il]);
                    if (is_g)
                      {
                        gor = rs_reg[il] + dcoef * (rs_reg[iu] - rs_reg[il]);
                      }
                  }
                if (is_w)
                  p_p[d_w] = press_w[il] + dcoef * (press_w[iu] - press_w[il]);
                if (is_g)
                  p_p[d_g] = press_g[il] + dcoef * (press_g[iu] - press_g[il]);
              }

            //calc saturation
            if (n_phases > 1)
              {
                scal_prop->process_init (i_cell, p_p, i_sat, &perm[i_cell * PLANE_ORIENTATION_TOTAL], poro[i_cell], s_p, pc_limit);
              }

            if (is_o)
              {
                p_oil += p_p[d_o];
                if (is_g)
                  {
                    (*gas_oil_ratio)[i_cell] += gor;
                  }
              }
            if (is_w)
              p_water += p_p[d_w];
            if (is_g)
              p_gas += p_p[d_g];
            if (n_phases > 1)
              {
                if (is_w)
                  s_water += s_p[ds_w];
                if (is_g)
                  s_gas += s_p[ds_g];
              }
          }

        if (is_o)
          {
            p_oil /= (item_t)n_layer;
            if (is_g)
              {
                (*gas_oil_ratio)[i_cell] /= (item_t)n_layer;
              }
          }
        if (is_w)
          p_water /= (item_t)n_layer;
        if (is_g)
          p_gas /= (item_t)n_layer;
        if (n_phases > 1)
          {
            if (is_w)
              s_water /= (item_t)n_layer;
            if (is_g)
              s_gas /= (item_t)n_layer;
          }

        item_t smin  = 0, smax = 0, pmax_table = 0;
        if (n_phases > 1 && is_w && is_swatinit)
          {
            const scal_region &scal_water_data = scal_prop->get_water_data ()->get_region (i_sat);
            smin = scal_water_data.get_phase_sat_min ();
            smin = scal_prop->get_water_scale ()->get_sl (smin)[i_cell];
            smax = scal_water_data.get_phase_sat_max ();
            smax = scal_prop->get_water_scale ()->get_su (smax)[i_cell];
            pmax_table = scal_water_data.get_pcp_max ();
          }

        if (n_phases > 1 && is_w && is_swatinit && (p_oil - p_water > 0) && pmax_table > EPS_DIFF/*todo*/)
          {
            if ((*swatinit)[i_original_cell] < smin)
              s_water = smin;
            else if ((*swatinit)[i_original_cell] > smax)
              s_water = smax;
            else
              {
                s_water = (*swatinit)[i_original_cell];

                //calc pcw
                scal_prop->calc_pcp (i_cell, (*swatinit)[i_original_cell], i_sat, p_water - p_oil, pcw[i_cell]);
              }
            //recalc gas saturation
            if (n_phases == 3 && (s_water + s_gas > 1.0))
              s_gas = (1.0 - s_water > 0.) ? 1.0 - s_water : 0.;
          }
        //check sum
        else if (n_phases == 3 && (s_water + s_gas > 1.0))
          scal_prop->calc_gas_water_zone (i_cell, i_sat, &perm[i_cell * PLANE_ORIENTATION_TOTAL], poro[i_cell],
                                                  p_gas - p_water, s_water, s_gas);

        //check pressures
        if (n_phases > 1)
          {
            if (is_w && (p_oil - p_water) < pc_limit[ds_w])
              p_oil = p_water + pc_limit[ds_w];
            if (is_g && (p_gas - p_oil) > pc_limit[ds_g])
              p_oil = p_gas - pc_limit[ds_g];
          }

        //set pressure
        if (n_phases > 1 || is_o)
          (*pressure)[i_cell] = p_oil;
        else if (is_w)
          (*pressure)[i_cell] = p_water;
        else
          (*pressure)[i_cell] = p_gas;

        //set saturation
        if (n_phases > 1)
          {
            if (is_w)
              {
                if (fabs (s_water - (item_t) 1.0) < (item_t) 0.001)
                  {
                    s_water = 0.999f;
                  }
                (*saturation_3p)[n_phases * i_cell + d_w] = (float)s_water;
              }
            if (is_g)
              (*saturation_3p)[n_phases * i_cell + d_g] = (float)s_gas;

            (*saturation_3p)[n_phases * i_cell + d_o] = 1. - (float)s_water - (float)s_gas;
          }

        //set Rs by new pressure and set main_variable
        if (is_o && is_g)
          {
            if ((*saturation_3p)[n_phases * i_cell + d_g] < EPS_DIFF && (*saturation_3p)[n_phases * i_cell + d_o] < EPS_DIFF)
              {
                //gas_oil_ratio[i_cell] = 0;
                main_variable[i_cell] = FI_MOMG_VAR;
              }
            else
              {
                i_pvt = (*pvt_regions)[i_cell];
                //gas_oil_ratio[i_cell] = pvt_oil_array[i_pvt]->interpolate_and_fix (pressure[i_cell]);

                if ((*saturation_3p)[n_phases * i_cell + d_g] > EPS_DIFF)
                  main_variable[i_cell] = FI_SG_VAR;
                else
                  main_variable[i_cell] = FI_RO_VAR;
              }
          }
      }

    if (pcw.size ())
      scal_prop->get_water_scale ()->insert_pcp (pcw);

    return 0;
  }

  //template int
  //calc_model<base_strategy_di>::calc_equil (const sp_idata_t &input_data, const sp_mesh_iface_t &mesh);

  //template int
  //calc_model<base_strategy_fi>::calc_equil (const sp_idata_t &input_data, const sp_mesh_iface_t &mesh);

  //template int
  //calc_model<base_strategy_mixi>::calc_equil (const sp_idata_t &input_data, const sp_mesh_iface_t &mesh);

  //template int
  //calc_model<base_strategy_di>::equil_calc_pressure (item_t prev_press, item_t cur_d, item_t h, index_t phase, index_t i_pvt,
  //    double rs_type, item_t depth_goc, item_t rs_dat,
  //    val_vs_depth *rsvd, val_vs_depth *pbvd,
  //    item_t &p, item_t *rs);

  //template int
  //calc_model<base_strategy_fi>::equil_calc_pressure (item_t prev_press, item_t cur_d, item_t h, index_t phase, index_t i_pvt,
  //    double rs_type, item_t depth_goc, item_t rs_dat,
  //    val_vs_depth *rsvd, val_vs_depth *pbvd,
  //    item_t &p, item_t *rs);

  //template int
  //calc_model<base_strategy_mixi>::equil_calc_pressure (item_t prev_press, item_t cur_d, item_t h, index_t phase, index_t i_pvt,
  //    double rs_type, item_t depth_goc, item_t rs_dat,
  //    val_vs_depth *rsvd, val_vs_depth *pbvd,
  //    item_t &p, item_t *rs);


}//ns bs
