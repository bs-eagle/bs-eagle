/**
 *       \file  equil_model.cpp
 *      \brief  EQUIL Initialization model
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  03.05.2011
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "stdafx.h"
//#include "equil_model.hpp"
#include "equil_model_depth.h"
#include "rs_mesh_iface.h"
#include "calc_model.h"
#include "jfunction.h"

#define EQUIL_PI_METRIC 0.000096816838

namespace blue_sky 
{
  void
  equil_calc_pressure_2 (bool is_g, BS_SP (pvt_3p_iface) pvt_prop,
      t_double prev_press, t_double cur_d, t_double h, t_long phase, t_long i_pvt,
      double rs_type, t_double depth_goc, t_double rs_dat,
      val_vs_depth *rsvd, val_vs_depth *pbvd,
      /*output*/ t_double &p, t_double *rs = 0)
  {
    t_double p0, p1, rho, pb;
    t_double invers_fvf = 0;
    t_double rs_max;
    t_double gravity_constant = EQUIL_PI_METRIC;

    if (phase == FI_PHASE_OIL && is_g)
      {
        if (!rs)
          bs_throw_exception ("no rs");

        *rs = rs_dat;
      }

    p0 = 0.;
    p1 = prev_press;
    while (fabs (p1 - p0) > EPS_DIFF)
      {
        p0 = p1;

        if (phase == FI_PHASE_WATER)  //water phase
          {
            pvt_prop->get_pvt_water (i_pvt)->calc (0.5 * (prev_press + p0), &invers_fvf, 0, 0, 0, 0, 0);
            rho = pvt_prop->get_pvt_water (i_pvt)->get_surface_density () * invers_fvf;
          }
        else if (phase == FI_PHASE_GAS)  //gas phase
          {
            pvt_prop->get_pvt_gas (i_pvt)->calc (0.5 * (prev_press + p0), &invers_fvf, 0, 0, 0, 0, 0);
            rho = pvt_prop->get_pvt_gas (i_pvt)->get_surface_density () * invers_fvf;
          }
        else  //oil phase
          {
            pvt_prop->get_pvt_oil (i_pvt)->calc (is_g, !is_g || (cur_d - 0.5 * h < depth_goc) ? FI_SG_VAR : FI_RO_VAR,
                                        0.5 * (prev_press + p0), rs ? *rs : 0, &invers_fvf,
                                        0, 0, 0, 0, 0, 0, 0, 0,
                                        &rs_max, 0);

            rho = pvt_prop->get_pvt_oil (i_pvt)->get_surface_density () * invers_fvf;

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
                            *rs = pvt_prop->get_pvt_oil (i_pvt)->interpolate_and_fix (pb);
                          }
                        else
                          {
                            bs_throw_exception ("Error: you should specify either RSVD or PBVD keyword to use option 7 in keyword EQUIL");
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

                rho += *rs * pvt_prop->get_pvt_gas (i_pvt)->get_surface_density () * invers_fvf;
              }
          }

        p1 = prev_press + rho * gravity_constant * h;
      }

    p = (p0 + p1) * 0.5;
  }

  // initialize initial conditions
  void
  calc_equil_vs_depth (bool is_o, bool is_g, bool is_w,
                       const phase_d_t &phase_d, const sat_d_t &sat_d,
					             t_long n_phases, t_long phases, t_long sat_counter,
					             t_long n_depth,
                       //t_long ds_w, t_long ds_g,
                       const t_long n_eql,               //!< number of equil regions      
                       const stdv_long &sat_regions,      //!< (n_eql) scal region number for current equil region 
                       const stdv_long &pvt_regions,      //!< (n_eql) pvt region number for current equil region 
                       BS_SP (pvt_3p_iface) pvt_prop,    //!< PVT properties
                       BS_SP (scal_3p_iface) scal_prop,  //!< SCAL properties
                       spv_float equil,                 //!< main equil data (from EQUIL keyword)
                       idata::vval_vs_depth *rsvd,
                       idata::vval_vs_depth *pbvd,
                       const stdv_float &min_depth,     //!< n_eql
                       const stdv_float &max_depth,     //!< n_eql
                       const stdv_float &perm,          //!< (n_eql) permeability
                       const stdv_float &poro,          //!< (n_eql) poro
                       spv_double pressure,          //!< (n_phases * n_eql * n_depth)
                       spv_double saturation         //!< (n_phases * n_eql * n_depth)
                      )
  {
    const t_long n_layer = 10;
    t_long ds_w, ds_g;
    t_long i_depth;
    t_long i_eql, i_pvt, i_sat = 0;
    t_double d_depth, h;
    t_double depth_dat, press_dat, rs_dat, depth_woc, press_woc, depth_goc, press_goc;
    t_long il, iu, im;
    t_double prev_d, prev_p, cur_d;
    t_long main_phase;
    stdv_double depth, press, rs;      // array of depthes, phases pressure, gas oil ratios
    t_double p_oil, p_water, p_gas;
    t_double p_p[3], s_p[2], pc_limit[2];
	
    ds_w = sat_d[FI_PHASE_WATER];
    ds_g = sat_d[FI_PHASE_GAS];
    t_long d_o = phase_d[FI_PHASE_OIL];
    t_long d_w = phase_d[FI_PHASE_WATER];
    t_long d_g = phase_d[FI_PHASE_GAS];
    
    const t_long N_equil = 2;

    t_double s_water, s_gas;
    stdv_double press_w, press_g, press_o;
    spv_float pcw = BS_KERNEL.create_object (v_float::bs_type ());

    depth.resize(n_eql * n_depth, 0);
    press.resize(n_eql * n_phases * n_depth, 0);

    if (is_o && is_g)
      rs.resize(n_eql * n_depth, 0);

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////// CALCULATE Pp AND Rs BY DEPTH /////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    //loop through all eql regions
    for (i_eql = 0; i_eql < n_eql; i_eql++)
      {
        i_pvt = pvt_regions[i_eql];    // pvt region for current equil region

        depth_dat = (*equil)[EQUIL_TOTAL * i_eql + EQUIL_DAT_DEPTH];
        depth_woc = (*equil)[EQUIL_TOTAL * i_eql + EQUIL_WOC_DEPTH];
        depth_goc = (*equil)[EQUIL_TOTAL * i_eql + EQUIL_GOC_DEPTH];
        press_dat = (*equil)[EQUIL_TOTAL * i_eql + EQUIL_DAT_PRESS];

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
        d_depth = (max_depth[i_eql] - min_depth[i_eql]) / (t_double)(n_depth - 1);

        //set dep array values
        for (i_depth = 0; i_depth < n_depth; i_depth++)
          depth[i_depth + i_eql * n_depth] = min_depth[i_eql] + i_depth * d_depth;

        // calc Rs at datum pressure from PVTO table
        rs_dat = 0.;
        if (is_o && is_g)
          rs_dat = pvt_prop->get_pvt_oil (i_pvt)->interpolate_and_fix (press_dat);

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
        for (t_long up_down = -1; up_down < 2; up_down += 2)
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

                equil_calc_pressure_2 (is_g, pvt_prop, prev_p, cur_d, h, main_phase, i_pvt,
                                         is_g ? (*equil)[EQUIL_TOTAL * i_eql + EQUIL_RS_TYPE] : 0,
                                         is_g ? depth_goc : 0, is_g ? rs_dat : 0,
                                         (is_g && rsvd->size ()) ? &((*rsvd)[i_eql]) : 0, (is_g && pbvd->size ()) ? &(*pbvd)[i_eql] : 0,
                                         press[i_depth + (i_eql * n_phases + phase_d[main_phase]) * n_depth],
                                         rs.size () ? &rs[i_depth + i_eql * n_depth] : 0);

                prev_p = press[i_depth + (i_eql * n_phases + phase_d[main_phase]) * n_depth];
                prev_d = cur_d;

                i_depth += up_down;
              }
          }

        if (n_phases > 1)
          {
            t_double p_temp, d_temp, rs_temp;

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

                    equil_calc_pressure_2 (is_g, pvt_prop, prev_p, d_temp, h, main_phase, i_pvt,
                      0, 0, 0, 0, 0,
                      p_temp);
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
                  p_temp += (*equil)[EQUIL_TOTAL * i_eql + EQUIL_WOC_PRESS];
                else
                  p_temp -= (*equil)[EQUIL_TOTAL * i_eql + EQUIL_GOC_PRESS];

                //---- calc pressure table of oil -----
                for (t_long up_down = -1; up_down < 2; up_down += 2)
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

                        equil_calc_pressure_2 (is_g, pvt_prop, prev_p, cur_d, h, FI_PHASE_OIL, i_pvt,
                                                 is_g ? (*equil)[EQUIL_TOTAL * i_eql + EQUIL_RS_TYPE] : 0,
                                                 is_g ? depth_goc : 0, is_g ? rs_dat : 0,
                                                 (is_g && rsvd->size ()) ? &(*rsvd)[i_eql] : 0, (is_g && pbvd->size ()) ? &(*pbvd)[i_eql] : 0,
                                                 press[i_depth + (i_eql * n_phases + phase_d[FI_PHASE_OIL]) * n_depth],
                                                 rs.size () ? &rs[i_depth + i_eql * n_depth] : 0);

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
                    equil_calc_pressure_2 (is_g, pvt_prop, prev_p, depth_woc, h, FI_PHASE_OIL, i_pvt,
                                             is_g ? (*equil)[EQUIL_TOTAL * i_eql + EQUIL_RS_TYPE] : 0,
                                             is_g ? depth_goc : 0, is_g ? rs_dat : 0,
                                             (is_g && rsvd->size ()) ? &(*rsvd)[i_eql] : 0, (is_g && pbvd->size ()) ? &(*pbvd)[i_eql] : 0,
                                             press_woc, &rs_temp);
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
                press_woc -= (*equil)[EQUIL_TOTAL * i_eql + EQUIL_WOC_PRESS];

                //---- calc pressure table of water -----
                for (t_long up_down = -1; up_down < 2; up_down += 2)
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

                        equil_calc_pressure_2 (is_g, pvt_prop, prev_p, cur_d, h, FI_PHASE_WATER, i_pvt,
                                                 0, 0, 0, 0, 0,
                                                 press[i_depth + (i_eql * n_phases + phase_d[FI_PHASE_WATER]) * n_depth]);

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
                    equil_calc_pressure_2 (is_g, pvt_prop, prev_p, depth_goc, h, FI_PHASE_OIL, i_pvt,
                                             is_g ? (*equil)[EQUIL_TOTAL * i_eql + EQUIL_RS_TYPE] : 0,
                                             is_g ? depth_goc : 0, is_g ? rs_dat : 0,
                                             (is_g && rsvd->size ()) ? &(*rsvd)[i_eql] : 0, (is_g && pbvd->size ()) ? &(*pbvd)[i_eql] : 0,
                                             press_goc, &rs_temp);
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
                press_goc += (*equil)[EQUIL_TOTAL * i_eql + EQUIL_GOC_PRESS];

                //---- calc pressure table of gas -----
                for (t_long up_down = -1; up_down < 2; up_down += 2)
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

                        equil_calc_pressure_2 (is_g, pvt_prop, prev_p, cur_d, h, FI_PHASE_GAS, i_pvt,
                                                 0, 0, 0, 0, 0,
                                                 press[i_depth + (i_eql * n_phases + phase_d[FI_PHASE_GAS]) * n_depth]);

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
    pressure->resize (n_phases * n_eql * n_depth, 0);
    saturation->resize (n_phases * n_eql * n_depth, 0); 
   
	for (i_eql = 0; i_eql < n_eql; ++i_eql)
      {
        p_oil = 0.;
        p_water = 0.;
        p_gas = 0.;
        if (n_phases > 1 && is_w)
          pc_limit[ds_w] = 0.;
        if (n_phases > 1 && is_g)
          pc_limit[ds_g] = 0.;
      
        i_sat = sat_regions[i_eql];
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
        
        for (i_depth = 0; i_depth < n_depth; ++i_depth)
          {
            if (is_o)
              p_p[d_o] = p_oil = press_o[i_depth];
            if (is_w)
              p_p[d_w] = p_water = press_w[i_depth];
            if (is_g)
              p_p[d_g] = p_gas = press_g[i_depth];
        
            if (phases > 1)
              {
                scal_prop->process_init_2 (p_p, i_sat, perm[i_eql], poro[i_eql], s_p, pc_limit);
              }
            t_double s_water = s_p[ds_w];
			t_double s_gas = 0.0;
			if (is_g)
				s_gas = s_p[ds_g];

            if (n_phases == 3 && (s_water + s_gas > 1.0))
              scal_prop->calc_gas_water_zone_2 (i_sat, perm[i_eql], poro[i_eql],
                                              p_gas - p_water, s_water, s_gas);
       
/* 
            //check pressures
            if (n_phases > 1)
              {
                if (is_w && (p_oil - p_water) < pc_limit[ds_w])
                  p_oil = p_water + pc_limit[ds_w];
                if (is_g && (p_gas - p_oil) > pc_limit[ds_g])
                  p_oil = p_gas - pc_limit[ds_g];
              }
*/            
            
            //set pressure
            if (n_phases > 1 || is_o)
              (*pressure)[d_o * n_eql * n_depth + i_eql * n_depth + i_depth] = p_oil;
            if (is_w)
              (*pressure)[d_w * n_eql * n_depth + i_eql * n_depth + i_depth] = p_water;
            if (is_g)
              (*pressure)[d_g * n_eql * n_depth + i_eql * n_depth + i_depth] = p_gas;
        
            //set saturation
            if (n_phases > 1)
              {
                if (is_w)
                  {
                    if (fabs (s_water - (t_double) 1.0) < (t_double) 0.001)
                      {
                        s_water = 0.999f;
                      }
                  }
                if (is_g)
                  (*saturation)[d_g * n_eql * n_depth + i_eql * n_depth + i_depth] = (float)s_gas;
                if (is_w)
                  (*saturation)[d_w * n_eql * n_depth + i_eql * n_depth + i_depth] = (float)s_water;
                  
                (*saturation)[d_o * n_eql * n_depth + i_eql * n_depth + i_depth] = 1. - (float)s_water - (float)s_gas;
              }
          }    
      }
  }


  equil_model_depth::equil_model_depth(bs_type_ctor_param)
  {
    pressure = BS_KERNEL.create_object (v_double::bs_type ());
    saturation = BS_KERNEL.create_object (v_double::bs_type ());
  }

  equil_model_depth::equil_model_depth(const equil_model_depth &x)
	  :bs_refcounter(x),
	   pressure (BS_KERNEL.create_object (v_double::bs_type ())),
	   saturation (BS_KERNEL.create_object (v_double::bs_type ()))
  {
  }

  void
  equil_model_depth::py_calc_equil(bool is_o, bool is_g, bool is_w,
		        BS_SP(scal_3p_iface) scal_props, 
				BS_SP(pvt_3p_iface) pvt_props,
				spv_float equil,
				const stdv_float &min_depth,
				const stdv_float &max_depth,
				const stdv_float &perm,
				const stdv_float &poro,
				t_long n_depth,
				//stdv_double density,
				sp_jfunction jfunc_water,
				sp_jfunction jfunc_oil)
  {
	  t_long n_eql = 1;
	  stdv_long sat_regions, pvt_regions;
	  idata::vval_vs_depth val_dummy;
	  
	  phase_d_t phase_d;
	  sat_d_t sat_d;
	  t_long n_phases = 0;
      t_long phases = 0;
	  t_long sat_counter = 0;    
      if (is_w)
      {
        phase_d[0] = n_phases++;
        phases |= 1 << FI_PHASE_WATER;
      }
      else 
        phase_d[0] = -1;

      if (is_g)
      {
        phase_d[1] = n_phases++;
        phases |= 1 << FI_PHASE_GAS;
      }
      else 
        phase_d[1] = -1;

      if (is_o)
      {
        phase_d[2] = n_phases++;
        phases |= 1 << FI_PHASE_OIL;
      }
      else 
        phase_d[2] = -1;

      for (size_t i = 0, sat_counter = 0; i < FI_PHASE_TOT; ++i)
      {
        if ((phases & (1 << i)) && (sat_counter < n_phases ))
          sat_d[i] = sat_counter++;
        else
          sat_d[i] = -1;
      }

	  /*pvt_props->fill_pvt_arrays(is_o, is_g, is_w,
		                       1.0, 0.1, 1000.0, 100,
							   density);*/
	  //scal_props->init_from_scal();
	  sat_regions.push_back(0);
	  pvt_regions.push_back(0);

	  calc_equil_vs_depth(is_o, is_g, is_w, 
		                    phase_d, sat_d,
                        n_phases, phases, sat_counter,
                        n_depth,
                        n_eql, 
                        sat_regions, pvt_regions, pvt_props, scal_props,
                        equil, 
                        &val_dummy, &val_dummy, 
                        min_depth, max_depth, 
                        perm, poro, 
                        this->pressure, 
                        this->saturation);
  }

  spv_double
  equil_model_depth::get_pressure()
  {
	  return pressure;
  }

  spv_double
  equil_model_depth::get_saturation()
  {
	  return saturation;
  }

  BLUE_SKY_TYPE_STD_CREATE (equil_model_depth);
  BLUE_SKY_TYPE_STD_COPY (equil_model_depth);

  BLUE_SKY_TYPE_IMPL(equil_model_depth,  objbase, "equil_model_depth", "equil_model_depth", "equil_model_depth");

} //namespace