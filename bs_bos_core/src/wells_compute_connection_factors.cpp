/**
 *       \file  wells_compute_connection_factors.cpp
 *      \brief  Computes connection factors
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  08.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"

#include "calc_well.h"
#include "wells_compute_connection_factors.h"
#include "wells_common_const.h"
#include "fi_params.h"
#include "well_connection.h"

namespace blue_sky
  {
  namespace wells
    {
    namespace compute_factors
      {

      template <typename strategy_t>
      void
      peaceman_model<strategy_t>::compute (connection_t &con,
                                           const physical_constants &internal_constants,
                                           const sp_params_t &params,
                                           const sp_mesh_iface_t &mesh,
                                           const item_array_t &perm,
                                           const item_array_t &ntg,
                                           bool ro_calc_flag)
      {
        BS_ASSERT (con.n_block_ >= 0) (con.n_block_);

        int dir_orth = con.dir_;
        int dir_plane1 = (dir_orth + 1) % 3;
        int dir_plane2 = (dir_orth + 2) % 3;

        item_t d1, d2, d_ort;
        item_t d[] = {0, 0, 0};
        mesh->get_element_size (con.n_block (), d[0], d[1], d[2]);
        d1 = d[dir_plane1];
        d2 = d[dir_plane2];
        d_ort = d[dir_orth];

        item_t perm1 = (perm[3 * con.n_block () + dir_plane1]);
        item_t perm2 = (perm[3 * con.n_block () + dir_plane2]);
        item_t coef1 = fabs (perm1) < 1.0e-7 ? sqrt (perm2) : sqrt (perm2 / perm1);
        item_t coef2 = fabs (perm2) < 1.0e-7 ? sqrt (perm1) : sqrt (perm1 / perm2);

        if (con.R0_ < MIN_ZERO_DIFF || ro_calc_flag)
          {
            item_t denomin = sqrt (coef1) + sqrt (coef2);
            if (denomin > 1.0e-7)
              con.R0_ = 0.28 * sqrt (coef1 * d1 * d1 + coef2 * d2 * d2) / denomin;
            else
              {
                BS_ASSERT (denomin > 1.0e-7) (denomin) (coef1) (coef2);
                throw bs_exception ("well::connection::compute_factors_by_peaceman_model", "Connection factor is equal to zero");
              }
          }

        if (::log (2.0 * con.R0_ / con.diam_) + con.skin_ < MIN_ZERO_DIFF)
          {
            if (params->get_bool (fi_params::USE_LOW_SKIN_TRANS_MULT))
              {
                item_t D      = item_t (0.14036);
                item_t alpha  = d1 > d2 ? (d1 / d2) : (d2 / d1);
                item_t dP     = 1. / (2. * PI) * ::log (sqrt (1. / ((alpha * alpha + 1) * 4 * D * D)
                                                      + (alpha - 1) / (alpha + 1))+ 1. / (2 * D * sqrt (alpha * alpha + 1)));;

                if (con.skin_ < ::log (con.diam_ / (2 * con.R0_)) - 2 * PI * dP)
                  {
                    BS_ASSERT (con.skin_ >= ::log (con.diam_ / (2 * con.R0_)) - 2 * PI * dP) (::log (con.diam_ / (2 * con.R0_)) - 2 * PI * dP) (dP);
                    throw bs_exception ("well::connection::compute_factors_by_peaceman_model", "Wrong value of skin factor for connection");
                  }

                con.skin_ += 2 * PI * dP;
              }
            else
              {
                BS_ASSERT (false && "FI_PARAMS_B_USE_LOW_SKIN_TRANS_MULT");
                throw bs_exception ("well::connection::compute_factors_by_peaceman_model", "Wrong value of skin factor for connection");
              }
          }

        if (con.kh_ <= 0 || ro_calc_flag)
          {
            con.kh_ = sqrt (perm1 * perm2) * d_ort;
          }

        item_t cdarcy = 2 * PI * internal_constants.darcy_constant;
        con.fact_ = cdarcy * con.mult_ * con.kh_ / ((::log (2. * con.R0_ / con.diam_) + con.skin_));

        if (dir_orth == direction_z)
          {
            con.fact_ *= ntg[con.n_block_];

            if (con.connection_type_ != CONNECTION_USUAL)
              {
#ifndef _USE_GRP_FRACTURE_CORRELATION_
                con.fact_ *= compute_grp_pi_mult (con);
#else
                item_t fracture_mult_corr = fracture_mult_correlat (sqrt (d1 * d2), 2.0 * con.fracture_half_length_);
                if (fracture_mult_corr < MIN_FRACTURE_MULT_CORR)
                  {
                    BOSWARN (section::wells, level::warning) << "Couldn't correlate multiplier for fracture in cell for well" << bs_end;
                  }
                else
                  {
                    con.fact_ *= fracture_mult_corr;
                  }
#endif
              }
          }

        if (con.fact_ < -MIN_ZERO_DIFF)
          {
            bs_throw_exception (boost::format ("Negative value of connection factor %f for connection of well") % con.fact_);
          }
      }

      template <typename strategy_t>
      typename strategy_t::item_t
      peaceman_model<strategy_t>::compute_grp_pi_mult (connection_t &con)
      {
        return (::log (4.0 * con.fracture_half_length_ / con.diam_) /*+  con.skin_ */) / (::log (4.0) /*+ con.skin_ */);
      }

      //////////////////////////////////////////////////////////////////////////
      ///**
      // * \brief computes connection factor by baby and odeh model.
      // *          d1 equals to msh->elements[k].d[dir_plane]
      // *          d2 equals to msh->elements[k].d[dir_orth]
      // *          d3 equals to msh->elements[k].d[2]
      // * */
      //void
      //baby_odeh_model::compute (well::connection &con,
      //                          physical_constants *internal_constants,
      //                          item_t d1, item_t d2, item_t d3,
      //                          item_t perm1, item_t perm2, item_t perm3,
      //                          const well::item_array_t &ntg,
      //                          bool ro_calc_flag)
      //{
      //  item_t kr_z     = perm1;
      //  item_t k_orth   = perm2;
      //  item_t kr_plane = perm3;
      //  item_t c        = d1;
      //  item_t d        = d2;
      //  item_t h        = d3;

      //  item_t center_coordinate_xy = c * 0.5;
      //  item_t center_coordinate_z  = h * 0.5;

      //  item_t coef1    = sqrt (kr_z / kr_plane);
      //  item_t coef2    = c / h;
      //  item_t coef3    = center_coordinate_xy / c;

      //  // compute log (C_h)
      //  item_t log_Ch   = 6.28 * coef2 * coef1 * (1. / 3. - coef3 + coef3 * coef3)
      //                  - log (fabs (sin (PI * center_coordinate_z / h)))
      //                  - 0.5 * log (coef2 * coef1)
      //                  - 1.088;

      //  item_t full_skin  = 0.0;
      //  const item_t &L   = con.d_length_;
      //  item_t y_mid      = 0.0;

      //  if (L < DIFF_EPSILON)
      //    {
      //      con.fact_ = 0.0;
      //      return ;
      //    }

      //  // case 1: connection fully penetrates well grid block
      //  if (L >= d)
      //    {
      //      full_skin = con.skin_;
      //    }
      //  // case 2:
      //  else if (c / sqrt (kr_plane) >= 0.75 * d / sqrt (k_orth))
      //    {
      //      y_mid         = (con.con_end_distance_ + con.con_start_distance_) * 0.5;
      //      item_t p_xyz  = (d / L - 1.0) * (log (h / con.diam_) + 0.25 * log (kr_plane / kr_z)
      //                    - log (sin (PI * center_coordinate_z / h)) - 1.84);
      //      item_t p_xy   = (2.0 * d * d / (L * h)) * sqrt (kr_z / k_orth)
      //                    * (F_table (L, d, y_mid, 0) + 0.5 * (F_table (L, d, y_mid, 1) - F_table (L, d, y_mid, 2)));
      //      full_skin     = con.skin_ + p_xyz + p_xy;
      //    }
      //  // case 3:
      //  else
      //    {
      //      y_mid         = (con.con_end_distance_ + con.con_start_distance_) * 0.5;
      //      item_t p_xyz  = (d / L - 1.0) * (log (h / con.diam_) + 0.25 * log (kr_plane / kr_z)
      //                    - log (sin (PI * center_coordinate_z / h)) - 1.84);
      //      item_t p_y    = (6.28 * d * d / (c * h)) * sqrt (kr_z * kr_plane) / k_orth
      //                    * (1.0 / 3.0 - y_mid / d + y_mid * y_mid / (d * d) + L / (24.0 * d) * (L / d - 3.0));
      //      item_t p_xy   = (d / L - 1.0) * 6.28 * (c / h) * coef1 * (1.0 / 3.0 - coef3 + coef3 * coef3);
      //      full_skin     = con.skin_ + p_xyz + p_y + p_xy;
      //    }

      //  item_t cdarcy = 2 * PI * internal_constants->darcy_constant;
      //  con.fact_ = cdarcy * con.mult_ * sqrt (kr_z * kr_plane) * d
      //            / ((log (sqrt (c * h) / con.diam_) + log_Ch + full_skin - 0.75));
      //}

      ////! definition of constant #1
      //#define CONST_1 0.145
      ////! definition of constant #2
      //#define CONST_2 0.137
      ///*!
      //\brief prototype of table function
      //*/
      //inline item_t f_prototype (item_t x)
      //{
      //  return x * (CONST_1 + log (x) - CONST_2 * x * x);
      //}
      ///*!
      //\brief Function for horizontal wells for Babu model
      //\param L -- length of connection
      //\param d -- block size
      //\param y_mid -- middle point
      //\param arg_type -- type of arguments of function F
      //*/
      //item_t baby_odeh_model::F_table (item_t L, item_t d, item_t y_mid, int arg_type)
      //{
      //  item_t s;
      //  item_t s1;
      //  item_t s2;
      //
      //  s = L / (2. * d);
      //  s1 = (4. * y_mid + L) / (2. * d);
      //  s2 = (4. * y_mid - L) / (2. * d);
      //  item_t result = 0.0;
      //
      //  if (arg_type == 1)
      //    {
      //      if (s1 <= 1)
      //        {
      //          result = (-1.0) * f_prototype (s1);
      //        }
      //      else // s1 > 1
      //        {
      //          item_t v = 2.0 - s1;
      //          result = f_prototype (v);
      //        }
      //    }
      //  else if (arg_type == 2)
      //    {
      //      if (s2 <= 1)
      //        {
      //          result = (-1.0) * f_prototype (s2);
      //        }
      //      else
      //        {
      //          item_t v = 2.0 - s2;
      //          result = f_prototype (v);
      //        }
      //    }
      //  else // arg)type == 0
      //    {
      //      result = (-1.0) * f_prototype (s);
      //    }
      //  return result;
      //}

      template struct peaceman_model <base_strategy_fi>;
      template struct peaceman_model <base_strategy_di>;
      template struct peaceman_model <base_strategy_mixi>;

    } // namespace compute_factors
  } // namespace wells
} // namespace blue_sky

