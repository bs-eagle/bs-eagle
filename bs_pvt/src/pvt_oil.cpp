/**
 * \file pvt_oil.cpp
 * \brief pvt properties for oil phase
 * \author Miryanov Sergey
 * \date 06.05.2008
 */
#include "bs_pvt_stdafx.h"

#include "pvt_oil.h"
#include "pvt_interpolator.h"
#include "scal_interpolate.h"

using namespace blue_sky::pvt;

namespace blue_sky
  {

  pvt_oil::pvt_oil (bs_type_ctor_param)
  {
  }

  pvt_oil::pvt_oil (const pvt_oil &pvt)
  : bs_refcounter (pvt), pvt_dead_oil (pvt)
  {
    if (this != &pvt)
      {
        *this = pvt;
        BS_ASSERT (false && "NOT IMPL");
      }
  }

  void
  pvt_oil::insert_vector (const v_double &vec)
  {
    const int elem_count = 4;
    BS_ASSERT (!(vec.size() % elem_count)) (vec.size ()) (elem_count);

    if (!this->init_dependent)
      {
        //main_gpr_.clear ();
        //main_pressure_.clear ();
        //main_fvf_.clear ();
        //main_visc_.clear ();

        base_t::init_dependent = true;
      }

    t_int n_points = (t_int) vec.size () / elem_count;
    
    if (pvt_input_props->init (n_points, PVT_OIL_INPUT_TOTAL))
      {
        throw bs_exception ("pvt_oil::insert_vector in table", "Error: initializing table of properties");
      }

    pvt_input_props->set_col_name (PVT_OIL_INPUT_GPR, "gor");
    pvt_input_props->set_col_name (PVT_OIL_INPUT_PRESSURE, "pressure");
    pvt_input_props->set_col_name (PVT_OIL_INPUT_FVF, "fvf");   
    pvt_input_props->set_col_name (PVT_OIL_INPUT_VISC, "visc");

    vector_t &main_gpr_          = pvt_input_props->get_col_vector (PVT_OIL_INPUT_GPR);
    vector_t &main_pressure_     = pvt_input_props->get_col_vector (PVT_OIL_INPUT_PRESSURE);
    vector_t &main_fvf_          = pvt_input_props->get_col_vector (PVT_OIL_INPUT_FVF);
    vector_t &main_visc_         = pvt_input_props->get_col_vector (PVT_OIL_INPUT_VISC);

    for (t_int i = 0; i < n_points; ++i)
      {
        main_gpr_[i]      = (vec[i * elem_count + 0]);
        main_pressure_[i] = (vec[i * elem_count + 1]);
        main_fvf_[i]      = (vec[i * elem_count + 2]);
        main_visc_[i]     = (vec[i * elem_count + 3]);
      }
  }

  void
  pvt_oil::build (t_double atm_p, t_double min_p, t_double max_p, t_long n_intervals)
  {
    int n_points = base_t::build_internal (atm_p, min_p, max_p, n_intervals, true);

    build_compressibility (n_points);
  }

  void
  pvt_oil::build_compressibility (t_long n_intervals)
  {
    vector_t &main_gpr_          = pvt_input_props->get_col_vector (PVT_OIL_INPUT_GPR);
    vector_t &main_pressure_     = pvt_input_props->get_col_vector (PVT_OIL_INPUT_PRESSURE);
    vector_t &main_fvf_          = pvt_input_props->get_col_vector (PVT_OIL_INPUT_FVF);
    vector_t &main_visc_         = pvt_input_props->get_col_vector (PVT_OIL_INPUT_VISC);
  
    //vector_t &pressure_     = pvt_props_table->get_col_vector (PVT_OIL_PRESSURE);
    //vector_t &inv_fvf_      = pvt_props_table->get_col_vector (PVT_OIL_INV_FVF);
    //vector_t &inv_visc_     = pvt_props_table->get_col_vector (PVT_OIL_INV_VISC);
    vector_t &gor_          = pvt_props_table->get_col_vector (PVT_OIL_GOR);

    if (main_pressure_.empty ())
      return ;


    if (pvt_props_table->add_col_vector ("compress_fvf") != PVT_OIL_COMPRESS_FVF) 
      {
        BS_ASSERT (false);
        throw bs_exception ("pvt_oil::build_compressibility", "Error: initializing vector of fvf compressibility");
        return;
      }  
    if (pvt_props_table->add_col_vector ("compress_visc") != PVT_OIL_COMPRESS_VISC)
      {
        BS_ASSERT (false);
        throw bs_exception ("pvt_oil::build_compressibility", "Error: initializing vector of viscosity compressibility");
        return;
      }  
    
    vector_t &compress_fvf_     = pvt_props_table->get_col_vector (PVT_OIL_COMPRESS_FVF);
    vector_t &compress_visc_    = pvt_props_table->get_col_vector (PVT_OIL_COMPRESS_VISC);
    
    //compress_fvf_.resize (n_intervals);
    //compress_visc_.resize (n_intervals);
    for (t_long i = 0; i < n_intervals; ++i)
      {
        t_long j1, j2, end_j1, end_j2;
        get_compressibility_interval (gor_[i], j1, j2, end_j1, end_j2);

        t_double diff = 0;
        if (j1 < 0)
          {
            diff              = 1.0 / (main_pressure_[end_j2] - main_pressure_[j2]);
            compress_fvf_[i]  = pvt::wtf (main_fvf_[j2], main_fvf_[end_j2]) * 2.0 * diff;
            compress_visc_[i] = pvt::wtf (main_visc_[j2], main_visc_[end_j2]) * 2.0 * diff;
          }
        else if (j2 < 0)
          {
            diff                = 1.0 / (main_pressure_[end_j1] - main_pressure_[j1]);
            compress_fvf_[i]    = pvt::wtf (main_fvf_[j1], main_fvf_[end_j1]) * 2.0 * diff;
            compress_visc_[i]   = pvt::wtf (main_visc_[j1], main_visc_[end_j1]) * 2.0 * diff;
          }
        else
          {
            diff            = 1.0 / (main_pressure_[end_j1] - main_pressure_[j1]);
            t_double c_fvf_1  = pvt::wtf (main_fvf_[j1], main_fvf_[end_j1]) * 2.0 * diff;
            t_double c_visc_1 = pvt::wtf (main_visc_[j1], main_visc_[end_j1]) * 2.0 * diff;

            diff            = 1.0 / (main_pressure_[end_j2] - main_pressure_[j2]);
            t_double c_fvf_2  = pvt::wtf (main_fvf_[j2], main_fvf_[end_j2]) * 2.0 * diff;
            t_double c_visc_2 = pvt::wtf (main_visc_[j2], main_visc_[end_j2]) * 2.0 * diff;

            if (fabs (main_gpr_[j2] - main_gpr_[j1]) >  EPS_DIFF)
              {
                compress_fvf_[i]  = pvt::interpolate (gor_[i], main_gpr_[j1], main_gpr_[j2], c_fvf_1, c_fvf_2, diff);
                compress_visc_[i] = pvt::interpolate (c_visc_1, c_visc_2, diff);
              }
            else
              {
                compress_fvf_[i]  = c_fvf_1;
                compress_visc_[i] = c_visc_1;
              }
          }
      }
  }



  void
  pvt_oil::get_compressibility_interval (t_double gor, t_long &j1, t_long &j2, t_long &end_j1, t_long &end_j2)
  {
    vector_t &main_gpr_          = pvt_input_props->get_col_vector (PVT_OIL_INPUT_GPR);
    vector_t &main_pressure_     = pvt_input_props->get_col_vector (PVT_OIL_INPUT_PRESSURE);
  
    j1 = j2 = -1;
    end_j1 = end_j2 = -2;
    for (t_long i = 0, cnt = (t_long)main_pressure_.size () ; i < cnt - 1; i++)
      {
        if (main_gpr_[i] > -0.5 && main_gpr_[i + 1] < -0.5)
          {
            // find end of slop
            int end_j = i + 1;
            for (; end_j < cnt && main_gpr_[end_j] < -0.5; ++end_j)
              ;
            if (end_j >= cnt || main_gpr_[end_j] > -0.5)
              --end_j;

            if (main_gpr_[i] <= gor)
              {
                j1 = i;
                end_j1 = end_j;
              }
            if (main_gpr_[i] >= gor)
              {
                j2 = i;
                end_j2 = end_j;

                break;
              }
          }
      }

    if (j1 < 0 && j2 < 0)
      throw bs_exception ("", "invalid interval");
  }

  void
  pvt_oil::check_oil ()
  {
    vector_t &main_gpr_          = pvt_input_props->get_col_vector (PVT_OIL_INPUT_GPR);
    vector_t &main_pressure_     = pvt_input_props->get_col_vector (PVT_OIL_INPUT_PRESSURE);
    vector_t &main_fvf_          = pvt_input_props->get_col_vector (PVT_OIL_INPUT_FVF);
    vector_t &main_visc_         = pvt_input_props->get_col_vector (PVT_OIL_INPUT_VISC);
  
    BS_ASSERT (main_pressure_.size ());

    base_t::check_common ();

    check_oil_common (main_pressure_, main_fvf_, main_visc_);

    if (main_gpr_.front () < -0.5)
      {
        // TODO: LOG
        BS_ASSERT (false);
        throw bs_exception ("", "first point in PVTO should be a saturated point");
      }

    t_long i_prev_sat = 0;
    t_long i_prev_unsat = 0;
    t_long state = 0;
    for (t_long i = 1, cnt = (t_long)main_pressure_.size (); i < cnt; ++i)
      {
        if (main_gpr_[i] > -0.5 && state == 1)
          {
            state = 0;
          }
        else if (main_gpr_[i] < -0.5 && state == 0)
          {
            state = 1;
            i_prev_unsat = i_prev_sat;
          }

        //for saturated curve
        if (state == 0)
          {
            if (main_pressure_[i] - main_pressure_[i_prev_sat] <= EPS_DIFF)
              {
                // TODO: LOG
                BS_ASSERT (false);
                throw bs_exception ("", "pressure curve in PVTO should be monotonically increasing function");
              }
            if (main_gpr_[i] - main_gpr_[i_prev_sat] <= EPS_DIFF)
              {
                // TODO: LOG
                BS_ASSERT (false);
                throw bs_exception ("", "Gas oil ratio curve in PVTO should be monotonically increasing function");
              }
            if (main_fvf_[i] - main_fvf_[i_prev_sat] <= 0)
              {
                // TODO: LOG
                BS_ASSERT (false);
                throw bs_exception ("", "FVF curve in PVTO should be monotonically increasing function for saturated oil");
              }
            if (main_visc_[i] - main_visc_[i_prev_sat] >= 0)
              {
                // TODO: LOG
                BS_ASSERT (false);
                throw bs_exception ("", "Viscosity curve in PVTO should be monotonically decreasing function for saturated oil");
              }

            i_prev_sat = i;
          }
        else
          {
            if (main_pressure_[i] - main_pressure_[i_prev_unsat] <= EPS_DIFF)
              {
                // TODO: LOG
                BS_ASSERT (false);
                throw bs_exception ("", "pressure curve in PVTO should be monotonically increasing function");
              }
            if (main_fvf_[i] - main_fvf_[i_prev_unsat] >= 0)
              {
                // TODO: LOG
                BS_ASSERT (false);
                throw bs_exception ("", "FVF curve in PVTO should be monotonically decreasing function for undersaturated oil");
              }
            if (main_visc_[i] - main_visc_[i_prev_unsat] <= 0)
              {
                // TODO: LOG
                BS_ASSERT (false);
                throw bs_exception ("", "Viscosity curve in PVTO should be monotonically increasing function for undersaturated oil");
              }

            i_prev_unsat = i;
          }
      }
  }

  bool
  pvt_oil::calc (const bool is_g, const int main_var, const t_double p, const t_double gor,
                             t_double *inv_fvf, t_double *d_inv_fvf, t_double *gor_d_inv_fvf,
                             t_double *inv_visc, t_double *d_inv_visc, t_double *gor_d_inv_visc,
                             t_double *inv_visc_fvf, t_double *d_inv_visc_fvf, t_double *gor_d_inv_visc_fvf,
                             t_double *gas_oil_ratio, t_double *d_gas_oil_ratio,
                             const t_double drsdt /* = -1.0 */, const t_double dt /* = 0 */,
                             const t_double old_gas_oil_ratio /* = 0 */) const
    {
      if (is_g && main_var == FI_RO_VAR)
        {
          return calc_undersaturated_oil (p, gor, inv_fvf, d_inv_fvf, gor_d_inv_fvf,
                                          inv_visc, d_inv_visc, gor_d_inv_visc,
                                          inv_visc_fvf, d_inv_visc_fvf, gor_d_inv_visc_fvf,
                                          gas_oil_ratio, d_gas_oil_ratio, drsdt, dt, old_gas_oil_ratio);
        }
      else
        {
          return base_t::calc_saturated_oil (is_g, main_var, p, gor, inv_fvf, d_inv_fvf, gor_d_inv_fvf,
                                             inv_visc, d_inv_visc, gor_d_inv_visc,
                                             inv_visc_fvf, d_inv_visc_fvf, gor_d_inv_visc_fvf,
                                             gas_oil_ratio, d_gas_oil_ratio, drsdt, dt, old_gas_oil_ratio);
        }
    }

  inline bool
  pvt_oil::calc_undersaturated_oil (const t_double p, const t_double gor,
      t_double *inv_fvf, t_double *d_inv_fvf, t_double *gor_d_inv_fvf,
      t_double *inv_visc, t_double *d_inv_visc, t_double *gor_d_inv_visc,
      t_double *inv_visc_fvf, t_double *d_inv_visc_fvf, t_double *gor_d_inv_visc_fvf,
      t_double *gas_oil_ratio, t_double *d_gas_oil_ratio,
      const t_double /*drsdt = -1.0 */, const t_double /*dt  = 0 */,
      const t_double /*old_gas_oil_ratio = 0 */) const
    {
      vector_t &pressure_     = pvt_props_table->get_col_vector (PVT_OIL_PRESSURE);
      vector_t &inv_fvf_      = pvt_props_table->get_col_vector (PVT_OIL_INV_FVF);
      vector_t &inv_visc_     = pvt_props_table->get_col_vector (PVT_OIL_INV_VISC);
      vector_t &gor_          = pvt_props_table->get_col_vector (PVT_OIL_GOR);
      vector_t &compress_fvf_     = pvt_props_table->get_col_vector (PVT_OIL_COMPRESS_FVF);
      vector_t &compress_visc_    = pvt_props_table->get_col_vector (PVT_OIL_COMPRESS_VISC);
      
      size_t i = binary_search (gor, gor_, std::less <t_double> ());
      if (i == 0)
        i++;

      BS_ASSERT (i < (size_t)pressure_.size ()) (i) (pressure_.size ());
      BS_ASSERT (i < (size_t)compress_fvf_.size ()) (i) (compress_fvf_.size ());

      if (i >= pressure_.size ())
        {
          BOSERR (section::pvt, level::error) << "pvt_oil::calc_undersaturated_oil: index_out_of_range (group): " << i << " [" << gor << "]" << bs_end;
          // TODO: BUG:
          return false;
        }
      if (i >= compress_fvf_.size ())
        {
          BOSERR (section::pvt, level::error) << "pvt_oil::calc_undersaturated_oil: index_out_of_range (comp): " << i << " [" << gor << "]" << bs_end;
          // TODO: BUG:
          return false;
        }

      t_double dp             = pressure_[i] - pressure_[i - 1];
      t_double dgor           = (gor_[i] - gor_[i - 1]);

      t_double pbub_, d_pbub, d_gor_ifvf_, d_gor_ivisc_, d_cfvf_, d_cvisc_, cfvf_, cvisc_, ifvf_ro_, ivisc_ro_;

      if (fabs (dgor) < EPS_DIFF)
        {
          d_pbub        = 0;
          d_gor_ifvf_   = 0;
          d_gor_ivisc_  = 0;
          d_cfvf_       = 0;
          d_cvisc_      = 0;

          cfvf_         = compress_fvf_[i - 1];
          cvisc_        = compress_visc_[i - 1];
          pbub_         = pressure_[i - 1];
          ifvf_ro_      = inv_fvf_[i - 1];
          ivisc_ro_     = inv_visc_[i - 1];
        }
      else
        {
          t_double idgor  = 1.0 / dgor;
          d_pbub        = dp * idgor;

          d_gor_ifvf_   = (inv_fvf_[i] - inv_fvf_[i - 1]) * idgor;
          d_gor_ivisc_  = (inv_visc_[i] - inv_visc_[i - 1]) * idgor;
          d_cfvf_       = (compress_fvf_[i] - compress_fvf_[i - 1]) * idgor;
          d_cvisc_      = (compress_visc_[i] - compress_visc_[i - 1]) * idgor;

          t_double diff_gor = gor - gor_[i - 1];

          cfvf_         = compress_fvf_[i - 1] + d_cfvf_ * diff_gor;
          cvisc_        = compress_visc_[i - 1] + d_cvisc_ * diff_gor;
          pbub_         = pressure_[i - 1] + d_pbub * diff_gor;

          ifvf_ro_      = (inv_fvf_[i - 1] + d_gor_ifvf_ * diff_gor);     // FVF for saturated oil
          ivisc_ro_     = (inv_visc_[i - 1] + d_gor_ivisc_ * diff_gor);   // VISC for saturated oil
        }

      t_double exp_fvf_   = 1 - cfvf_ * (p - pbub_);                      // FVF for undesaturated oil
      t_double ifvf_      = ifvf_ro_ * exp_fvf_;

      t_double exp_visc_  = exp (-cvisc_ * (p - pbub_));                 // VISC for undesaturated oil
      t_double ivisc_     = ivisc_ro_ * exp_visc_;

      set_pvt_pointer (inv_fvf, ifvf_);
      set_pvt_pointer (inv_visc, ivisc_);
      set_pvt_pointer (inv_visc_fvf, ifvf_ * ivisc_);

      //if (d_inv_fvf || d_inv_visc || d_inv_visc_fvf)
        {
          t_double d_inv_fvf_   = -cfvf_  * ifvf_ro_;
          //t_double d_inv_visc_  = -cvisc_ * ivisc_ro_;
          t_double d_inv_visc_  = -cvisc_ * ivisc_;

          set_pvt_pointer (d_inv_fvf, d_inv_fvf_);
          set_pvt_pointer (d_inv_visc, d_inv_visc_);
          set_pvt_pointer (d_inv_visc_fvf, d_inv_fvf_ * ivisc_ + d_inv_visc_ * ifvf_);
        }

      //if (gor_d_inv_fvf || gor_d_inv_visc || gor_d_inv_visc_fvf)
        {
          t_double diff_pbub = p - pbub_;
          d_gor_ifvf_   = d_gor_ifvf_ * (1 - cfvf_ * diff_pbub) - ifvf_ro_ * (diff_pbub * d_cfvf_ - cfvf_ * d_pbub);
          d_gor_ivisc_  = d_gor_ivisc_ * exp_visc_ + ivisc_ * (- d_cvisc_ * diff_pbub + cvisc_ * d_pbub);

          set_pvt_pointer (gor_d_inv_fvf, d_gor_ifvf_);
          set_pvt_pointer (gor_d_inv_visc, d_gor_ivisc_);
          set_pvt_pointer (gor_d_inv_visc_fvf, d_gor_ifvf_ * ivisc_ + d_gor_ivisc_ * ifvf_);
        }

      set_pvt_pointer (gas_oil_ratio, gor);
      set_pvt_pointer (d_gas_oil_ratio, 0);

      return true;
    }

  void
  pvt_oil::print () const
  {
    vector_t &pressure_     = pvt_props_table->get_col_vector (PVT_OIL_PRESSURE);
    vector_t &inv_fvf_      = pvt_props_table->get_col_vector (PVT_OIL_INV_FVF);
    vector_t &inv_visc_     = pvt_props_table->get_col_vector (PVT_OIL_INV_VISC);
    vector_t &inv_visc_fvf_ = pvt_props_table->get_col_vector (PVT_OIL_INV_VISC_FVF);
    vector_t &gor_          = pvt_props_table->get_col_vector (PVT_OIL_GOR);
    vector_t &compress_fvf_     = pvt_props_table->get_col_vector (PVT_OIL_COMPRESS_FVF);
    vector_t &compress_visc_    = pvt_props_table->get_col_vector (PVT_OIL_COMPRESS_VISC);

    BS_ASSERT (pressure_.size () == inv_fvf_.size ());
    BS_ASSERT (inv_fvf_.size ()  == inv_visc_.size ());
    BS_ASSERT (inv_visc_.size () == inv_visc_fvf_.size ());
    BS_ASSERT (inv_visc_fvf_.size () == compress_fvf_.size ());
    BS_ASSERT (compress_fvf_.size () == compress_visc_.size ());
    BS_ASSERT (compress_visc_.size () == gor_.size ());
    
    BOSOUT (section::pvt, level::medium) << bs_end;
    BOSOUT (section::pvt, level::medium) << "****************************************************** PVTO **************************************************************" << bs_end;
    BOSOUT (section::pvt, level::medium) 
      << boost::format ("*%9s%13s%13s%11s%12s%12s%12s%12s%12s%12s%2s\n")
      % "P"
      % "RS" 
      % "FVF"
      % "1/FVF"
      % "VISC"
      % "1/VISC"
      % "VISC*FVF"
      % "1/FVF/VISC"
      % "C_FVF"
      % "C_VISC"
      % "*"
      << bs_end;
    BOSOUT (section::pvt, level::medium) << "*************************************************************************************************************************" << bs_end;
    for (size_t j = 0, jcnt = pressure_.size ();j < jcnt; ++j)
      {
        BOSOUT (section::pvt, level::medium)
          << boost::format (
            "*%9.6f\t\t%9.6f\t\t%9.6f\t\t%9.6f\t\t%9.6f\t\t%9.6f\t\t%9.6f\t\t%9.6f\t\t%9.6f\t\t%9.6f%2s"
          )
          % (pressure_[j])
          % (gor_[j])
          % (1.0 / inv_fvf_[j])
          % (inv_fvf_[j])
          % (1.0 / inv_visc_[j])
          % (inv_visc_[j])
          % (1.0 / inv_visc_fvf_[j])
          % (inv_visc_fvf_[j])
          % (compress_fvf_[j])
          % (compress_visc_[j])
          % "*"
          << bs_end;
      }
    BOSOUT (section::pvt, level::medium) << "*************************************************************************************************************************" << bs_end;
  }

  BLUE_SKY_TYPE_STD_CREATE (pvt_oil);
  BLUE_SKY_TYPE_STD_COPY (pvt_oil);

  BLUE_SKY_TYPE_IMPL(pvt_oil,  pvt_dead_oil, "pvt_oil", "Oil PVT calculation class", "Oil PVT calculation");

} // namespace blue_sky

