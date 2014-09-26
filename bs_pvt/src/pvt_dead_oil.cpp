/**
 * \file pvt_dead_oil.cpp
 * \brief
 * \author Miryanov Sergey
 * \date 07.05.2008
 */

#include "bs_common.h"
#include BS_FORCE_PLUGIN_IMPORT()
#include "bos_report.h"
#include BS_STOP_PLUGIN_IMPORT()

#include "pvt_dead_oil.h"
#include "pvt_interpolator.h"
#include "scal_interpolate.h"

using namespace blue_sky::pvt;

namespace blue_sky
  {

  pvt_dead_oil::pvt_dead_oil (bs_type_ctor_param)
  {
    if (pvt_input_props->init (0, PVT_OIL_INPUT_TOTAL))
      {
        bs_throw_exception ("Error: initializing table of properties");
      }
    if (pvt_props_table->init (0, PVT_OIL_TOTAL))
      {
        bs_throw_exception ("Error: initializing table of properties");
      }
  }

  pvt_dead_oil::pvt_dead_oil (const pvt_dead_oil &pvt)
  : bs_refcounter (pvt)
  {
    if (this != &pvt)
      {
        // TODO: LOG
        *this = pvt;
        BS_ASSERT (false && "NOT IMPL");
      }
  }


  void
  pvt_dead_oil::insert_vector (const v_double &vec)
  {
    const int elem_count = 3;
    BS_ASSERT (!(vec.size() % elem_count)) (vec.size ()) (elem_count);

    if (!base_t::init_dependent)
      {
        base_t::init_dependent = true;

        //pvt_input_props->clear ();

        //main_gpr_.clear ();
        //main_pressure_.clear ();
        //main_fvf_.clear ();
        //main_visc_.clear ();
      }
    t_int n_points = (t_int) vec.size () / elem_count;

    if (pvt_input_props->init (n_points, PVT_OIL_INPUT_TOTAL))
      {
        throw bs_exception ("pvt_dead_oil::insert_vector in table", "Error: initializing table of properties");
      }

    int n_cols = pvt_input_props->get_n_cols ();
    //n_cols==3 -- 2 phase model, 3 columns in PVT: pressure, fvf, viscosity
    int shift = 0;
    if (n_cols == 3)
      shift = 1;

    pvt_input_props->set_col_name (PVT_OIL_INPUT_GPR, L"gor");
    pvt_input_props->set_col_name (PVT_OIL_INPUT_PRESSURE - shift, L"pressure");
    pvt_input_props->set_col_name (PVT_OIL_INPUT_FVF - shift, L"fvf");
    pvt_input_props->set_col_name (PVT_OIL_INPUT_VISC - shift, L"visc");

    vector_t &main_gpr_        = pvt_input_props->get_col_vector (PVT_OIL_INPUT_GPR);
    vector_t &main_pressure_     = pvt_input_props->get_col_vector (PVT_OIL_INPUT_PRESSURE - shift);
    vector_t &main_fvf_          = pvt_input_props->get_col_vector (PVT_OIL_INPUT_FVF - shift);
    vector_t &main_visc_         = pvt_input_props->get_col_vector (PVT_OIL_INPUT_VISC - shift);

    for (t_int i = 0; i < n_points; ++i)
      {
        main_gpr_[i]      = (0.0);
        main_pressure_[i] = (vec[i * elem_count + 0]);
        main_fvf_[i]      = (vec[i * elem_count + 1]);
        main_visc_[i]     = (vec[i * elem_count + 2]);
      }
  }

  int
  pvt_dead_oil::build_internal (t_double /*atm_p*/, t_double min_p, t_double max_p, t_long /*n_intervals*/, bool is_pvto)
  {
    int n_cols = pvt_input_props->get_n_cols ();
    int shift = 0;
    if (n_cols == 3)
      shift = 1;
    vector_t &main_gpr_          = pvt_input_props->get_col_vector (PVT_OIL_INPUT_GPR);
    vector_t &main_pressure_     = pvt_input_props->get_col_vector (PVT_OIL_INPUT_PRESSURE - shift);
    vector_t &main_fvf_          = pvt_input_props->get_col_vector (PVT_OIL_INPUT_FVF - shift);
    vector_t &main_visc_         = pvt_input_props->get_col_vector (PVT_OIL_INPUT_VISC - shift);

    check_oil ();

    if (main_pressure_.empty ())
      return 0;

    t_long n_points = (t_long)main_pressure_.size ();
    if (is_pvto)
      {
        n_points = 0;
        for (t_long i = 0, cnt = (t_long)main_gpr_.size (); i < cnt; ++i)
          {
            if (main_gpr_[i] > -0.5)
              ++n_points;
          }
      }

    bool is_min         = false,      is_max    = false;
    t_long i_min       = 0,          i_max     = (t_long)main_pressure_.size () - 1;
    t_long i_next_min  = i_min + 1,  i_pre_max = i_max - 1;

    if (is_pvto)
      {
        while (main_gpr_[i_min] < 0.)
          ++i_min;

        if (min_p < main_pressure_[i_min])
          {
            is_min = true;
            n_points += (is_min != false);

            i_next_min = i_min + 1;
            while (main_gpr_[i_next_min] < 0.)
              ++i_next_min;

            if (i_next_min > (t_long)main_pressure_.size () - 1)
              {
                BS_ASSERT (i_next_min <= (t_long)main_pressure_.size () - 1) (i_next_min) (main_pressure_.size ());
                throw bs_exception ("pvt_dead_oil::build", "Error: number of Rs values in PVTO should be greater than 1");
              }
          }

        while (main_gpr_[i_max] < 0.)
          --i_max;

        if (max_p > main_pressure_[i_max])
          {
            is_max = true;
            n_points += (is_max != false);

            i_pre_max = i_max - 1;
            while (main_gpr_[i_pre_max] < 0.)
              --i_pre_max;

            if (i_pre_max < 0)
              {
                BS_ASSERT (i_pre_max >= 0) (i_pre_max);
                throw bs_exception ("pvt_dead_oil::build", "Error: number of Rs values in PVTO should be greater than 1");
              }
          }
      }
    else
      {
        is_min = min_p < main_pressure_.front ();
        is_max = max_p > main_pressure_.back ();
        n_points += (is_min != false);
        n_points += (is_max != false);
      }

    if (pvt_props_table->init (n_points, PVT_OIL_TOTAL))
      {
        throw bs_exception ("pvt_dead_oil::init table", "Error: initializing table of properties");
      }

    pvt_props_table->set_col_name (PVT_OIL_PRESSURE, L"pressure");
    pvt_props_table->set_col_name (PVT_OIL_INV_FVF, L"inv_fvf");
    pvt_props_table->set_col_name (PVT_OIL_INV_VISC, L"inv_visc");
    pvt_props_table->set_col_name (PVT_OIL_INV_VISC_FVF, L"inv_visc_fvf");
    pvt_props_table->set_col_name (PVT_OIL_GOR, L"gor");

    vector_t &pressure_     = pvt_props_table->get_col_vector (PVT_OIL_PRESSURE);
    vector_t &inv_fvf_      = pvt_props_table->get_col_vector (PVT_OIL_INV_FVF);
    vector_t &inv_visc_     = pvt_props_table->get_col_vector (PVT_OIL_INV_VISC);
    vector_t &inv_visc_fvf_ = pvt_props_table->get_col_vector (PVT_OIL_INV_VISC_FVF);
    vector_t &gor_          = pvt_props_table->get_col_vector (PVT_OIL_GOR);

    if (main_pressure_.size () < 2)
      {
        if (is_pvto)
          {
            BS_ASSERT (main_pressure_.size () >= 2) (main_pressure_.size ());
            throw bs_exception ("pvt_dead_oil::build", "Error: number of rows in PVTO should be greater than 1");
          }

        for (t_long i = 0; i < n_points; ++i)
          {
            if (i == 0 && is_min)
              pressure_.front () = min_p;
            else if (i == n_points - 1 && is_max)
              pressure_.back ()  = max_p;
            else
              pressure_[i] = main_pressure_.front ();

            gor_[i]           = 0.0;
            inv_fvf_[i]       = 1.0 / main_fvf_.front ();
            inv_visc_[i]      = 1.0 / main_visc_.front ();
            inv_visc_fvf_[i]  = inv_fvf_[i] * inv_visc_[i];
          }
      }
    else
      {
        t_long i_rs = 0;
        for (t_long i = 0; i < n_points; ++i)
          {
            if ((i == 0 && is_min) || (i == n_points - 1 && is_max))
              {
                pressure_[i] = i == 0 ? min_p : max_p;

                t_double diff = 0;

                if (i == 0)
                  {
                    diff = (min_p - main_pressure_[i_min]) / (main_pressure_[i_next_min] - main_pressure_[i_min]);
                    inv_fvf_[i] = 1.0 / (main_fvf_[i_min] + (main_fvf_[i_next_min] - main_fvf_[i_min]) * diff);
                    inv_visc_[i] = 1.0 / (main_visc_[i_min] + (main_visc_[i_next_min] - main_visc_[i_min]) * diff);

                    if (is_pvto)
                      gor_[i] = main_gpr_[i_min] + (main_gpr_[i_next_min] - main_gpr_[i_min]) * diff;

                    if (!is_pvto || gor_[i] < 0)
                      gor_[i] = 0;
                  }
                else
                  {
                    diff = (max_p - main_pressure_[i_max]) / (main_pressure_[i_max] - main_pressure_[i_pre_max]);
                    inv_fvf_[i] = 1.0 / (main_fvf_[i_max] + (main_fvf_[i_max] - main_fvf_[i_pre_max]) * diff);
                    inv_visc_[i] = 1.0 / (main_visc_[i_max] + (main_visc_[i_max] - main_visc_[i_pre_max]) * diff);

                    if (inv_fvf_[i] < 0)
                      inv_fvf_[i] = inv_fvf_[i - 1];
                    if (inv_visc_[i] < 0)
                      inv_visc_[i] = inv_visc_[i - 1];

                    if (is_pvto)
                      gor_[i] = main_gpr_[i_max] + (main_gpr_[i_max] - main_gpr_[i_pre_max]) * diff;

                    if (!is_pvto || gor_[i] < 0)
                      gor_[i] = 0;
                  }

                inv_visc_fvf_[i] = inv_fvf_[i] * inv_visc_[i];
              }
            else
              {
                if (is_pvto)
                  {
                    while (main_gpr_[i_rs] < -0.5)
                      ++i_rs;
                  }

                pressure_[i]      = main_pressure_[i_rs];
                gor_[i]           = main_gpr_[i_rs];
                inv_fvf_[i]       = 1.0 / main_fvf_[i_rs];
                inv_visc_[i]      = 1.0 / main_visc_[i_rs];
                inv_visc_fvf_[i]  = inv_fvf_[i] * inv_visc_[i];
                ++i_rs;
              }
          }
      }

    return n_points;
  }

  void
  pvt_dead_oil::build (t_double atm_p, t_double min_p, t_double max_p, t_long n_intervals)
  {
    build_internal (atm_p, min_p, max_p, n_intervals, false);
  }

  void
  pvt_dead_oil::check_oil ()
  {
    int n_cols = pvt_input_props->get_n_cols ();
    int shift = 0;
    if (n_cols == 3)
      shift = 1;
    vector_t &main_pressure_     = pvt_input_props->get_col_vector (PVT_OIL_INPUT_PRESSURE - shift);
    vector_t &main_fvf_          = pvt_input_props->get_col_vector (PVT_OIL_INPUT_FVF - shift);
    vector_t &main_visc_         = pvt_input_props->get_col_vector (PVT_OIL_INPUT_VISC - shift);

    base_t::check_common ();

    check_oil_common (main_pressure_, main_fvf_, main_visc_);
    check_gas_common (main_pressure_, main_fvf_, main_visc_);
  }

  bool
  pvt_dead_oil::calc (const bool is_g, const int main_var, const t_double p, const t_double gor,
                                  t_float *inv_fvf, t_float *d_inv_fvf, t_float *gor_d_inv_fvf,
                                  t_float *inv_visc, t_float *d_inv_visc, t_float *gor_d_inv_visc,
                                  t_float *inv_visc_fvf, t_float *d_inv_visc_fvf, t_float *gor_d_inv_visc_fvf,
                                  t_float *gas_oil_ratio, t_float *d_gas_oil_ratio,
                                  const t_double drsdt /* = -1.0 */, const t_double dt /* = 0 */,
                                  const t_double old_gas_oil_ratio /* = 0  */) const
    {
      return calc_saturated_oil (is_g, main_var, p, gor, inv_fvf, d_inv_fvf, gor_d_inv_fvf,
                                 inv_visc, d_inv_visc, gor_d_inv_visc,
                                 inv_visc_fvf, d_inv_visc_fvf, gor_d_inv_visc_fvf,
                                 gas_oil_ratio, d_gas_oil_ratio, drsdt, dt, old_gas_oil_ratio);
    }

  bool
  pvt_dead_oil::calc_saturated_oil (const bool /*is_g*/, const int /*main_var*/, const t_double p, const t_double /*gor*/,
      t_float *inv_fvf, t_float *d_inv_fvf, t_float *gor_d_inv_fvf,
      t_float *inv_visc, t_float *d_inv_visc, t_float *gor_d_inv_visc,
      t_float *inv_visc_fvf, t_float *d_inv_visc_fvf, t_float *gor_d_inv_visc_fvf,
      t_float *gas_oil_ratio, t_float *d_gas_oil_ratio,
      const t_double /*drsdt = -1.0 */, const t_double /*dt = 0 */,
      const t_double /*old_gas_oil_ratio = 0 */) const
    {
      vector_t &pressure_     = pvt_props_table->get_col_vector (PVT_OIL_PRESSURE);
      vector_t &inv_fvf_      = pvt_props_table->get_col_vector (PVT_OIL_INV_FVF);
      vector_t &inv_visc_     = pvt_props_table->get_col_vector (PVT_OIL_INV_VISC);
      vector_t &inv_visc_fvf_ = pvt_props_table->get_col_vector (PVT_OIL_INV_VISC_FVF);
      vector_t &gor_          = pvt_props_table->get_col_vector (PVT_OIL_GOR);

      size_t i = binary_search (p, pressure_, std::less <t_double> ());
      if (i == 0)
        ++i;

      if (i >= pressure_.size ())
        {
          BOSERR (section::pvt, level::error) << "pvt_dead_oil::calc_saturated_oil: index_out_of_range: " << i << " [" << p << "]" << bs_end;
          return false;
        }

      t_double dp           = pressure_[i] - pressure_[i - 1];
      t_double idp          = 1.0 / dp;
      t_double diff_p       = p - pressure_[i - 1];

      t_double d_ifvf       = (inv_fvf_[i]          - inv_fvf_[i - 1]) * idp;
      t_double ifvf         = (inv_fvf_[i - 1]      + d_ifvf * diff_p);

      t_double d_ivisc      = (inv_visc_[i]         - inv_visc_[i - 1]) * idp;
      t_double ivisc        = (inv_visc_[i - 1]     + d_ivisc * diff_p);

      t_double d_ivisc_fvf  = (inv_visc_fvf_[i]     - inv_visc_fvf_[i - 1]) * idp;
      t_double ivisc_fvf    = (inv_visc_fvf_[i - 1] + d_ivisc_fvf * diff_p);

      t_double d_gor        = (gor_[i] - gor_[i - 1]) * idp;
      t_double cur_gor      = (gor_[i - 1] + d_gor * diff_p);

      set_pvt_pointer (inv_fvf, ifvf);
      set_pvt_pointer (inv_visc, ivisc);
      set_pvt_pointer (inv_visc_fvf, ivisc_fvf);
      set_pvt_pointer (d_inv_fvf, d_ifvf);
      set_pvt_pointer (d_inv_visc, d_ivisc);
      set_pvt_pointer (d_inv_visc_fvf, d_ivisc_fvf);

      set_pvt_pointer (gas_oil_ratio, cur_gor);
      set_pvt_pointer (d_gas_oil_ratio, d_gor);

      set_pvt_pointer (gor_d_inv_fvf, 0);
      set_pvt_pointer (gor_d_inv_visc, 0);
      set_pvt_pointer (gor_d_inv_visc_fvf, 0);

      return true;
    }

  t_double
  pvt_dead_oil::interpolate_and_fix (t_double cell_pbub) const
    {
      vector_t &pressure_     = pvt_props_table->get_col_vector (PVT_OIL_PRESSURE);
      //vector_t &inv_fvf_      = pvt_props_table->get_col_vector (PVT_OIL_INV_FVF);
      //vector_t &inv_visc_     = pvt_props_table->get_col_vector (PVT_OIL_INV_VISC);
      //vector_t &inv_visc_fvf_ = pvt_props_table->get_col_vector (PVT_OIL_INV_VISC_FVF);
      vector_t &gor_          = pvt_props_table->get_col_vector (PVT_OIL_GOR);

      size_t l = binary_search (cell_pbub, pressure_, std::less <t_double> ());
      size_t n = pressure_.size ();
      if (n == 1)
        {
          if (l == n)
            return gor_.back ();
          else if (l == 0)
            return gor_.front ();
          else
            throw bs_exception ("pvt_oil::interpolate_and_fix", "invalid index");
        }
      else if (n > 0)
        {
          if (l == n)
            l--;
          else if (l == 0)
            l++;

          t_double g2 = gor_[l - 1], g1 = gor_[l];
          t_double p2 = pressure_[l - 1], p1 = pressure_[l];
          return g2 + (cell_pbub - p2) * (g1 - g2)	/ (p1 - p2);
        }
      else
        {
          bs_throw_exception (boost::format ("Unexpected error: n (%lu) out of range") % n);
        }
    }

  t_double
  pvt_dead_oil::get_gor_for_pressure (t_double pressure_data) const
    {
      vector_t &pressure_     = pvt_props_table->get_col_vector (PVT_OIL_PRESSURE);
      vector_t &gor_          = pvt_props_table->get_col_vector (PVT_OIL_GOR);

      if (pressure_data < pressure_.front ())
        {
          throw bs_exception ("pvt_oil::get_gor_for_pressure", "Invalid pressure data value");
        }

      t_long i_rs = (t_long) (((pressure_data) - pressure_.front ()) / this->get_p_step ());

      return gor_[i_rs] + (gor_[i_rs + 1] - gor_[i_rs]) / this->get_p_step () * (pressure_data - pressure_[i_rs]);
    }

  void
  pvt_dead_oil::print () const
  {
    vector_t &pressure_     = pvt_props_table->get_col_vector (PVT_OIL_PRESSURE);
    vector_t &inv_fvf_      = pvt_props_table->get_col_vector (PVT_OIL_INV_FVF);
    vector_t &inv_visc_     = pvt_props_table->get_col_vector (PVT_OIL_INV_VISC);
    vector_t &inv_visc_fvf_ = pvt_props_table->get_col_vector (PVT_OIL_INV_VISC_FVF);
    vector_t &gor_          = pvt_props_table->get_col_vector (PVT_OIL_GOR);

    BS_ASSERT (pressure_.size () == inv_fvf_.size ());
    BS_ASSERT (inv_fvf_.size ()  == inv_visc_.size ());
    BS_ASSERT (inv_visc_.size () == inv_visc_fvf_.size ());

    BOSOUT (section::pvt, level::medium) << bs_end;
    BOSOUT (section::pvt, level::medium) << "*************************************** PVDO ****************************************" << bs_end;
    BOSOUT (section::pvt, level::medium)
      << boost::format ("*%9s%13s%12s%12s%12s%12s%12s%12s%2s\n")
      % "P"
      % "RS"
      % "FVF"
      % "1/FVF"
      % "VISC"
      % "1/VISC"
      % "VISC*FVF"
      % "1/FVF/VISC"
      % "*"
      << bs_end;
    BOSOUT (section::pvt, level::medium) << "*************************************************************************************" << bs_end;
    for (size_t j = 0, jcnt = pressure_.size ();j < jcnt; ++j)
      {
        BOSOUT (section::pvt, level::medium)
          << boost::format (
            "*%9.6f\t\t%9.6f\t\t%9.6f\t\t%9.6f\t\t%9.6f\t\t%9.6f\t\t%9.6f\t\t%9.6f%2s"
          )
          % (pressure_[j])
          % (gor_[j])
          % (1.0 / inv_fvf_[j])
          % (inv_fvf_[j])
          % (1.0 / inv_visc_[j])
          % (inv_visc_[j])
          % (1.0 / (inv_visc_fvf_[j]))
          % (inv_visc_fvf_[j])
          % "*"
          << bs_end;
      }
    BOSOUT (section::pvt, level::medium) << "*************************************************************************************" << bs_end;
  }
  BLUE_SKY_TYPE_STD_CREATE (pvt_dead_oil);
  BLUE_SKY_TYPE_STD_COPY (pvt_dead_oil);

  BLUE_SKY_TYPE_IMPL(pvt_dead_oil,  pvt_base, "pvt_dead_oil", "Dead Oil PVT calculation class", "Dead Oil PVT calculation");

} // namespace blue_sky

