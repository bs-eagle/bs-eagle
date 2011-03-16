/**
 * \file pvt_base.cpp
 * \brief
 * \author Miryanov Sergey
 * \date 07.05.2008
 */
#include "bs_pvt_stdafx.h"

#include "pvt_base.h"

namespace blue_sky
  {

  pvt_base::pvt_base ( )
  {
    p_step = 0;
    surface_density = 0;
    init_dependent = true;
    pvt_input_props = BS_KERNEL.create_object ("table");
    pvt_props_table = BS_KERNEL.create_object ("table"); 
  }

  t_double
  pvt_base::interpolate_and_fix (t_double /*cell_pbub*/) const
    {
      BS_ASSERT (false && "BASE METHOD CALL");
      return 0;
    }

  t_double
  pvt_base::get_p_step () const
    {
      return p_step;
    }

  t_double
  pvt_base::get_surface_density () const
    {
      return surface_density;
    }

  void pvt_base::set_surface_density (t_double density)
  {
    surface_density = density;
  }

  void
  pvt_base::set_density(t_double density, t_double md)
  {
    surface_density = density;
    molar_density = md;
  }

  t_double
  pvt_base::get_gor_for_pressure (t_double /*pressure_data*/) const
    {
      BS_ASSERT (false && "BASE METHOD CALL");
      return 0;
    }

  void pvt_base::check_pressure_interval (t_double min_p, t_double max_p)
  {
    if (max_p - min_p < 1.0e-12)
      {
        // TODO: LOG
        BS_ASSERT (false) (max_p) (min_p);
        throw bs_exception ("", "invalid pressure interval");
      }
  }

  void pvt_base::check_interval_numbers (t_long &n_intervals)
  {
    if (n_intervals < 1)
      {
        // TODO: LOG
        n_intervals = 1;
      }
  }

  void pvt_base::check_common ()
  {
    if (surface_density < 0)
      {
        // TODO: LOG
        BS_ASSERT (false) (surface_density);
        throw bs_exception ("", "surface density lower than 0");
      }
  }

  void pvt_base::check_gas_common (const vector_t &pressure, const vector_t &fvf, const vector_t &visc)
  {
    for (t_long i = 1, cnt = (t_long)pressure.size (); i < cnt; ++i)
      {
        if (pressure[i] - pressure[i - 1] < EPS_DIFF)
          {
            throw bs_exception ("", "pressure curve should be monotonically increasing function");
          }
        if (fvf[i] - fvf[i - 1] >= 0)
          {
            BOSOUT (section::pvt, level::critical) << "gas: fvf" << bs_end;
            for (t_long j = 0; j < cnt; ++j)
              {
                BOSOUT (section::pvt, level::critical) << fvf[j] << bs_end;
              }

            throw bs_exception ("", "FVF curve should be monotonically decreasing function");
          }
        if (visc[i] - visc[i - 1] <= 0)
          {
            BOSOUT (section::pvt, level::critical) << "gas: visc" << bs_end;
            for (t_long j = 0; j < cnt; ++j)
              {
                BOSOUT (section::pvt, level::critical) << visc[j] << bs_end;
              }

            throw bs_exception ("", "Viscosity curve should be monotonically increasing function");
          }
      }
  }

  void 
  pvt_base::check_oil_common (const vector_t &pressure, const vector_t &fvf, const vector_t &visc)
  {
    for (t_long i = 0, cnt = (t_long)pressure.size (); i < cnt; ++i)
      {
        if (pressure[i] < 0)
          {
            // TODO: LOG
            BS_ASSERT (false) (pressure[i]);
            throw bs_exception ("", "pressure should be greater than 0");
          }
        if (fvf[i] < 0)
          {
            // TODO:LOG
            BS_ASSERT (false) (fvf[i]);
            throw bs_exception ("", "fvf should be greater than 0");
          }
        if (visc[i] < 0)
          {
            // TODO: LOG
            BS_ASSERT (false) (visc[i]);
            throw bs_exception ("", "viscosity should be greater than 0");
          }
      }
  }

} // namespace blue_sky
