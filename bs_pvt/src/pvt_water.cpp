/**
 * \file pvt_water.cpp
 * \brief
 * \author Miryanov Sergey
 * \date 07.05.2008
 */
#include "bs_pvt_stdafx.h"

#include "pvt_water.h"

using namespace blue_sky::pvt_detail;

namespace blue_sky
  {

  pvt_water::pvt_water (bs_type_ctor_param)
  {
  }

  pvt_water::pvt_water (const pvt_water &pvt)
  : bs_refcounter (pvt)
  {
    if (this != &pvt)
      {
        *this = pvt;
        BS_ASSERT (false && "NOT IMPL");
      }
  }

  void
  pvt_water::insert_vector (const v_double &vec)
  {
    const int elem_count = 4;
    BS_ASSERT (!(vec.size() % elem_count)) (vec.size ()) (elem_count);

    if (!this->init_dependent)
      {
        //main_pressure_.clear ();
        //main_compressibility_.clear ();
        //main_fvf_.clear ();
        //main_viscosibility_.clear ();
        //main_visc_.clear ();
        //main_gpr_.clear ();

        this->init_dependent = true;
      }

    t_int n_points = (t_int) vec.size () / elem_count;
    
    if (pvt_input_props->init (n_points, PVT_WATER_INPUT_TOTAL))
      {
        throw bs_exception ("pvt_water::insert_vector in table", "Error: initializing table of properties");
      }

    pvt_input_props->set_col_name (PVT_WATER_INPUT_PRESSURE, "pressure");
    pvt_input_props->set_col_name (PVT_WATER_INPUT_FVF, "fvf");   
    pvt_input_props->set_col_name (PVT_WATER_INPUT_COMPRESSIBILITY, "compressibility");   
    pvt_input_props->set_col_name (PVT_WATER_INPUT_VISC, "visc");
    pvt_input_props->set_col_name (PVT_WATER_INPUT_VISCOSIBILITY, "viscosibility");
    pvt_input_props->set_col_name (PVT_WATER_INPUT_GPR, "gor");

    vector_t &main_pressure_        = pvt_input_props->get_col_vector (PVT_WATER_INPUT_PRESSURE);
    vector_t &main_fvf_             = pvt_input_props->get_col_vector (PVT_WATER_INPUT_FVF);
    vector_t &main_compressibility_ = pvt_input_props->get_col_vector (PVT_WATER_INPUT_COMPRESSIBILITY);
    vector_t &main_visc_            = pvt_input_props->get_col_vector (PVT_WATER_INPUT_VISC);
    vector_t &main_viscosibility_   = pvt_input_props->get_col_vector (PVT_WATER_INPUT_VISCOSIBILITY);
    vector_t &main_gpr_             = pvt_input_props->get_col_vector (PVT_WATER_INPUT_GPR);
    
    for (t_int i = 0; i < n_points; ++i)
      {
        main_pressure_[i]        = (vec[i * elem_count + 0]);
        main_fvf_[i]             = (vec[i * elem_count + 1]);
        main_compressibility_[i] = (vec[i * elem_count + 2]);
        main_visc_[i]            = (vec[i * elem_count + 3]);
        main_viscosibility_[i]   = (0);
        main_gpr_[i]             = (vec[i * elem_count + 2]);
      }
  }

  void
  pvt_water::build (t_double /*atm_p*/, t_double min_p, t_double max_p, t_long n_intervals)
  {
    vector_t &main_pressure_        = pvt_input_props->get_col_vector (PVT_WATER_INPUT_PRESSURE);
    vector_t &main_fvf_             = pvt_input_props->get_col_vector (PVT_WATER_INPUT_FVF);
    vector_t &main_compressibility_ = pvt_input_props->get_col_vector (PVT_WATER_INPUT_COMPRESSIBILITY);
    vector_t &main_visc_            = pvt_input_props->get_col_vector (PVT_WATER_INPUT_VISC);
    vector_t &main_viscosibility_   = pvt_input_props->get_col_vector (PVT_WATER_INPUT_VISCOSIBILITY);
    //vector_t &main_gpr_             = pvt_input_props->get_col_vector (PVT_WATER_INPUT_GPR);

    BS_ASSERT (n_intervals > 0) (n_intervals);
    if (base_t::init_dependent)
      {
        if (pvt_props_table->init (n_intervals, PVT_WATER_TOTAL))
          {
            throw bs_exception ("pvt_water::init table", "Error: initializing table of properties");
          }
        
        pvt_props_table->set_col_name (PVT_WATER_PRESSURE, "pressure");
        pvt_props_table->set_col_name (PVT_WATER_INV_FVF, "inv_fvf");   
        pvt_props_table->set_col_name (PVT_WATER_INV_VISC, "inv_visc");
        pvt_props_table->set_col_name (PVT_WATER_INV_VISC_FVF, "inv_visc_fvf");

        base_t::init_dependent = false;
      }


    vector_t &pressure_     = pvt_props_table->get_col_vector (PVT_WATER_PRESSURE);
    vector_t &inv_fvf_      = pvt_props_table->get_col_vector (PVT_WATER_INV_FVF);
    vector_t &inv_visc_     = pvt_props_table->get_col_vector (PVT_WATER_INV_VISC);
    vector_t &inv_visc_fvf_ = pvt_props_table->get_col_vector (PVT_WATER_INV_VISC_FVF);

    check_water ();
    check_pressure_interval (min_p, max_p);
    this->check_interval_numbers (n_intervals);

    base_t::p_step = (max_p - min_p) / (t_double)n_intervals;

    t_double pressure = min_p;
    for (t_long i = 0; i <= n_intervals; ++i, pressure += this->p_step)
      {
        t_double diff       = pressure - main_pressure_.front ();

        inv_fvf_[i]       = exp (main_compressibility_.front () * diff) / main_fvf_.front ();
        inv_visc_[i]      = exp (main_viscosibility_.front () * diff) / main_visc_.front ();

        inv_visc_fvf_[i]  = inv_visc_[i] * inv_fvf_[i];
        pressure_[i]      = pressure;
      }
  }

  void
  pvt_water::check_water ()
  {
    vector_t &main_pressure_        = pvt_input_props->get_col_vector (PVT_WATER_INPUT_PRESSURE);
    vector_t &main_fvf_             = pvt_input_props->get_col_vector (PVT_WATER_INPUT_FVF);
    vector_t &main_compressibility_ = pvt_input_props->get_col_vector (PVT_WATER_INPUT_COMPRESSIBILITY);
    vector_t &main_visc_            = pvt_input_props->get_col_vector (PVT_WATER_INPUT_VISC);
    vector_t &main_viscosibility_   = pvt_input_props->get_col_vector (PVT_WATER_INPUT_VISCOSIBILITY);
    //vector_t &main_gpr_             = pvt_input_props->get_col_vector (PVT_WATER_INPUT_GPR);
    this->check_common ();

    if (main_compressibility_.front () < 0)
      {
        // TODO: LOG
        throw bs_exception ("pvt_water::check_water", "compressibility in PVTW should be greter or equal than 0");
      }
    if (main_visc_.front () < 0)
      {
        // TODO: LOG
        throw bs_exception ("pvt_water::check_water", "viscosity at reference pressure in PVTW should be greter or equal than 0");
      }
    if (main_viscosibility_.front () < 0)
      {
        // TODO: LOG
        throw bs_exception ("pvt_water::check_water", "viscosibility in PVTW should be greter or equal than 0");
      }
    if (main_pressure_.front () < 0)
      {
        // TODO: LOG
        throw bs_exception ("pvt_water::check_water", "reference pressure in PVTW should be greter than 0");
      }
    if (main_fvf_.front () < 0)
      {
        // TODO: LOG
        throw bs_exception ("pvt_water::check_water", "FVF at reference pressure in PVTW should be greter than 0");
      }
  }

  void
  pvt_water::calc (const t_double p,
                               t_double *inv_fvf, t_double *d_inv_fvf,
                               t_double *inv_visc, t_double *d_inv_visc,
                               t_double *inv_visc_fvf, t_double *d_inv_visc_fvf) const
    {
      vector_t &pressure_     = pvt_props_table->get_col_vector (PVT_WATER_PRESSURE);
      vector_t &inv_fvf_      = pvt_props_table->get_col_vector (PVT_WATER_INV_FVF);
      vector_t &inv_visc_     = pvt_props_table->get_col_vector (PVT_WATER_INV_VISC);
      vector_t &inv_visc_fvf_ = pvt_props_table->get_col_vector (PVT_WATER_INV_VISC_FVF);
    
      if (p < pressure_.front () || p >= pressure_.back ())
        {
          t_long i = 0;
          if (p >= pressure_.back ())
            i = (t_long)pressure_.size () - 1;

          set_pvt_pointer (inv_fvf, inv_fvf_[i]);
          set_pvt_pointer (inv_visc, inv_visc_[i]);
          //set_pvt_pointer (inv_visc_fvf, inv_fvf_[i] * inv_visc_[i]);
          set_pvt_pointer (inv_visc_fvf, inv_visc_fvf_[i]);
          set_pvt_pointer (d_inv_fvf, 0.0);
          set_pvt_pointer (d_inv_visc, 0.0);
          set_pvt_pointer (d_inv_visc_fvf, 0.0);
        }
      else
        {
          t_double ddp              = 1.0 / this->p_step;
          t_long i               = (t_long)((p - pressure_.front ()) * ddp);
          if ((size_t)i >= pressure_.size ())
            throw bs_exception ("pvt_water::calc", "i out of range");

          t_double diff             = (p - pressure_[i]);

          t_double d_inv_fvf_       = (inv_fvf_[i+1] - inv_fvf_[i]) * ddp;
          t_double ifvf_            = inv_fvf_[i] + d_inv_fvf_ * diff;

          t_double d_inv_visc_      = (inv_visc_[i+1] - inv_visc_[i]) * ddp;
          t_double ivisc_           = inv_visc_[i] + d_inv_visc_ * diff;

          t_double d_inv_visc_fvf_  = (inv_visc_fvf_[i+1] - inv_visc_fvf_[i]) * ddp;
          t_double ivisc_fvf_       = inv_visc_fvf_[i] + d_inv_visc_fvf_ * diff;

          set_pvt_pointer (d_inv_fvf, d_inv_fvf_);
          set_pvt_pointer (d_inv_visc, d_inv_visc_);
          set_pvt_pointer (inv_fvf, ifvf_);
          set_pvt_pointer (inv_visc, ivisc_);
          set_pvt_pointer (inv_visc_fvf, ivisc_fvf_);
          set_pvt_pointer (d_inv_visc_fvf, d_inv_visc_fvf_);
        }
    }

  void
  pvt_water::print () const
  {
    vector_t &pressure_     = pvt_props_table->get_col_vector (PVT_WATER_PRESSURE);
    vector_t &inv_fvf_      = pvt_props_table->get_col_vector (PVT_WATER_INV_FVF);
    vector_t &inv_visc_     = pvt_props_table->get_col_vector (PVT_WATER_INV_VISC);
    vector_t &inv_visc_fvf_ = pvt_props_table->get_col_vector (PVT_WATER_INV_VISC_FVF);
  
    BS_ASSERT (pressure_.size () == inv_fvf_.size ());
    BS_ASSERT (inv_fvf_.size ()  == inv_visc_.size ());
    BS_ASSERT (inv_visc_.size () == inv_visc_fvf_.size ());
    
    BOSOUT (section::pvt, level::medium) << bs_end;
    BOSOUT (section::pvt, level::medium) << "************************************** PVTW *****************************************" << bs_end;
    BOSOUT (section::pvt, level::medium) 
      << boost::format ("*%9s%13s%12s%12s%12s%12s%12s%2s")
      % "P"
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
            "*%9.6f\t\t%9.6f\t\t%9.6f\t\t%9.6f\t\t%9.6f\t\t%9.6f\t\t%9.6f%2s"
          )
          % (pressure_[j])
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
  BLUE_SKY_TYPE_STD_CREATE (pvt_water);
  BLUE_SKY_TYPE_STD_COPY (pvt_water);

  BLUE_SKY_TYPE_IMPL(pvt_water,  pvt_base, "pvt_water", "Water PVT calculation class", "Water PVT calculation");

} // namespace blue_sky

