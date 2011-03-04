/**
 * \file pvt_gas.cpp
 * \brief
 * \author Miryanov Sergey
 * \date 07.05.2008
 */
#include "bs_pvt_stdafx.h"

#include "pvt_gas.h"
#include "pvt_interpolator.h"
#include "scal_interpolate.h"

using namespace blue_sky::pvt;

namespace blue_sky
  {

  template <typename item_t, typename main_pressure_t>
  size_t
  get_interval (item_t pressure, const main_pressure_t &main_pressure)
  {
    return binary_search (pressure, main_pressure, std::less <item_t> ());
  }

  pvt_gas::pvt_gas (bs_type_ctor_param)
  {
  }

  pvt_gas::pvt_gas (const pvt_gas &pvt)
  : bs_refcounter (pvt)
  {
    if (this != &pvt)
      {
        *this = pvt;
        BS_ASSERT (false && "NOT IMPL");
      }
  }

  void
  pvt_gas::insert_vector (const input_vector_t &vec)
  {
    const int elem_count = 3;
    BS_ASSERT (!(vec.size() % elem_count)) (vec.size ()) (elem_count);

    if (!this->init_dependent)
      {
        main_pressure_.clear ();
        main_fvf_.clear ();
        main_visc_.clear ();

        this->init_dependent = true;
      }

    for (int i = 0, cnt = (int)(vec.size() / elem_count); i < cnt; ++i)
      {
        main_pressure_.push_back  (vec[i * elem_count + 0]);
        main_fvf_.push_back       (vec[i * elem_count + 1]);
        main_visc_.push_back      (vec[i * elem_count + 2]);
      }
  }

  void
  pvt_gas::build (item_t atm_p, item_t min_p, item_t max_p, t_int n_intervals)
  {
    if (this->init_dependent)
      {
        pressure_.assign (n_intervals + 1, 0);
        inv_fvf_.assign (n_intervals + 1, 0);
        inv_visc_.assign (n_intervals + 1, 0);
        inv_visc_fvf_.assign (n_intervals + 1, 0);

        this->init_dependent = false;
      }

    check_gas ();
    check_pressure_interval (min_p, max_p);
    this->check_interval_numbers (n_intervals);

    this->p_step = (max_p - min_p) / (item_t)n_intervals;

    if (main_pressure_.size () < 2)
      {
        item_t pressure = min_p;
        for (int i = 0; i <= n_intervals; ++i, pressure += this->p_step)
          {
            pressure_[i]      = pressure;
            inv_fvf_[i]       = 1.0 / main_fvf_.front ();
            inv_visc_[i]      = 1.0 / main_visc_.front ();
            inv_visc_fvf_[i]  = inv_fvf_[i] * inv_visc_[i];
          }
      }
    else
      {
        item_t pressure = min_p;
        for (int i = 0; i <= n_intervals; ++i, pressure += this->p_step)
          {
            size_t j = get_interval (pressure, main_pressure_);

            item_t diff             = 0;
            pressure_[i]            = pressure;

            if (j == 0)
              {
                inv_fvf_[i]         = 1.0 / main_fvf_.front ();
                inv_visc_[i]        = 1.0 / main_visc_.front ();
              }
            else if (j == main_pressure_.size ())
              {
                inv_fvf_[i]         = pvt::interpolate_x2 (pressure, main_pressure_[j-1], main_pressure_[j-2], main_fvf_[j-1], main_fvf_[j-2], diff);
                inv_visc_[i]        = 1.0 / pvt::interpolate_x3 (main_visc_[j-1], main_visc_[j-2], diff);
              }
            else
              {
                inv_fvf_[i]         = pvt::interpolate_x (pressure, main_pressure_[j], main_pressure_[j-1], main_fvf_[j], main_fvf_[j-1], diff);
                inv_visc_[i]        = pvt::interpolate_x (main_visc_[j], main_visc_[j-1], diff);
              }

            inv_visc_fvf_[i]        = inv_visc_[i] * inv_fvf_[i];
          }
      }
  }


  void
  pvt_gas::check_gas ()
  {
    this->check_common ();

    if (main_pressure_.empty ())
      {
        bs_throw_exception ("Main_pressure is empty");
      }

    check_oil_common (main_pressure_, main_fvf_, main_visc_);
    check_gas_common (main_pressure_, main_fvf_, main_visc_);
  }

  void
  pvt_gas::calc (const item_t p,
                             item_t *inv_fvf, item_t *d_inv_fvf,
                             item_t *inv_visc, item_t *d_inv_visc,
                             item_t *inv_visc_fvf, item_t *d_inv_visc_fvf) const
    {
      item_t idp = 1.0 / this->p_step;
      index_t i = (index_t)((p - pressure_.front ()) * idp);

      if ((size_t)(i + 1) >= pressure_.size () || p < pressure_.front ())
        {
          index_t l = 0;
          if (i + 1 >= (index_t)pressure_.size ())
            l = (index_t)pressure_.size () - 1;

          set_pvt_pointer (inv_fvf, inv_fvf_[l]);
          set_pvt_pointer (inv_visc, inv_visc_[l]);
          set_pvt_pointer (inv_visc_fvf, inv_visc_[l] * inv_fvf_[l]);
          set_pvt_pointer (d_inv_fvf, 0.0);
          set_pvt_pointer (d_inv_visc, 0.0);
          set_pvt_pointer (d_inv_visc_fvf, 0.0);
        }
      else
        {
          item_t diff         = (p - pressure_[i]);

          item_t d_inv_fvf_   = (inv_fvf_[i+1] - inv_fvf_[i]) * idp;
          item_t ifvf_        = inv_fvf_[i] + d_inv_fvf_ * diff;

          item_t d_inv_visc_  = (inv_visc_[i+1] - inv_visc_[i]) * idp;
          item_t ivisc_       = inv_visc_[i] + d_inv_visc_ * diff;

          set_pvt_pointer (inv_fvf, ifvf_);
          set_pvt_pointer (inv_visc, ivisc_);
          set_pvt_pointer (inv_visc_fvf, ivisc_ * ifvf_);

          set_pvt_pointer (d_inv_fvf, d_inv_fvf_);
          set_pvt_pointer (d_inv_visc, d_inv_visc_);
          set_pvt_pointer (d_inv_visc_fvf, ifvf_ * d_inv_visc_ + ivisc_ * d_inv_fvf_);
        }
    }

  void
  pvt_gas::print () const
  {
    BS_ASSERT (pressure_.size () == inv_fvf_.size ());
    BS_ASSERT (inv_fvf_.size ()  == inv_visc_.size ());
    BS_ASSERT (inv_visc_.size () == inv_visc_fvf_.size ());
    
    BOSOUT (section::pvt, level::medium) << bs_end;
    BOSOUT (section::pvt, level::medium) << "************************************** PVTG *****************************************" << bs_end;
    BOSOUT (section::pvt, level::medium) 
      << boost::format ("*%9s%13s%12s%12s%12s%12s%12s%2s\n")
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

  BLUE_SKY_TYPE_STD_CREATE (pvt_gas);
  BLUE_SKY_TYPE_STD_COPY (pvt_gas);

  BLUE_SKY_TYPE_IMPL(pvt_gas,  pvt_base, "pvt_gas", "Gas PVT calculation class", "Gas PVT calculation");

} // namespace blue_sky

