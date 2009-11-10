/**
 * \file pvt_water.cpp
 * \brief
 * \author Miryanov Sergey
 * \date 07.05.2008
 */
#include "bs_pvt_stdafx.h"

#include "pvt_water.h"

using namespace blue_sky::pvt;

namespace blue_sky
  {

  template <typename strategy_t>
  pvt_water<strategy_t>::pvt_water (bs_type_ctor_param)
  {
  }

  template <typename strategy_t>
  pvt_water<strategy_t>::pvt_water (const pvt_water &pvt)
  : bs_refcounter (pvt)
  {
    if (this != &pvt)
      {
        *this = pvt;
        BS_ASSERT (false && "NOT IMPL");
      }
  }

  template <typename strategy_t>
  void
  pvt_water<strategy_t>::insert_vector (const input_vector_t &vec)
  {
    const int elem_count = 4;
    BS_ASSERT (!(vec.size() % elem_count)) (vec.size ()) (elem_count);

    if (!this->init_dependent)
      {
        main_pressure_.clear ();
        main_compressibility_.clear ();
        main_fvf_.clear ();
        main_viscosibility_.clear ();
        main_visc_.clear ();
        main_gpr_.clear ();

        this->init_dependent = true;
      }

    for (int i = 0, cnt = (int)(vec.size() / elem_count); i < cnt; ++i)
      {
        main_pressure_.push_back  (vec[i * elem_count + 0]);
        main_compressibility_.push_back (vec[i * elem_count + 2]);
        main_fvf_.push_back (vec[i * elem_count + 1]);
        //main_viscosibility_.push_back (vec[i * elem_count + 3]);
        main_viscosibility_.push_back (0);
        main_visc_.push_back (vec[i * elem_count + 3]);
        main_gpr_.push_back (vec[i * elem_count + 2]);
      }
  }

  template <typename strategy_t>
  void
  pvt_water<strategy_t>::build (item_t atm_p, item_t min_p, item_t max_p, index_t n_intervals)
  {
    BS_ASSERT (n_intervals > 0) (n_intervals);
    if (base_t::init_dependent)
      {
        pressure_.assign (n_intervals + 1, 0);
        inv_fvf_.assign (n_intervals + 1, 0);
        inv_visc_.assign (n_intervals + 1, 0);
        inv_visc_fvf_.assign (n_intervals + 1, 0);

        base_t::init_dependent = false;
      }

    check_water ();
    check_pressure_interval (min_p, max_p);
    this->check_interval_numbers (n_intervals);

    base_t::p_step = (max_p - min_p) / (item_t)n_intervals;

    item_t pressure = min_p;
    for (index_t i = 0; i <= n_intervals; ++i, pressure += this->p_step)
      {
        item_t diff       = pressure - main_pressure_.front ();

        inv_fvf_[i]       = exp (main_compressibility_.front () * diff) / main_fvf_.front ();
        inv_visc_[i]      = exp (main_viscosibility_.front () * diff) / main_visc_.front ();

        inv_visc_fvf_[i]  = inv_visc_[i] * inv_fvf_[i];
        pressure_[i]      = pressure;
      }
  }

  template <typename strategy_t>
  void
  pvt_water<strategy_t>::check_water ()
  {
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

  template <typename strategy_t>
  void
  pvt_water<strategy_t>::calc (const item_t p,
                               item_t *inv_fvf, item_t *d_inv_fvf,
                               item_t *inv_visc, item_t *d_inv_visc,
                               item_t *inv_visc_fvf, item_t *d_inv_visc_fvf) const
    {
      if (p < pressure_.front () || p >= pressure_.back ())
        {
          index_t i = 0;
          if (p >= pressure_.back ())
            i = (index_t)pressure_.size () - 1;

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
          item_t ddp              = 1.0 / this->p_step;
          index_t i               = (index_t)((p - pressure_.front ()) * ddp);
          if ((size_t)i >= pressure_.size ())
            throw bs_exception ("pvt_water::calc", "i out of range");

          item_t diff             = (p - pressure_[i]);

          item_t d_inv_fvf_       = (inv_fvf_[i+1] - inv_fvf_[i]) * ddp;
          item_t ifvf_            = inv_fvf_[i] + d_inv_fvf_ * diff;

          item_t d_inv_visc_      = (inv_visc_[i+1] - inv_visc_[i]) * ddp;
          item_t ivisc_           = inv_visc_[i] + d_inv_visc_ * diff;

          item_t d_inv_visc_fvf_  = (inv_visc_fvf_[i+1] - inv_visc_fvf_[i]) * ddp;
          item_t ivisc_fvf_       = inv_visc_fvf_[i] + d_inv_visc_fvf_ * diff;

          set_pvt_pointer (d_inv_fvf, d_inv_fvf_);
          set_pvt_pointer (d_inv_visc, d_inv_visc_);
          set_pvt_pointer (inv_fvf, ifvf_);
          set_pvt_pointer (inv_visc, ivisc_);
          set_pvt_pointer (inv_visc_fvf, ivisc_fvf_);
          set_pvt_pointer (d_inv_visc_fvf, d_inv_visc_fvf_);
        }
    }

  template <typename strategy_t>
  void
  pvt_water <strategy_t>::print () const
  {
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

  template class pvt_water <base_strategy_fi>;
  template class pvt_water <base_strategy_di>;
  template class pvt_water <base_strategy_mixi>;
} // namespace blue_sky

