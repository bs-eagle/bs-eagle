/**
 *       \file  well_controller.cpp
 *      \brief  Implementation of well_controller
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  16.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"
#include "well_controller.h"
#include "calc_well.h"
#include "calc_model.h"
#include "well_connection.h"

#include "well_rate_control.h"

#include "reservoir.h"
#include "facility_manager.h"
#include "default_well_rate_control_factory.h"
#define WELL_ACCUMULATION_TERM 1e-6

namespace blue_sky
  {
  namespace wells
    {

    ///////////////////////////////////////////////////////////////////////////
    bool
    is_oil_rate_value (rate_value_type type)
    {
      return (type & oil_rate_value) == oil_rate_value;
    }
    bool
    is_water_rate_value (rate_value_type type)
    {
      return (type & water_rate_value) == water_rate_value;
    }
    bool
    is_gas_rate_value (rate_value_type type)
    {
      return (type & gas_rate_value) == gas_rate_value;
    }
    bool
    is_liquid_rate_value (rate_value_type type)
    {
      return (type & liquid_rate_value) == liquid_rate_value;
    }
    rate_value_type
    rate_value_cast (const std::string &str)
    {
      if (str == "ORAT")
        return oil_rate_value;
      else if (str == "WRAT")
        return water_rate_value;
      else if (str == "GRAT")
        return gas_rate_value;
      else if (str == "BHP")
        return bhp_value;
      else if (str == "LRAT")
        return liquid_rate_value;
      else if (str == "RESV")
        {
          BS_ASSERT (false && "Unsupported rate_value type") (str);
          throw bs_exception ("rate_value_cast", "Unsupported rate_value type");
        }
      else if (str == "THP")
        {
          BS_ASSERT (false && "Unsupported rate_value type") (str);
          throw bs_exception ("rate_value_cast", "Unsupported rate_value type");
        }
      else
        {
          BS_ASSERT (false && "Unsupported rate_value type") (str);
          throw bs_exception ("rate_value_cast", "Unsupported rate_value type");
        }
    }
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    injection_type
    injection_type_cast (const std::string &str)
    {
      if (str == "WATER")
        return injection_water;
      else if (str == "GAS")
        return injection_gas;
      else if (str == "OIL")
        return injection_oil;
      else if (str == "NONE")
        return injection_none;
      else
        {
          BS_ASSERT (false && "Unsupported injection type") (str);
          throw bs_exception ("injection_type_cast", "Unsupported injection type");
        }
    }
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    template <typename strategy_t>
    typename well_controller<strategy_t>::sp_rate_control_t well_controller<strategy_t>::dummy_control_ = 0;

    template <typename strategy_t>
    well_controller<strategy_t>::well_controller (bs_type_ctor_param param /* = NULL */)
    : bhp_ (0)
    , bhp_history_ (0)
    , injection_type_ (injection_none)
    {
      current_control_ = 0;
    }

    template <typename strategy_t>
    well_controller<strategy_t>::well_controller (const well_controller<strategy_t> &controller)
    : bs_refcounter (controller), objbase (controller)
    {
      if (this != &controller)
        {
          *this = controller;
          //this_ = this;
        }
    }

    template <typename strategy_t>
    const typename well_controller<strategy_t>::rate_data_t &
    well_controller<strategy_t>::rate () const
      {
        return rate_;
      }
    template <typename strategy_t>
    const typename well_controller<strategy_t>::item_t &
    well_controller<strategy_t>::bhp () const
      {
        return bhp_;
      }
    template <typename strategy_t>
    const typename well_controller<strategy_t>::item_t &
    well_controller<strategy_t>::bhp_history () const
      {
        return bhp_history_;
      }
    ///////////////////////////////////////////////////////////////////////////
    template <typename strategy_t>
    void
    well_controller<strategy_t>::add_bhp_control  (const sp_rate_control_t &control)
    {
      BS_ASSERT (control);
      bhp_control_ = control;
    }
    template <typename strategy_t>
    void
    well_controller<strategy_t>::add_rate_control (const sp_rate_control_t &control)
    {
      BS_ASSERT (control);
      rate_control_ = control;
    }

    template <typename strategy_t>
    void
    well_controller<strategy_t>::set_main_control (const sp_well_t &well, rate_control_type type)
    {
      current_control_ = 0;
      BS_ASSERT (bhp_control_);
      BS_ASSERT (rate_control_);

      if (is_bhp_control (type))
        {
          BOSOUT (section::wells, level::low) << "[" << well->name () << "] set control to bhp" << bs_end;

          current_control_ = bhp_control_;
          well->set_bhp (bhp_);
        }
      else if (is_rate_control (type))
        {
          BOSOUT (section::wells, level::low) << "[" << well->name () << "] set control to rate" << bs_end;

          current_control_ = rate_control_;
        }

      if (!current_control_)
        {
          bs_throw_exception ("Current control is null");
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    template <typename strategy_t>
    void
    well_controller <strategy_t>::set_bhp (item_t value)
    {
      bhp_ = value;
    }
    template <typename strategy_t>
    void
    well_controller <strategy_t>::set_bhp_history (item_t value)
    {
      bhp_history_ = value;
    }
    template <typename strategy_t>
    void
    well_controller <strategy_t>::clear_rate ()
    {
      rate_ = 0;
    }
    template <typename strategy_t>
    void
    well_controller<strategy_t>::set_rate (rate_value_type type, item_t value)
    {
      if (!is_production ())
        {
          BS_ASSERT (injection_type_ != injection_none);

          if (type == rate_value)
            {
              if (injection_type_ == injection_oil)
                {
                  rate_.inj.oil = value;
                  rate_.inj.liquid = rate_.inj.oil + rate_.inj.water;
                }
              else if (injection_type_ == injection_water)
                {
                  rate_.inj.water = value;
                  rate_.inj.liquid = rate_.inj.oil + rate_.inj.water;
                }
              else if (injection_type_ == injection_gas)
                {
                  rate_.inj.gas = value;
                }
            }
          else if (type == liquid_inner_rate_value)
            {
              rate_.liquid_inner = value;
            }
          else
            {
              throw bs_exception ("well_controller::set_rate", "invalid control value type");
            }
        }
      else
        {
          switch (type)
            {
            case oil_rate_value:
              rate_.prod.oil = value;
              rate_.prod.liquid = rate_.prod.oil + rate_.prod.water;
              break;
            case water_rate_value:
              rate_.prod.water = value;
              rate_.prod.liquid = rate_.prod.oil + rate_.prod.water;
              break;
            case gas_rate_value:
              rate_.prod.gas = value;
              break;
            case liquid_rate_value:
              rate_.prod.liquid = value;
              break;
            case liquid_inner_rate_value:
              rate_.liquid_inner = value;
              break;
            default:
              throw bs_exception ("well_controller::set_rate", "invalid control value type");
            }
        }
    }

    template <typename strategy_t>
    void
    well_controller<strategy_t>::save_control ()
    {
      BS_ASSERT (current_control_);
      saved_control_ = current_control_;
    }
    template <typename strategy_t>
    bool
    well_controller<strategy_t>::restore_control ()
    {
      BS_ASSERT (current_control_);
      bool flag = current_control_ != saved_control_;
      current_control_ = saved_control_;

      return flag;
    }
    template <typename strategy_t>
    void
    well_controller<strategy_t>::save_niter_control ()
    {
      BS_ASSERT (current_control_);
      saved_niter_control_ = current_control_;
    }
    template <typename strategy_t>
    bool
    well_controller<strategy_t>::restore_niter_control ()
    {
      BS_ASSERT (current_control_);
      bool flag = current_control_ != saved_niter_control_;
      current_control_ = saved_niter_control_;

      return flag;
    }

    ///////////////////////////////////////////////////////////////////////////
    template <typename strategy_t>
    bool
    well_controller<strategy_t>::check (sp_well_t &well)
    {
      BS_ASSERT (current_control_);
      BS_ASSERT (rate_control_);
      BS_ASSERT (bhp_control_);

      //BOSOUT (section::wells, level::low) << "[" << well->name () << "] " << "limit_bhp: "     << bhp_          << "\twell_bhp: "     << well->bhp ()           << bs_end;
      //BOSOUT (section::wells, level::low) << "[" << well->name () << "] " << "limit_liquid: "  << rate_.liquid   << "\twell_liquid: "  << -well->liquid_rate ()  << bs_end;
      //BOSOUT (section::wells, level::low) << "[" << well->name () << "] " << "limit_oil: "     << rate_.oil      << "\twell_oil: "     << -well->oil_rate ()     << bs_end;
      //BOSOUT (section::wells, level::low) << "[" << well->name () << "] " << "limit_water: "   << rate_.water    << "\twell_water: "   << -well->water_rate ()   << bs_end;
      //BOSOUT (section::wells, level::low) << "[" << well->name () << "] " << "limit_gas: "     << rate_.gas      << "\twell_gas: "     << -well->gas_rate ()     << bs_end;

      if (current_control_ == rate_control_)
        {
          bool do_switch = false;
          if (is_production ())
            {
              do_switch = bhp_ > well->bhp ();
            }
          else
            {
              do_switch = bhp_ < well->bhp ();
            }

          if (do_switch)
            {
              BOSOUT (section::wells, level::medium) << "[" << well->name () << "] switch from rate to bhp (is_prod: " << is_production() << ") " << bs_end;

              current_control_ = bhp_control_;
              well->set_bhp (bhp_);
              return true;
            }
        }
      else if (current_control_ == bhp_control_)
        {
          rate_control_type control_type = rate_control_->get_control_type ();
          bool do_switch = false;

          if (is_production ())
            {
              switch (control_type)
              {
                case liquid_rate_control:
                  do_switch = rate_.prod.liquid < -well->rate ().prod.liquid;
                  break;
                case oil_rate_control:
                  do_switch = rate_.prod.oil < -well->rate ().prod.oil;
                  break;
                case water_rate_control:
                  do_switch = rate_.prod.water < -well->rate ().prod.water;
                  break;
                case gas_rate_control:
                  do_switch = rate_.prod.gas < -well->rate ().prod.gas;
                  break;
                default:
                  bs_throw_exception (boost::format ("control_type is unknown: %d") % control_type);
              }
            }
          else
            {
              switch (injection_type_)
              {
                case injection_oil:
                  do_switch = rate_.inj.oil < well->rate ().inj.oil;
                  break;
                case injection_water:
                  do_switch = rate_.inj.water < well->rate ().inj.water;
                  break;
                case injection_gas:
                  do_switch = rate_.inj.gas < well->rate ().inj.gas;
                  break;
                default:
                  throw bs_exception ("well_controller::check", "injection_type is unknown");
              }
            }

          if (do_switch)
            {
              BOSOUT (section::wells, level::medium) << "[" << well->name () << "] switch from bhp to rate (is_prod: " << is_production() << ") " << bs_end;

              current_control_ = rate_control_;
              //if (is_production ())
              //  {
              //    well->bhp_ = bhp_ + 1;
              //  }
              //else
              //  {
              //    well->bhp_ = bhp_ - 1;
              //  }
              return true;
            }
          else
            {
              //return false; // in this case well should shutdown himself
            }
        }

      return false;
    }

    //////////////////////////////////////////////////////////////////////////
    template <typename strategy_t>
    bool
    well_controller<strategy_t>::is_bhp () const
      {
        BS_ASSERT (current_control_);
        return current_control_->is_bhp ();
      }
    template <typename strategy_t>
    bool
    well_controller<strategy_t>::is_rate () const
      {
        BS_ASSERT (current_control_);
        return current_control_->is_rate ();
      }
    template <typename strategy_t>
    bool
    well_controller<strategy_t>::is_production () const
      {
        BS_ASSERT (current_control_);
        return current_control_->is_production ();
      }

    template <typename strategy_t>
    bool
    well_controller <strategy_t>::is_valid_connection_bhp (item_t pressure, item_t bhp) const
      {
        item_t diff = pressure - bhp;
        if (is_production ())
          return diff > 0;
        else
          return diff < 0;
      }


    template <typename strategy_t>
    void
    well_controller<strategy_t>::switch_to_bhp (sp_well_t &well)
    {
      BS_ASSERT (current_control_);
      BS_ASSERT (bhp_control_);

      BOSOUT (section::wells, level::medium) << "[" << well->name () << "] switch to bhp" << bs_end;
      BOSOUT (section::wells, level::low)    << "\tlimit_bhp: " << bhp_ << "\twell_bhp: " << well->bhp () << bs_end;

      well->set_bhp (bhp_);
      current_control_ = bhp_control_;
    }

    template <typename strategy_t>
    const injection_type &
    well_controller<strategy_t>::injection () const
    {
      return injection_type_;
    }
    template <typename strategy_t>
    void
    well_controller<strategy_t>::set_injection_type (injection_type type)
    {
      injection_type_ = type;
    }

    template <typename strategy_t>
    void
    well_controller <strategy_t>::calc_rate (const sp_calc_model_t &calc_model, sp_well_t &well, sp_jmatrix_t &jmatrix) const
    {
      BS_ASSERT (current_control_);
      const smart_ptr <const this_t> sp_this (this);
      current_control_->compute_rate (calc_model, jmatrix, well, sp_this);
    }
    template <typename strategy_t>
    void
    well_controller <strategy_t>::calc_derivs (const sp_calc_model_t &calc_model, sp_well_t &well, sp_jmatrix_t &jmatrix) const
    {
      BS_ASSERT (current_control_);
      const smart_ptr <const this_t> sp_this (this);
      current_control_->compute_derivs (calc_model, jmatrix, well, sp_this);
    }

    template <typename strategy_t>
    rate_control_type
    well_controller <strategy_t>::get_control_type () const
      {
        BS_ASSERT (current_control_);
        return current_control_->get_control_type ();
      }

    //////////////////////////////////////////////////////////////////////////

    template <typename strategy_t>
    well_controller_factory<strategy_t>::well_controller_factory (bs_type_ctor_param param /* = NULL */)
    : well_rate_control_factory_ (BS_KERNEL.create_object (default_well_rate_control_factory <strategy_t>::bs_type ()))
    {

    }
    template <typename strategy_t>
    well_controller_factory<strategy_t>::well_controller_factory (const well_controller_factory &f)
    : bs_refcounter (f), objbase (f)
    {
      *this = f;
    }

    template <typename strategy_t>
    typename well_controller_factory<strategy_t>::sp_well_controller_t
    well_controller_factory<strategy_t>::create_controller () const
    {
      return BS_KERNEL.create_object (well_controller_t::bs_type (), true);
    }

    template <typename strategy_t>
    typename well_controller_factory<strategy_t>::sp_rate_control_t
    well_controller_factory<strategy_t>::create_control (rate_control_type control_type, bool is_prod, const sp_calc_model_t &calc_model) const
    {
      sp_rate_control_t control = BS_KERNEL.create_object (well_rate_control_t::bs_type (), true);
      if (!control)
        {
          bs_throw_exception ("Can't create well_rate_control");
        }

      control->set_is_bhp (is_bhp_control (control_type));
      control->set_is_prod (is_prod);
      control->set_control_type (control_type);
      control->set_impl (well_rate_control_factory_->create_control (control_type, is_bhp_control (control_type), is_prod, calc_model));

      return control;
    }

    template <typename strategy_t>
    void
    well_controller_factory <strategy_t>::set_rate_control_factory (const sp_well_rate_control_factory_t &f)
    {
      well_rate_control_factory_ = f;
    }

    //////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////
    BLUE_SKY_TYPE_STD_CREATE_T_DEF (well_controller, (class));
    BLUE_SKY_TYPE_STD_COPY_T_DEF (well_controller, (class));
    BLUE_SKY_TYPE_IMPL_T_EXT (1, (well_controller<base_strategy_fi>), 1, (objbase), "well_controller_seq_fi", "Base class for well controllers", "Base class for well controllers", false);
    BLUE_SKY_TYPE_IMPL_T_EXT (1, (well_controller<base_strategy_di>), 1, (objbase), "well_controller_seq_di", "Base class for well controllers", "Base class for well controllers", false);
    BLUE_SKY_TYPE_IMPL_T_EXT (1, (well_controller<base_strategy_mixi>), 1, (objbase), "well_controller_seq_mixi", "Base class for well controllers", "Base class for well controllers", false);

    BLUE_SKY_TYPE_STD_CREATE_T_DEF (well_controller_factory, (class));
    BLUE_SKY_TYPE_STD_COPY_T_DEF (well_controller_factory, (class));
    BLUE_SKY_TYPE_IMPL_T_EXT (1, (well_controller_factory<base_strategy_fi>), 1, (objbase), "well_controller_factory_seq_fi", "Base class for factory of well controllers", "Base class for factory of well controllers", false);
    BLUE_SKY_TYPE_IMPL_T_EXT (1, (well_controller_factory<base_strategy_di>), 1, (objbase), "well_controller_factory_seq_di", "Base class for factory of well controllers", "Base class for factory of well controllers", false);
    BLUE_SKY_TYPE_IMPL_T_EXT (1, (well_controller_factory<base_strategy_mixi>), 1, (objbase), "well_controller_factory_seq_mixi", "Base class for factory of well controllers", "Base class for factory of well controllers", false);
    //////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////
    bool
    well_controller_register_type (const blue_sky::plugin_descriptor &pd)
    {
      bool res = true;

      res &= BS_KERNEL.register_type (pd, well_controller<base_strategy_fi>::bs_type ());
      BS_ASSERT (res);
      res &= BS_KERNEL.register_type (pd, well_controller<base_strategy_di>::bs_type ());
      BS_ASSERT (res);
      res &= BS_KERNEL.register_type (pd, well_controller<base_strategy_mixi>::bs_type ());
      BS_ASSERT (res);

      return res;
    }
    bool
    well_controller_factory_register_type (const blue_sky::plugin_descriptor &pd)
    {
      bool res = true;

      res &= BS_KERNEL.register_type (pd, well_controller_factory<base_strategy_fi>::bs_type ());
      BS_ASSERT (res);
      res &= BS_KERNEL.register_type (pd, well_controller_factory<base_strategy_di>::bs_type ());
      BS_ASSERT (res);
      res &= BS_KERNEL.register_type (pd, well_controller_factory<base_strategy_mixi>::bs_type ());
      BS_ASSERT (res);

      return res;
    }
    //////////////////////////////////////////////////////////////////////////

  } // namespace wells
} // namespace blue_sky

