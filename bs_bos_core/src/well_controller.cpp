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

#include "reservoir.h"
#include "facility_manager.h"
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

    well_controller::well_controller (bs_type_ctor_param /*param = NULL */)
    : bhp_ (0)
    , bhp_history_ (0)
    , injection_type_ (injection_none)
    {
    }

    well_controller::well_controller (const well_controller &controller)
    : bs_refcounter (controller), objbase (controller)
    {
      if (this != &controller)
        {
          *this = controller;
          //this_ = this;
        }
    }

    const well_controller::rate_data_t &
    well_controller::rate () const
      {
        return rate_;
      }
    const well_controller::item_t &
    well_controller::bhp () const
      {
        return bhp_;
      }
    const well_controller::item_t &
    well_controller::bhp_history () const
      {
        return bhp_history_;
      }
    ///////////////////////////////////////////////////////////////////////////
    void
    well_controller::set_main_control (const sp_well_t &well, rate_control_type type, bool is_production)
    {
      if (is_bhp_control (type))
        {
          BOSOUT (section::wells, level::low) << "[" << well->name () << "] set control to bhp" << bs_end;

          current_control_.is_bhp         = true;
          current_control_.is_production  = is_production;
          current_control_.control_type   = type;
          well->set_bhp (bhp_);
        }
      else if (is_rate_control (type))
        {
          BOSOUT (section::wells, level::low) << "[" << well->name () << "] set control to rate" << bs_end;

          current_control_.is_bhp         = false;
          current_control_.is_production  = is_production;
          current_control_.control_type   = type;
        }
      else
        {
          bs_throw_exception ("Unsupported rate_control_type");
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    void
    well_controller::set_bhp (item_t value)
    {
      bhp_ = value;
    }
    void
    well_controller::set_bhp_history (item_t value)
    {
      bhp_history_ = value;
    }
    void
    well_controller::clear_rate ()
    {
      rate_ = 0;
    }
    void
    well_controller::set_rate (rate_value_type type, item_t value)
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
              bs_throw_exception (boost::format ("Invalid control value type (%d)") % type);
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
              bs_throw_exception (boost::format ("Invalid control value type (%d)") % type);
            }
        }
    }

    void
    well_controller::save_control ()
    {
      saved_control_ = current_control_;
    }
    bool
    well_controller::restore_control ()
    {
      bool flag = current_control_ != saved_control_;
      current_control_ = saved_control_;

      return flag;
    }
    void
    well_controller::save_niter_control ()
    {
      saved_niter_control_ = current_control_;
    }
    bool
    well_controller::restore_niter_control ()
    {
      bool flag = current_control_ != saved_niter_control_;
      current_control_ = saved_niter_control_;

      return flag;
    }

    ///////////////////////////////////////////////////////////////////////////
    bool
    well_controller::check (sp_well_t &well)
    {
      //BOSOUT (section::wells, level::low) << "[" << well->name () << "] " << "limit_bhp: "     << bhp_          << "\twell_bhp: "     << well->bhp ()           << bs_end;
      //BOSOUT (section::wells, level::low) << "[" << well->name () << "] " << "limit_liquid: "  << rate_.liquid   << "\twell_liquid: "  << -well->liquid_rate ()  << bs_end;
      //BOSOUT (section::wells, level::low) << "[" << well->name () << "] " << "limit_oil: "     << rate_.oil      << "\twell_oil: "     << -well->oil_rate ()     << bs_end;
      //BOSOUT (section::wells, level::low) << "[" << well->name () << "] " << "limit_water: "   << rate_.water    << "\twell_water: "   << -well->water_rate ()   << bs_end;
      //BOSOUT (section::wells, level::low) << "[" << well->name () << "] " << "limit_gas: "     << rate_.gas      << "\twell_gas: "     << -well->gas_rate ()     << bs_end;

      if (!current_control_.is_bhp)
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

              current_control_.is_bhp = true;
              well->set_bhp (bhp_);
              return true;
            }
        }
      else if (current_control_.is_bhp)
        {
          bool do_switch = false;
          if (current_control_.is_production)
            {
              switch (current_control_.control_type)
              {
                case bhp_control:
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
                  bs_throw_exception (boost::format ("Control_type is unknown: %d (well: %s)") % current_control_.control_type % well->name ());
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
                  bs_throw_exception (boost::format ("Injection_type is unknown: %d (well: %s)") % injection_type_ % well->name ());
              }
            }

          if (do_switch)
            {
              BOSOUT (section::wells, level::medium) << "[" << well->name () << "] switch from bhp to rate (is_prod: " << is_production() << ") " << bs_end;

              current_control_.is_bhp = false;
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
    bool
    well_controller::is_bhp () const
      {
        return current_control_.is_bhp;
      }
    bool
    well_controller::is_rate () const
      {
        return !current_control_.is_bhp;
      }
    bool
    well_controller::is_production () const
      {
        return current_control_.is_production;
      }

    bool
    well_controller::is_valid_connection_bhp (item_t pressure, item_t bhp) const
      {
        item_t diff = pressure - bhp;
        if (is_production ())
          return diff > 0;
        else
          return diff < 0;
      }


    void
    well_controller::switch_to_bhp (sp_well_t &well)
    {
      BOSOUT (section::wells, level::medium) << "[" << well->name () << "] switch to bhp" << bs_end;
      BOSOUT (section::wells, level::low)    << "\tlimit_bhp: " << bhp_ << "\twell_bhp: " << well->bhp () << bs_end;

      well->set_bhp (bhp_);
      current_control_.is_bhp = true;
    }

    const injection_type &
    well_controller::injection () const
    {
      return injection_type_;
    }
    void
    well_controller::set_injection_type (injection_type type)
    {
      injection_type_ = type;
    }

    rate_control_type
    well_controller::get_control_type () const
      {
        return current_control_.control_type;
      }

    //////////////////////////////////////////////////////////////////////////

    well_controller_factory::well_controller_factory (bs_type_ctor_param /*param = NULL */)
    {

    }
    well_controller_factory::well_controller_factory (const well_controller_factory &f)
    : bs_refcounter (f), objbase (f)
    {
      *this = f;
    }

    well_controller_factory::sp_well_controller_t
    well_controller_factory::create_controller () const
    {
      return BS_KERNEL.create_object (well_controller_t::bs_type (), true);
    }
    //////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////
    BLUE_SKY_TYPE_STD_CREATE (well_controller);
    BLUE_SKY_TYPE_STD_COPY (well_controller);
    BLUE_SKY_TYPE_IMPL (well_controller, objbase, "well_controller_seq", "Base class for well controllers", "Base class for well controllers");

    BLUE_SKY_TYPE_STD_CREATE (well_controller_factory);
    BLUE_SKY_TYPE_STD_COPY (well_controller_factory);
    BLUE_SKY_TYPE_IMPL (well_controller_factory, objbase, "well_controller_factory_seq", "Base class for factory of well controllers", "Base class for factory of well controllers");
    //////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////
    bool
    well_controller_register_type (const blue_sky::plugin_descriptor &pd)
    {
      bool res = true;

      res &= BS_KERNEL.register_type (pd, well_controller::bs_type ()); BS_ASSERT (res);

      return res;
    }
    bool
    well_controller_factory_register_type (const blue_sky::plugin_descriptor &pd)
    {
      bool res = true;

      res &= BS_KERNEL.register_type (pd, well_controller_factory::bs_type ()); BS_ASSERT (res);

      return res;
    }
    //////////////////////////////////////////////////////////////////////////

  } // namespace wells
} // namespace blue_sky

