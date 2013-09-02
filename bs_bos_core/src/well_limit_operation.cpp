/**
 * \file well_limit_operation.cpp
 * \brief impl of
 * \author Sergey Miryanov
 * \date 16.07.2008
 * */
#include "stdafx.h"

#include "well_limit_operation.h"

namespace blue_sky
  {
  namespace wells
    {

    limit_operation_type
    limit_operation_cast (const std::string &str)
    {
      if (str == "NONE")
        return operation_none;
      else if (str == "CON")
        return operation_con;
      else if (str == "WELL")
        return operation_well;
      else if (str == "")
        {
          BS_ASSERT (false && "limit_operation_type is empty. will be used operation_none");
          return operation_none;
        }
      else
        {
          BS_ASSERT (false && "Unsupported value of limit_operation type") (str);
          throw bs_exception ("limit_operation_cast", "Unsupported value of limit_operation type");
        }
    }

    well_limit_operation::well_limit_operation (bs_type_ctor_param /*param = NULL */)
        : min_oil_rate_ (0)
        , max_water_cut_ (0)
        , min_water_cut_ (0)
    {

    }
    well_limit_operation::well_limit_operation (const well_limit_operation &op)
    : bs_refcounter (op), objbase (op)
    {
      min_oil_rate_		= op.min_oil_rate_;
      max_water_cut_	= op.max_water_cut_;
      min_water_cut_	= op.min_water_cut_;
    }

    void
    well_limit_operation::set_min_oil_rate (item_t min_oil_rate)
    {
      min_oil_rate_ = min_oil_rate;
    }
    void
    well_limit_operation::set_max_water_cut (item_t max_water_cut)
    {
      max_water_cut_ = max_water_cut;
    }
    void
    well_limit_operation::set_min_water_cut (item_t min_water_cut)
    {
      min_water_cut_ = min_water_cut;
    }

    void
    well_limit_operation::set_value (int /*value_type*/, item_t /*value*/)
    {
      BS_ASSERT (false && "NOT IMPL YET");
    }

    bool
    well_limit_operation::fi_check_limits () const
      {
        BS_ASSERT (false && "NOT IMPL YET");
        return false;

        //if (current_status < 0)
        //  {
        //    if (current_rate.oil > 0.1 ||
        //        current_rate.water > 0.1 ||
        //        current_rate.gas > 0.1)
        //      {
        //        rep->print (LOG_ITERS_SECTION, LOG_ERR,
        //                    "***---------------------------------------------------------------------***\n");
        //        rep->print (LOG_ITERS_SECTION, LOG_ERR,
        //                    GET_TEXT ("Warning: production well %s could not operate under current conditions\n"),
        //                    get_name ());
        //        print_current_status ();
        //        fi_close_well ();
        //        rep->print (LOG_ITERS_SECTION, LOG_ERR,
        //                    "***---------------------------------------------------------------------***\n");
        //        return 1;
        //      }
        //  }
        //// injection well
        //else if (current_status > 0)
        //  {
        //    if (current_rate.oil_inj < -0.1 ||
        //        current_rate.water_inj < -0.1 ||
        //        current_rate.gas_inj < -0.1)
        //      {
        //        rep->print (LOG_ITERS_SECTION, LOG_ERR,
        //                    "***---------------------------------------------------------------------***\n");
        //        rep->print (LOG_ITERS_SECTION, LOG_ERR,
        //                    GET_TEXT ("Warning: injection well %s could not operate under current conditions\n"),
        //                    get_name ());
        //        print_current_status ();
        //        fi_close_well ();
        //        rep->print (LOG_ITERS_SECTION, LOG_ERR,
        //                    "***---------------------------------------------------------------------***\n");
        //        return 1;
        //      }
        //  }
        //return 0;

      }

    //////////////////////////////////////////////////////////////////////////

    well_limit_operation_factory::well_limit_operation_factory (bs_type_ctor_param /*param = NULL */)
    {

    }
    well_limit_operation_factory::well_limit_operation_factory (const well_limit_operation_factory &f)
    : bs_refcounter (f), objbase (f)
    {

    }

    well_limit_operation_factory::sp_well_limit_operation_t
    well_limit_operation_factory::create_limit (limit_operation_type /*type*/) const
      {
        BS_ASSERT (false && "BASE METHOD CALL");
        return BS_KERNEL.create_object (well_limit_operation::bs_type ());
      }

    //////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////
    BLUE_SKY_TYPE_STD_CREATE (well_limit_operation);
    BLUE_SKY_TYPE_STD_COPY (well_limit_operation);
    BLUE_SKY_TYPE_IMPL (well_limit_operation, objbase, "well_limit_operation", "Base class for well limit operations", "Base class for well limit operation");

    BLUE_SKY_TYPE_STD_CREATE (well_limit_operation_factory);
    BLUE_SKY_TYPE_STD_COPY (well_limit_operation_factory);
    BLUE_SKY_TYPE_IMPL (well_limit_operation_factory, objbase, "well_limit_operation_factory", "Base class for factory of well limit operations", "Base class for factory of well limit operations");
    //////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////
    bool
    well_limit_operation_register_type (const blue_sky::plugin_descriptor &pd)
    {
      bool res = true;

      res &= BS_KERNEL.register_type (pd, well_limit_operation::bs_type ());
      BS_ASSERT (res);

      return res;
    }
    bool
    well_limit_operation_factory_register_type (const blue_sky::plugin_descriptor &pd)
    {
      bool res = true;

      res &= BS_KERNEL.register_type (pd, well_limit_operation_factory::bs_type ());
      BS_ASSERT (res);

      return res;
    }
    //////////////////////////////////////////////////////////////////////////

  }	// namespace wells
}	// namespace blue_sky
