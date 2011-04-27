/**
 *       \file  well_events.cpp
 *      \brief  Constructors for well events
 *     \author  Morozov Andrey
 *       \date  07.06.2008
 *  \copyright  This source code is released under the terms of
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"

#include "well_events.h"
#include "calc_model.h"
#include "event_filter.h"
#include "reservoir.h"
#include "facility_manager.h"
#include "well_connection.h"

#define EREGF(class_name,name,descr,type) ereg(name, #name, descr, type)

namespace blue_sky
  {

  void
  well_event::apply (const sp_top &top, const sp_mesh_iface_t &mesh,
                           const sp_calc_model_t &calc_model, const smart_ptr <idata, true> &data) const
  {
    if (top->get_event_filter ()->accept_well (get_well_name ()))
      {
        apply_internal (top, mesh, calc_model, data);
      }
    else
      {
        BOSOUT (section::schedule, level::low) << "[" << get_well_name () << "] reject well event " << get_event_name () << bs_end;
      }
  }
  void
  well_event::apply_internal (const sp_top &top, const sp_mesh_iface_t &mesh,
                           const sp_calc_model_t &calc_model, const smart_ptr <idata, true> &data) const
  {
    BS_ASSERT (false && "BASE METHOD CALL");
  }

  std::string
  well_event::get_well_name () const
    {
      BS_ASSERT (false && "BASE METHOD CALL");
      return "BASE METHOD CALL";
    }
  std::string
  well_event::get_group_name () const
    {
      BS_ASSERT (false && "BASE METHOD CALL");
      return "BASE METHOD CALL";
    }
  std::string
  well_event::get_event_name () const
    {
      BS_ASSERT (false && "BASE METHOD CALL");
      return "BASE METHOD CALL";
    }


  //constructors
  WELSPECS_event::WELSPECS_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  WELSPECS_event::WELSPECS_event(const WELSPECS_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }

//========================================================================

  WELLCON_event::WELLCON_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  WELLCON_event::WELLCON_event(const WELLCON_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }

//========================================================================

  COMPDAT_event::COMPDAT_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  COMPDAT_event::COMPDAT_event(const COMPDAT_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }

//========================================================================

  WCONPROD_event::WCONPROD_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  WCONPROD_event::WCONPROD_event(const WCONPROD_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }

//========================================================================

  WCONHIST_event::WCONHIST_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  WCONHIST_event::WCONHIST_event(const WCONHIST_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }

//========================================================================

  WCONINJE_event::WCONINJE_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  WCONINJE_event::WCONINJE_event(const WCONINJE_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }

  //========================================================================

  WECON_event::WECON_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  WECON_event::WECON_event(const WECON_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }

  //========================================================================

  WECONINJ_event::WECONINJ_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  WECONINJ_event::WECONINJ_event(const WECONINJ_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }

  //========================================================================

  WEFAC_event::WEFAC_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  WEFAC_event::WEFAC_event(const WEFAC_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }

  //========================================================================

  WELTARG_event::WELTARG_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  WELTARG_event::WELTARG_event(const WELTARG_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }


//========================================================================

  WPIMULT_event::WPIMULT_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  WPIMULT_event::WPIMULT_event(const WPIMULT_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }

//========================================================================

  COMPENSATION_event::COMPENSATION_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  COMPENSATION_event::COMPENSATION_event(const COMPENSATION_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }

//========================================================================

  PERMFRAC_event::PERMFRAC_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  PERMFRAC_event::PERMFRAC_event(const PERMFRAC_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }



  //bs stuff
  BLUE_SKY_TYPE_STD_CREATE (WELSPECS_event)
  BLUE_SKY_TYPE_STD_COPY (WELSPECS_event)
  BLUE_SKY_TYPE_IMPL (WELSPECS_event, well_event, "WELSPECS", "WELSPECS", "BOS_Core WELSPECS_event class")

  BLUE_SKY_TYPE_STD_CREATE (WELLCON_event)
  BLUE_SKY_TYPE_STD_COPY (WELLCON_event)
  BLUE_SKY_TYPE_IMPL (WELLCON_event, well_event, "WELLCON", "WELLCON", "BOS_Core WELLCON_event class")

  BLUE_SKY_TYPE_STD_CREATE (COMPDAT_event)
  BLUE_SKY_TYPE_STD_COPY (COMPDAT_event)
  BLUE_SKY_TYPE_IMPL (COMPDAT_event, well_event, "COMPDAT", "COMPDAT", "BOS_Core COMPDAT_event class")

  BLUE_SKY_TYPE_STD_CREATE (WCONPROD_event)
  BLUE_SKY_TYPE_STD_COPY (WCONPROD_event)
  BLUE_SKY_TYPE_IMPL (WCONPROD_event, well_event, "WCONPROD", "WCONPROD", "BOS_Core WCONPROD_event class")

  BLUE_SKY_TYPE_STD_CREATE (WCONHIST_event)
  BLUE_SKY_TYPE_STD_COPY (WCONHIST_event)
  BLUE_SKY_TYPE_IMPL (WCONHIST_event, well_event, "WCONHIST", "WCONHIST", "BOS_Core WCONHIST_event class")

  BLUE_SKY_TYPE_STD_CREATE (WCONINJE_event)
  BLUE_SKY_TYPE_STD_COPY (WCONINJE_event)
  BLUE_SKY_TYPE_IMPL (WCONINJE_event, well_event, "WCONINJE", "WCONINJE", "BOS_Core WCONINJE_event class")

  BLUE_SKY_TYPE_STD_CREATE (WECON_event)
  BLUE_SKY_TYPE_STD_COPY (WECON_event)
  BLUE_SKY_TYPE_IMPL (WECON_event, well_event, "WECON", "WECON", "BOS_Core WECON_event class")

  BLUE_SKY_TYPE_STD_CREATE (WECONINJ_event)
  BLUE_SKY_TYPE_STD_COPY (WECONINJ_event)
  BLUE_SKY_TYPE_IMPL (WECONINJ_event, well_event, "WECONINJ", "WECONINJ", "BOS_Core WECONINJ_event class")

  BLUE_SKY_TYPE_STD_CREATE (WEFAC_event)
  BLUE_SKY_TYPE_STD_COPY (WEFAC_event)
  BLUE_SKY_TYPE_IMPL (WEFAC_event, well_event, "WEFAC", "WEFAC", "BOS_Core WEFAC_event class")

  BLUE_SKY_TYPE_STD_CREATE (WELTARG_event)
  BLUE_SKY_TYPE_STD_COPY (WELTARG_event)
  BLUE_SKY_TYPE_IMPL (WELTARG_event, well_event, "WELTARG", "WELTARG", "BOS_Core WELTARG_event class")

  BLUE_SKY_TYPE_STD_CREATE (WPIMULT_event)
  BLUE_SKY_TYPE_STD_COPY (WPIMULT_event)
  BLUE_SKY_TYPE_IMPL (WPIMULT_event, well_event, "WPIMULT", "WPIMULT", "BOS_Core WPIMULT_event class")

  BLUE_SKY_TYPE_STD_CREATE (COMPENSATION_event)
  BLUE_SKY_TYPE_STD_COPY (COMPENSATION_event)
  BLUE_SKY_TYPE_IMPL (COMPENSATION_event, well_event, "COMPENSATION", "COMPENSATION", "BOS_Core COMPENSATION_event class")

  BLUE_SKY_TYPE_STD_CREATE (PERMFRAC_event)
  BLUE_SKY_TYPE_STD_COPY (PERMFRAC_event)
  BLUE_SKY_TYPE_IMPL (PERMFRAC_event, well_event, "PERMFRAC", "PERMFRAC", "BOS_Core PERMFRAC_event class")


}//ns bs
