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

  template <typename strategy_t>
  void
  well_event <strategy_t>::apply (const sp_top &top, const sp_mesh_iface_t &mesh,
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
  template <typename strategy_t>
  void
  well_event <strategy_t>::apply_internal (const sp_top &top, const sp_mesh_iface_t &mesh,
                           const sp_calc_model_t &calc_model, const smart_ptr <idata, true> &data) const
  {
    BS_ASSERT (false && "BASE METHOD CALL");
  }

  template <typename strategy_t>
  std::string
  well_event <strategy_t>::get_well_name () const
    {
      BS_ASSERT (false && "BASE METHOD CALL");
      return "BASE METHOD CALL";
    }
  template <typename strategy_t>
  std::string
  well_event <strategy_t>::get_group_name () const
    {
      BS_ASSERT (false && "BASE METHOD CALL");
      return "BASE METHOD CALL";
    }
  template <typename strategy_t>
  std::string
  well_event <strategy_t>::get_event_name () const
    {
      BS_ASSERT (false && "BASE METHOD CALL");
      return "BASE METHOD CALL";
    }


  //constructors
  template <typename strategy_t>
  WELSPECS_event<strategy_t>::WELSPECS_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  template <typename strategy_t>
  WELSPECS_event<strategy_t>::WELSPECS_event(const WELSPECS_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }

//========================================================================

  template <typename strategy_t>
  WELLCON_event<strategy_t>::WELLCON_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  template <typename strategy_t>
  WELLCON_event<strategy_t>::WELLCON_event(const WELLCON_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }

//========================================================================

  template <typename strategy_t>
  COMPDAT_event<strategy_t>::COMPDAT_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  template <typename strategy_t>
  COMPDAT_event<strategy_t>::COMPDAT_event(const COMPDAT_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }

//========================================================================

  template <typename strategy_t>
  WCONPROD_event<strategy_t>::WCONPROD_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  template <typename strategy_t>
  WCONPROD_event<strategy_t>::WCONPROD_event(const WCONPROD_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }

//========================================================================

  template <typename strategy_t>
  WCONHIST_event<strategy_t>::WCONHIST_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  template <typename strategy_t>
  WCONHIST_event<strategy_t>::WCONHIST_event(const WCONHIST_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }

//========================================================================

  template <typename strategy_t>
  WCONINJE_event<strategy_t>::WCONINJE_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  template <typename strategy_t>
  WCONINJE_event<strategy_t>::WCONINJE_event(const WCONINJE_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }

  //========================================================================

  template <typename strategy_t>
  WECON_event<strategy_t>::WECON_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  template <typename strategy_t>
  WECON_event<strategy_t>::WECON_event(const WECON_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }

  //========================================================================

  template <typename strategy_t>
  WECONINJ_event<strategy_t>::WECONINJ_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  template <typename strategy_t>
  WECONINJ_event<strategy_t>::WECONINJ_event(const WECONINJ_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }

  //========================================================================

  template <typename strategy_t>
  WEFAC_event<strategy_t>::WEFAC_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  template <typename strategy_t>
  WEFAC_event<strategy_t>::WEFAC_event(const WEFAC_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }

  //========================================================================

  template <typename strategy_t>
  WELTARG_event<strategy_t>::WELTARG_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  template <typename strategy_t>
  WELTARG_event<strategy_t>::WELTARG_event(const WELTARG_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }


//========================================================================

  template <typename strategy_t>
  WPIMULT_event<strategy_t>::WPIMULT_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  template <typename strategy_t>
  WPIMULT_event<strategy_t>::WPIMULT_event(const WPIMULT_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }

//========================================================================

  template <typename strategy_t>
  COMPENSATION_event<strategy_t>::COMPENSATION_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  template <typename strategy_t>
  COMPENSATION_event<strategy_t>::COMPENSATION_event(const COMPENSATION_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }

//========================================================================

  template <typename strategy_t>
  PERMFRAC_event<strategy_t>::PERMFRAC_event(bs_type_ctor_param param)
  : main_params_ (BS_KERNEL.create_object (main_params_class::bs_type ()))
  {
  }
  template <typename strategy_t>
  PERMFRAC_event<strategy_t>::PERMFRAC_event(const PERMFRAC_event& src)
  : bs_refcounter (src)
  {
    *this = src;
  }



  //bs stuff
  BLUE_SKY_TYPE_STD_CREATE_T_DEF (WELSPECS_event, (class))
  BLUE_SKY_TYPE_STD_COPY_T_DEF (WELSPECS_event, (class))
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WELSPECS_event <base_strategy_fi>), 1, (well_event <base_strategy_fi>), "WELSPECS_seq_fi", "WELSPECS", "BOS_Core WELSPECS_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WELSPECS_event <base_strategy_di>), 1, (well_event <base_strategy_di>), "WELSPECS_seq_di", "WELSPECS", "BOS_Core WELSPECS_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WELSPECS_event <base_strategy_mixi>), 1, (well_event <base_strategy_mixi>), "WELSPECS_seq_mixi", "WELSPECS", "BOS_Core WELSPECS_event class", false)

  BLUE_SKY_TYPE_STD_CREATE_T_DEF (WELLCON_event, (class))
  BLUE_SKY_TYPE_STD_COPY_T_DEF (WELLCON_event, (class))
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WELLCON_event <base_strategy_fi>), 1, (well_event <base_strategy_fi>), "WELLCON_seq_fi", "WELLCON", "BOS_Core WELLCON_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WELLCON_event <base_strategy_di>), 1, (well_event <base_strategy_di>), "WELLCON_seq_di", "WELLCON", "BOS_Core WELLCON_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WELLCON_event <base_strategy_mixi>), 1, (well_event <base_strategy_mixi>), "WELLCON_seq_mixi", "WELLCON", "BOS_Core WELLCON_event class", false)

  BLUE_SKY_TYPE_STD_CREATE_T_DEF (COMPDAT_event, (class))
  BLUE_SKY_TYPE_STD_COPY_T_DEF (COMPDAT_event, (class))
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (COMPDAT_event <base_strategy_fi>), 1, (well_event <base_strategy_fi>), "COMPDAT_seq_fi", "COMPDAT", "BOS_Core COMPDAT_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (COMPDAT_event <base_strategy_di>), 1, (well_event <base_strategy_di>), "COMPDAT_seq_di", "COMPDAT", "BOS_Core COMPDAT_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (COMPDAT_event <base_strategy_mixi>), 1, (well_event <base_strategy_mixi>), "COMPDAT_seq_mixi", "COMPDAT", "BOS_Core COMPDAT_event class", false)

  BLUE_SKY_TYPE_STD_CREATE_T_DEF (WCONPROD_event, (class))
  BLUE_SKY_TYPE_STD_COPY_T_DEF (WCONPROD_event, (class))
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WCONPROD_event <base_strategy_fi>), 1, (well_event <base_strategy_fi>), "WCONPROD_seq_fi", "WCONPROD", "BOS_Core WCONPROD_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WCONPROD_event <base_strategy_di>), 1, (well_event <base_strategy_di>), "WCONPROD_seq_di", "WCONPROD", "BOS_Core WCONPROD_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WCONPROD_event <base_strategy_mixi>), 1, (well_event <base_strategy_mixi>), "WCONPROD_seq_mixi", "WCONPROD", "BOS_Core WCONPROD_event class", false)

  BLUE_SKY_TYPE_STD_CREATE_T_DEF (WCONHIST_event, (class))
  BLUE_SKY_TYPE_STD_COPY_T_DEF (WCONHIST_event, (class))
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WCONHIST_event <base_strategy_fi>), 1, (well_event <base_strategy_fi>), "WCONHIST_seq_fi", "WCONHIST", "BOS_Core WCONHIST_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WCONHIST_event <base_strategy_di>), 1, (well_event <base_strategy_di>), "WCONHIST_seq_di", "WCONHIST", "BOS_Core WCONHIST_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WCONHIST_event <base_strategy_mixi>), 1, (well_event <base_strategy_mixi>), "WCONHIST_seq_mixi", "WCONHIST", "BOS_Core WCONHIST_event class", false)

  BLUE_SKY_TYPE_STD_CREATE_T_DEF (WCONINJE_event, (class))
  BLUE_SKY_TYPE_STD_COPY_T_DEF (WCONINJE_event, (class))
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WCONINJE_event <base_strategy_fi>), 1, (well_event <base_strategy_fi>), "WCONINJE_seq_fi", "WCONINJE", "BOS_Core WCONINJE_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WCONINJE_event <base_strategy_di>), 1, (well_event <base_strategy_di>), "WCONINJE_seq_di", "WCONINJE", "BOS_Core WCONINJE_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WCONINJE_event <base_strategy_mixi>), 1, (well_event <base_strategy_mixi>), "WCONINJE_seq_mixi", "WCONINJE", "BOS_Core WCONINJE_event class", false)

  BLUE_SKY_TYPE_STD_CREATE_T_DEF (WECON_event, (class))
  BLUE_SKY_TYPE_STD_COPY_T_DEF (WECON_event, (class))
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WECON_event <base_strategy_fi>), 1, (well_event <base_strategy_fi>), "WECON_seq_fi", "WECON", "BOS_Core WECON_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WECON_event <base_strategy_di>), 1, (well_event <base_strategy_di>), "WECON_seq_di", "WECON", "BOS_Core WECON_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WECON_event <base_strategy_mixi>), 1, (well_event <base_strategy_mixi>), "WECON_seq_mixi", "WECON", "BOS_Core WECON_event class", false)

  BLUE_SKY_TYPE_STD_CREATE_T_DEF (WECONINJ_event, (class))
  BLUE_SKY_TYPE_STD_COPY_T_DEF (WECONINJ_event, (class))
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WECONINJ_event <base_strategy_fi>), 1, (well_event <base_strategy_fi>), "WECONINJ_seq_fi", "WECONINJ", "BOS_Core WECONINJ_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WECONINJ_event <base_strategy_di>), 1, (well_event <base_strategy_di>), "WECONINJ_seq_di", "WECONINJ", "BOS_Core WECONINJ_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WECONINJ_event <base_strategy_mixi>), 1, (well_event <base_strategy_mixi>), "WECONINJ_seq_mixi", "WECONINJ", "BOS_Core WECONINJ_event class", false)

  BLUE_SKY_TYPE_STD_CREATE_T_DEF (WEFAC_event, (class))
  BLUE_SKY_TYPE_STD_COPY_T_DEF (WEFAC_event, (class))
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WEFAC_event <base_strategy_fi>), 1, (well_event <base_strategy_fi>), "WEFAC_seq_fi", "WEFAC", "BOS_Core WEFAC_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WEFAC_event <base_strategy_di>), 1, (well_event <base_strategy_di>), "WEFAC_seq_di", "WEFAC", "BOS_Core WEFAC_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WEFAC_event <base_strategy_mixi>), 1, (well_event <base_strategy_mixi>), "WEFAC_seq_mixi", "WEFAC", "BOS_Core WEFAC_event class", false)

  BLUE_SKY_TYPE_STD_CREATE_T_DEF (WELTARG_event, (class))
  BLUE_SKY_TYPE_STD_COPY_T_DEF (WELTARG_event, (class))
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WELTARG_event <base_strategy_fi>), 1, (well_event <base_strategy_fi>), "WELTARG_seq_fi", "WELTARG", "BOS_Core WELTARG_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WELTARG_event <base_strategy_di>), 1, (well_event <base_strategy_di>), "WELTARG_seq_di", "WELTARG", "BOS_Core WELTARG_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WELTARG_event <base_strategy_mixi>), 1, (well_event <base_strategy_mixi>), "WELTARG_seq_mixi", "WELTARG", "BOS_Core WELTARG_event class", false)

  BLUE_SKY_TYPE_STD_CREATE_T_DEF (WPIMULT_event, (class))
  BLUE_SKY_TYPE_STD_COPY_T_DEF (WPIMULT_event, (class))
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WPIMULT_event <base_strategy_fi>), 1, (well_event <base_strategy_fi>), "WPIMULT_seq_fi", "WPIMULT", "BOS_Core WPIMULT_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WPIMULT_event <base_strategy_di>), 1, (well_event <base_strategy_di>), "WPIMULT_seq_di", "WPIMULT", "BOS_Core WPIMULT_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (WPIMULT_event <base_strategy_mixi>), 1, (well_event <base_strategy_mixi>), "WPIMULT_seq_mixi", "WPIMULT", "BOS_Core WPIMULT_event class", false)

  BLUE_SKY_TYPE_STD_CREATE_T_DEF (COMPENSATION_event, (class))
  BLUE_SKY_TYPE_STD_COPY_T_DEF (COMPENSATION_event, (class))
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (COMPENSATION_event <base_strategy_fi>), 1, (well_event <base_strategy_fi>), "COMPENSATION_seq_fi", "COMPENSATION", "BOS_Core COMPENSATION_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (COMPENSATION_event <base_strategy_di>), 1, (well_event <base_strategy_di>), "COMPENSATION_seq_di", "COMPENSATION", "BOS_Core COMPENSATION_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (COMPENSATION_event <base_strategy_mixi>), 1, (well_event <base_strategy_mixi>), "COMPENSATION_seq_mixi", "COMPENSATION", "BOS_Core COMPENSATION_event class", false)

  BLUE_SKY_TYPE_STD_CREATE_T_DEF (PERMFRAC_event, (class))
  BLUE_SKY_TYPE_STD_COPY_T_DEF (PERMFRAC_event, (class))
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (PERMFRAC_event <base_strategy_fi>), 1, (well_event <base_strategy_fi>), "PERMFRAC_seq_fi", "PERMFRAC", "BOS_Core PERMFRAC_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (PERMFRAC_event <base_strategy_di>), 1, (well_event <base_strategy_di>), "PERMFRAC_seq_di", "PERMFRAC", "BOS_Core PERMFRAC_event class", false)
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (PERMFRAC_event <base_strategy_mixi>), 1, (well_event <base_strategy_mixi>), "PERMFRAC_seq_mixi", "PERMFRAC", "BOS_Core PERMFRAC_event class", false)


}//ns bs
