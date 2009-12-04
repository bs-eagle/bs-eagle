/**
 *       \file  calc_well.cpp
 *      \brief  Impementation of base well class
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  23.06.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"
#include "calc_well.h"

#include "well_connection.h"
#include "well_controller.h"

#include "calc_model.h"

#include "calc_model_type_helper.h"

#include "reservoir.h"
#include "facility_manager.h"

#include "matrix_inverse.h"

#include "calc_well_pressure.h"
#include "calc_rho.h"
#include "wellbore_density_calc.h"
#include "calc_perf_bhp.h"

#include "exp_temp_mx.h"
#include "apply_wefac.h"
#include "default_rr_eliminator.h"

namespace blue_sky
  {

  ///////////////////////////////////////////////////////////////////////////
  well_state_type
  well_state_cast (const std::string &str)
  {
    if (str == "OPEN")
      return well_open;
    else if (str == "SHUT")
      return well_shut;
    else if (str == "")
      return well_shut;
    else
      {
        BS_ASSERT (false && "Unsupported value of well state") (str);
        throw bs_exception ("well_state_cast", "Unsupported value of well state");
      }
  }
  ///////////////////////////////////////////////////////////////////////////


  //////////////////////////////////////////////////////////////////////////
  template <typename strategy_t>
  well<strategy_t>::~well ()
  {

  }

  template <typename strategy_t>
  well<strategy_t>::well (bs_type_ctor_param /*param = NULL */)
      : calc_well_pressure_ (BS_KERNEL.create_object (calc_well_pressure <strategy_t>::bs_type (), true))
      , calc_rho_ (BS_KERNEL.create_object (calc_total_average_rho <strategy_t>::bs_type (), true))
      , calc_perf_density_ (BS_KERNEL.create_object (wellbore_density_calc <strategy_t>::bs_type (), true))
      , calc_perf_bhp_ (BS_KERNEL.create_object (calc_perf_bhp <strategy_t>::bs_type (), true))
  {
    clear_data ();
    bhp_ = 0;
  }

  template <typename strategy_t>
  well <strategy_t>::well (const std::string &well_name)
  : calc_well_pressure_ (BS_KERNEL.create_object (calc_well_pressure <strategy_t>::bs_type (), true))
  , calc_rho_ (BS_KERNEL.create_object (calc_total_average_rho <strategy_t>::bs_type (), true))
  , calc_perf_density_ (BS_KERNEL.create_object (wellbore_density_calc <strategy_t>::bs_type (), true))
  , calc_perf_bhp_ (BS_KERNEL.create_object (calc_perf_bhp <strategy_t>::bs_type (), true))
  {
    set_name (well_name);

    clear_data ();
    bhp_ = 0;
  }

  template <typename strategy_t>
  well<strategy_t>::well (const well &w)
        : bs_refcounter ()
  {
    *this = w;
  }

  template <typename strategy_t>
  shared_vector <typename strategy_t::item_t>
  well <strategy_t>::get_ww_value () 
  {
    return shared_array <item_t> (0, 0);
  }
  template <typename strategy_t>
  shared_vector <typename strategy_t::item_t>
  well <strategy_t>::get_bw_value () 
  {
    return shared_array <item_t> (0, 0);
  }

  template <typename strategy_t>
  const std::string &
  well <strategy_t>::name () const
  {
    return name_;
  }

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  template <typename strategy_t>
  void
  well<strategy_t>::set_coord (index_t i_coord, index_t j_coord)
  {
    i_coord_ = i_coord;
    j_coord_ = j_coord;
  }
  template <typename strategy_t>
  void
  well<strategy_t>::set_bhp_depth (item_t bhp_depth)
  {
    bhp_depth_ = bhp_depth;
  }
  template <typename strategy_t>
  void
  well<strategy_t>::set_state (well_state_type state, const sp_calc_model_t &calc_model)
  {
    well_state_.state = state;
    if (state == well_shut)
      {
        rate_ = 0;
        shut_well (calc_model);
      }

#ifdef _DEBUG
    if (state == well_shut)
      {
        BOSOUT (section::wells, level::debug) << "[" << name () << "] " << "well shut" << bs_end;
      }
    else
      {
        BOSOUT (section::wells, level::debug) << "[" << name () << "] " << "well open" << bs_end;
      }
#endif

    check_shut ();
  }
  template <typename strategy_t>
  void
  well<strategy_t>::set_exploitation_factor (item_t exploitation_factor)
  {
    exploitation_factor_ = exploitation_factor;
  }

  template <typename strategy_t>
  void
  well<strategy_t>::set_controller (sp_well_controller_t controller)
  {
    well_controller_ = controller;
  }
  template <typename strategy_t>
  void
  well<strategy_t>::set_limit_operation (sp_well_limit_operation_t limit_operation)
  {
    well_limit_operation_ = limit_operation;
  }
  template <typename strategy_t>
  typename well<strategy_t>::sp_well_controller_t
  well<strategy_t>::get_controller () const
    {
      return well_controller_;
    }
  template <typename strategy_t>
  typename well<strategy_t>::sp_well_limit_operation_t
  well<strategy_t>::get_limit_operation () const
    {
      return well_limit_operation_;
    }
#if 0
  template <typename strategy_t>
  typename well<strategy_t>::ijk_filter_t
  well<strategy_t>::filter_connections (index_t /*i_coord*/, index_t /*j_coord*/, index_t /*k*/, 
                                        index_t /*z1*/, index_t /*z2*/)
  {
    return ijk_filter_t (connection_list_, ijk_pred_t (0));  // TODO:
  }
  template <typename strategy_t>
  typename well<strategy_t>::ijk1k2_filter_t
  well<strategy_t>::filter_connections (index_t /*i_coord*/, index_t /*j_coord*/, index_t /*kw1*/, 
                                        index_t /*kw2*/)
  {
    return ijk1k2_filter_t (connection_list_, ijk1k2_pred_t (0));  // TODO:
  }
#endif //0
  template <typename strategy_t>
  const std::string &
  well<strategy_t>::get_name () const
    {
      return name_;
    }
  template <typename strategy_t>
  void
  well<strategy_t>::set_name (const std::string &name)
  {
    name_ = name;
  }

  template <typename strategy_t>
  void
  well<strategy_t>::compute_connection_factor (const physical_constants &internal_constants,
      const sp_params_t &params,
      const sp_mesh_iface_t &mesh,
      const item_array_t &perm,
      const item_array_t &ntg,
      bool ro_calc_flag)
  {
    BS_ASSERT (!is_no_connections ()) (name ());

    connection_iterator_t it = connections_begin (), e = connections_end ();
    for (; it != e; ++it)
      {
        const sp_connection_t &c (*it);

        if (!c->is_shut ())
          c->compute_factors (internal_constants, params, mesh, perm, ntg, ro_calc_flag);
      }
  }

  template <typename strategy_t>
  bool
  well<strategy_t>::fi_check_limits () const
    {
      BS_ASSERT (well_limit_operation_);
      if (!well_limit_operation_)
        {
          return false;
          //throw bs_exception ("well::fi_check_limits", "well_limit_operation is null");
        }

      return well_limit_operation_->fi_check_limits ();
    }

  template <typename strategy_t>
  bool
  well <strategy_t>::check_connections_bhp (const item_array_t &pressure) const
    {
      BS_ASSERT (!is_shut ()) (name ());
      if (is_no_connections ())
        {
          BOSOUT (section::wells, level::debug)
            << "[" << name_ << "] check_connections_bhp: connection list is empty"
            << bs_end;

          return false;
        }

      connection_iterator_t it = connections_begin (), e = connections_end ();
      for (; it != e; ++it)
        {
          const sp_connection_t &c (*it);
          if (c->is_shut ())
            continue;

#ifdef _DEBUG
          index_t n_block = c->n_block ();
          BOSOUT (section::wells, level::debug) << "[" << name_ << ": " << n_block << "] " << "bulkp: " << c->bulkp << " cur_bhp: " << c->cur_bhp << bs_end;
#endif
          if (well_controller_->is_valid_connection_bhp (c->bulkp, c->cur_bhp))
            {
              BOSOUT (section::wells, level::low) << "[" << name_ << "] " << " (work)" << bs_end;
              return true;
            }
        }

      BOSOUT (section::wells, level::low) << "[" << name_ << "] " << " (doesn't work)" << bs_end;
      return false;
    }

  template <typename strategy_t>
  void
  well<strategy_t>::shut_well (const sp_calc_model_t &calc_model)
  {
#ifdef _DEBUG
    BOSOUT (section::wells, level::debug) << boost::format ("[%s]: shut") % name () << bs_end;
#endif

    well_state_.state = well_shut;
    rate_ = 0;

    const item_array_t &pressure = calc_model->pressure;
    connection_iterator_t it = connections_begin (), e = connections_end ();
    for (; it != e; ++it)
      {
        const sp_connection_t &c (*it);
        index_t n_block = c->n_block ();
        item_t bulkp = pressure[n_block];

        c->set_cur_bhp (bulkp);
        c->set_bulkp (bulkp);
      }
  }

  template <typename strategy_t>
  bool
  well <strategy_t>::is_shut () const
    {
      return well_state_.state == well_shut;
    }

  template <typename strategy_t>
  typename well<strategy_t>::item_t
  well<strategy_t>::get_reference_depth (const sp_mesh_iface_t &mesh) const
  {
    BS_ASSERT (!is_no_connections ()) (name ());

    if (is_no_primary_connections ())
      return input_reference_depth_;

    const sp_connection_t &c = get_first_connection ();
    wells::connection_direction_type dir = c->get_dir ();
    item_t dtop = 0;
    if (dir == wells::direction_x || dir == wells::direction_y)
      {
        dtop = c->get_connection_depth ();
      }
    else
      {
        dtop = mesh->get_element_dtop (c->n_block ());
      }

    return input_reference_depth_ > 0 ? input_reference_depth_ : dtop;
  }

  template <typename strategy_t>
  typename well<strategy_t>::item_t
  well<strategy_t>::get_reference_pressure () const
    {
      return bhp_;
    }

  template <typename strategy_t>
  void
  well<strategy_t>::set_bhp (item_t bhp)
  {
    bhp_ = bhp;
  }

  template <typename strategy_t>
  bool
  well<strategy_t>::is_bhp () const
    {
      BS_ASSERT (well_controller_);
      return well_controller_->is_bhp ();
    }
  template <typename strategy_t>
  bool
  well<strategy_t>::is_rate () const
    {
      BS_ASSERT (well_controller_);
      return well_controller_->is_rate ();
    }

  template <typename strategy_t>
  const typename well<strategy_t>::sp_well_controller_t &
  well<strategy_t>::get_well_controller () const
    {
      return well_controller_;
    }

  //////////////////////////////////////////////////////////////////////////

  template <typename strategy_t>
  typename well<strategy_t>::item_t
  well <strategy_t>::get_input_rate () const
    {
      BS_ASSERT (well_controller_);
      using namespace wells;
      rate_control_type control_type = well_controller_->get_control_type ();
      injection_type inj_type = well_controller_->injection ();
      bool is_prod = well_controller_->is_production ();

      if (is_prod)
        {
          if (control_type == liquid_rate_control)
            return -well_controller_->rate_.prod.liquid;
          else if (control_type == water_rate_control)
            return -well_controller_->rate_.prod.water;
          else if (control_type == gas_rate_control)
            return -well_controller_->rate_.prod.gas;
          else if (control_type == oil_rate_control)
            return -well_controller_->rate_.prod.oil;
        }
      else
        {
          if (is_rate ())
            {
              if (inj_type == injection_water)
                return well_controller_->rate_.inj.water;
              else if (inj_type == injection_gas)
                return well_controller_->rate_.inj.gas;
              else if (inj_type == injection_oil)
                return well_controller_->rate_.inj.oil;
            }
          else
            return 0;
        }

      BS_ASSERT (false && "get_input_rate: return 0") (is_prod) (control_type) (inj_type);
      return 0;
    }

  template <typename strategy_t>
  void
  well <strategy_t>::reset_init_approx ()
  {
    init_approx_is_calc_ = false;
  }

  template <typename strategy_t>
  void
  well<strategy_t>::process (bool /*is_start*/, double /*dt*/, const sp_calc_model_t &/*calc_model*/, const sp_mesh_iface_t &/*mesh*/, sp_jmatrix_t &/*jmatrix*/)
  {
    bs_throw_exception ("PURE_CALL");
  }
  template <typename strategy_t>
  void
  well <strategy_t>::clear_data ()
  {
    rate_     = 0;
    rate_rc_  = 0;
    gor_      = 0;
  }

  template <typename strategy_t>
  void 
  well <strategy_t>::restore_solution (double /*dt*/, const item_array_t &/*p_sol*/, const item_array_t &/*s_sol*/, index_t /*block_size*/)
  {
    bs_throw_exception ("PURE CALL");
  }

  template <typename strategy_t>
  void
  well <strategy_t>::custom_init (const sp_calc_model_t &/*mdl*/)
  {
  }

  template <typename strategy_t>
  void
  well <strategy_t>::pre_large_step (const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh)
  {
    compute_connection_factor (calc_model->internal_constants, 
      calc_model->ts_params, 
      mesh,
      calc_model->rock_grid_prop->permeability,
      calc_model->rock_grid_prop->net_to_gros,
      false);

    reset_init_approx ();
    check_shut ();
    custom_init (calc_model);
  }
  template <typename strategy_t>
  void
  well<strategy_t>::pre_small_step ()
  {
    BS_ASSERT (well_controller_);

    saved_well_state_ = well_state_;
    well_controller_->save_control ();
  }
  template <typename strategy_t>
  void
  well<strategy_t>::pre_newton_step ()
  {
    BS_ASSERT (well_controller_);

    saved_niter_well_state_ = well_state_;
    well_controller_->save_niter_control ();
  }
  template <typename strategy_t>
  void
  well<strategy_t>::restart_small_step ()
  {
    BS_ASSERT (well_controller_);

    well_state_ = saved_well_state_;
    init_approx_is_calc_ = well_controller_->restore_control ();
  }
  template <typename strategy_t>
  void
  well<strategy_t>::restart_newton_step ()
  {
    BS_ASSERT (well_controller_);

    well_state_ = saved_niter_well_state_;
    init_approx_is_calc_ = well_controller_->restore_niter_control ();
  }

  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE_T_DEF (well, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF (well, (class));
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (well<base_strategy_fi>), 1, (facility_base<base_strategy_fi>), "well_seq_fi", "well_seq_fi", "well_seq_fi", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (well<base_strategy_di>), 1, (facility_base<base_strategy_di>), "well_seq_di", "well_seq_di", "well_seq_di", false);
  BLUE_SKY_TYPE_IMPL_T_EXT (1, (well<base_strategy_mixi>), 1, (facility_base<base_strategy_mixi>), "well_seq_mixi", "well_seq_mixi", "well_seq_mixi", false);
  //////////////////////////////////////////////////////////////////////////


  //////////////////////////////////////////////////////////////////////////
  bool
  calc_well_register_types (const blue_sky::plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, wells::connection<base_strategy_fi>::bs_type ());
    BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, wells::connection<base_strategy_di>::bs_type ());
    BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, wells::connection<base_strategy_mixi>::bs_type ());
    BS_ASSERT (res);

    res &= BS_KERNEL.register_type (pd, well<base_strategy_fi>::bs_type ());
    BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, well<base_strategy_di>::bs_type ());
    BS_ASSERT (res);

    return res;
  }

} // namespace blue_sky
