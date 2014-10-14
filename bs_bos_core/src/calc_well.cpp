/**
 *       \file  calc_well.cpp
 *      \brief  Impementation of base well class
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  23.06.2008
 *  \copyright  This source code is released under the terms of
 *              the BSD License. See LICENSE for more details.
 * */
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
  well::~well ()
  {

  }

  well::well (bs_type_ctor_param /*param = NULL */)
  : well_events_init_ (this)
  , i_coord_ (-1)
  , j_coord_ (-1)
  , bhp_depth_ (0)
  , exploitation_factor_ (1)
  , init_approx_is_calc_ (false)
  , input_reference_depth_ (0)
  , calc_well_pressure_ (BS_KERNEL.create_object (calc_well_pressure::bs_type (), true))
  , calc_rho_ (BS_KERNEL.create_object (calc_total_average_rho::bs_type (), true))
  , calc_perf_density_ (BS_KERNEL.create_object (wellbore_density_calc::bs_type (), true))
  , calc_perf_bhp_ (BS_KERNEL.create_object (calc_perf_bhp::bs_type (), true))
  {
    clear_data ();
    bhp_ = 0;
  }

  well::well (const std::string &well_name)
  : well_events_init_ (this)
  , i_coord_ (-1)
  , j_coord_ (-1)
  , bhp_depth_ (0)
  , exploitation_factor_ (1)
  , init_approx_is_calc_ (false)
  , input_reference_depth_ (0)
  , calc_well_pressure_ (BS_KERNEL.create_object (calc_well_pressure::bs_type (), true))
  , calc_rho_ (BS_KERNEL.create_object (calc_total_average_rho::bs_type (), true))
  , calc_perf_density_ (BS_KERNEL.create_object (wellbore_density_calc::bs_type (), true))
  , calc_perf_bhp_ (BS_KERNEL.create_object (calc_perf_bhp::bs_type (), true))
  {
    set_name (well_name);

    clear_data ();
    bhp_ = 0;
  }

  well::well (const well &w)
        : bs_refcounter ()
        , well_events_init_ (w.well_events_init_)
  {
    *this = w;
  }

  const std::string &
  well::name () const
  {
    return name_;
  }

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  void
  well::set_coord (index_t i_coord, index_t j_coord)
  {
    i_coord_ = i_coord;
    j_coord_ = j_coord;
  }
  void
  well::set_bhp_depth (item_t bhp_depth)
  {
    bhp_depth_ = bhp_depth;
  }
  void
  well::set_state (well_state_type state, const sp_calc_model_t &calc_model)
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
  void
  well::set_exploitation_factor (item_t exploitation_factor)
  {
    exploitation_factor_ = exploitation_factor;
  }

  void
  well::set_controller (sp_well_controller_t controller)
  {
    well_controller_ = controller;
  }
  void
  well::set_limit_operation (sp_well_limit_operation_t limit_operation)
  {
    well_limit_operation_ = limit_operation;
  }
  well::sp_well_controller_t
  well::get_controller () const
    {
      return well_controller_;
    }
  well::sp_well_limit_operation_t
  well::get_limit_operation () const
    {
      return well_limit_operation_;
    }
#if 0
  well::ijk_filter_t
  well::filter_connections (index_t /*i_coord*/, index_t /*j_coord*/, index_t /*k*/,
                                        index_t /*z1*/, index_t /*z2*/)
  {
    return ijk_filter_t (connection_list_, ijk_pred_t (0));  // TODO:
  }
  well::ijk1k2_filter_t
  well::filter_connections (index_t /*i_coord*/, index_t /*j_coord*/, index_t /*kw1*/,
                                        index_t /*kw2*/)
  {
    return ijk1k2_filter_t (connection_list_, ijk1k2_pred_t (0));  // TODO:
  }
#endif //0
  const std::string &
  well::get_name () const
    {
      return name_;
    }
  void
  well::set_name (const std::string &name)
  {
    name_ = name;
  }

  void
  well::compute_connection_factor (const physical_constants &internal_constants,
      const sp_params_t &params,
      const sp_mesh_iface_t &mesh,
      const stdv_float &perm,
      const stdv_float &ntg,
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

  bool
  well::fi_check_limits () const
    {
      BS_ASSERT (well_limit_operation_);
      if (!well_limit_operation_)
        {
          return false;
          //throw bs_exception ("well::fi_check_limits", "well_limit_operation is null");
        }

      return well_limit_operation_->fi_check_limits ();
    }

  bool
  well::check_connections_bhp (const spv_double & /*pressure*/) const
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

  void
  well::shut_well (const sp_calc_model_t &calc_model)
  {
#ifdef _DEBUG
    BOSOUT (section::wells, level::debug) << boost::format ("[%s]: shut") % name () << bs_end;
#endif

    well_state_.state = well_shut;
    rate_ = 0;

    const t_double *pressure = &(*calc_model->pressure)[0];
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

  bool
  well::is_shut () const
    {
      return well_state_.state == well_shut;
    }

  well::item_t
  well::get_reference_depth (const sp_mesh_iface_t &mesh) const
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

  well::item_t
  well::get_reference_pressure () const
    {
      return bhp_;
    }

  void
  well::set_bhp (item_t bhp)
  {
    bhp_ = bhp;
  }

  bool
  well::is_bhp () const
    {
      BS_ASSERT (well_controller_);
      return well_controller_->is_bhp ();
    }
  bool
  well::is_rate () const
    {
      BS_ASSERT (well_controller_);
      return well_controller_->is_rate ();
    }

  const well::sp_well_controller_t &
  well::get_well_controller () const
    {
      return well_controller_;
    }

  //////////////////////////////////////////////////////////////////////////

  well::item_t
  well::get_input_rate () const
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

  void
  well::process_impl (bool, double, const sp_calc_model_t &, const sp_mesh_iface_t &, BS_SP (jacobian) &)
  {
    bs_throw_exception ("PURE CALL");
  }
  void
  well::process (bool is_start, double dt, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, BS_SP (jacobian) &jacobian)
  {
    if (is_shut ())
      {
        BOSOUT (section::wells, level::debug) << "[" << name () << "] is_shut" << bs_end;
        return ;
      }

#ifdef _DEBUG
    if (!is_shut () && is_no_connections ())
      {
        bs_throw_exception (boost::format ("[%s]: not shut but connection list is empty") % name ());
      }
#endif

    pre_process_facilities ();
    process_impl (is_start, dt, calc_model, mesh, jacobian);
    post_process_facilities ();
  }
  void
  well::clear_data ()
  {
    rate_     = 0;
    rate_rc_  = 0;
    gor_      = 0;
  }

  void
  well::restore_solution (double /*dt*/, const spv_double &/*p_sol*/, const spv_double &/*s_sol*/, index_t /*block_size*/)
  {
    bs_throw_exception ("PURE CALL");
  }

  void
  well::custom_init (const sp_calc_model_t &/*mdl*/)
  {
  }

  void
  well::pre_large_step (const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh)
  {
    compute_connection_factor (calc_model->internal_constants,
      calc_model->ts_params,
      mesh,
      calc_model->rock_grid_prop->permeability,
      calc_model->rock_grid_prop->net_to_gros,
      false);

    check_shut ();
    custom_init (calc_model);
  }
  void
  well::pre_small_step ()
  {
    BS_ASSERT (well_controller_);

    saved_well_state_ = well_state_;
    well_controller_->save_control ();
  }
  void
  well::pre_newton_step ()
  {
    BS_ASSERT (well_controller_);

    saved_niter_well_state_ = well_state_;
    well_controller_->save_niter_control ();
  }
  void
  well::restart_small_step ()
  {
    BS_ASSERT (well_controller_);

    well_state_ = saved_well_state_;
    init_approx_is_calc_ = well_controller_->restore_control ();
  }
  void
  well::restart_newton_step ()
  {
    BS_ASSERT (well_controller_);

    well_state_ = saved_niter_well_state_;
    init_approx_is_calc_ = well_controller_->restore_niter_control ();
  }

  void
  well::pre_process_facilities ()
  {
    for (size_t i = 0, cnt = well_facility_list_.size (); i < cnt; ++i)
      {
        well_facility_list_[i]->pre_process (this);
      }
  }

  void
  well::post_process_facilities ()
  {
    for (size_t i = 0, cnt = well_facility_list_.size (); i < cnt; ++i)
      {
        well_facility_list_[i]->post_process (this);
      }
  }

  void
  well::add_well_facility (const sp_well_facility_t &facility)
  {
    well_facility_list_.push_back (facility);
  }

  void
  well::delete_well_facility (well_facility_list_t::iterator iter)
  {
    well_facility_list_.erase (iter);
  }

  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE (well);
  BLUE_SKY_TYPE_STD_COPY (well);
  BLUE_SKY_TYPE_IMPL (well, facility_base, "well_seq", "well_seq", "well_seq");
  //////////////////////////////////////////////////////////////////////////


  //////////////////////////////////////////////////////////////////////////
  bool
  calc_well_register_types (const blue_sky::plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, wells::connection::bs_type ());   BS_ASSERT (res);
    res &= BS_KERNEL.register_type (pd, well::bs_type ());                BS_ASSERT (res);

    return res;
  }

} // namespace blue_sky
