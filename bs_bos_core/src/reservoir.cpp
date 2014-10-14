/**
 *       \file  reservoir.cpp
 *      \brief  implementation of reservoir class
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  16.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "reservoir.h"
#include "well_connection.h"
#include "facility_manager.h"
#include "jacobian.h"
#include "calc_model.h"
#include "event_filter.h"

#include "closure.h"
#include "for_each_well.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "scal_3p_iface.hpp"
#include BS_STOP_PLUGIN_IMPORT ()

#include "data_saver.h"

namespace blue_sky
  {

  reservoir::~reservoir ()
  {

  }

  //! constructor
  /**
   * \brief  copy-ctor for reservoir
   * \param  cl reservoir instance to be copied
   * */
  reservoir::reservoir(const reservoir & cl)
  : bs_refcounter (cl), objbase (cl)
  {
    if (&cl != this)
      *this = cl;
  }


  /**
   * \brief  'default' reservoir ctor
   * \param  param additional ctor params
   * */
  reservoir::reservoir (bs_type_ctor_param /*param*/)
  {
    facility_list_                = BS_KERNEL.create_object(facility_manager_t::bs_type(), true);
    well_factory_                 = BS_KERNEL.create_object(well_factory_t::bs_type(), true);
    well_controller_factory_      = BS_KERNEL.create_object(controller_factory_t::bs_type(), true);
    well_limit_operation_factory_ = BS_KERNEL.create_object(limit_operation_factory_t::bs_type(), true);
    event_filter_                 = BS_KERNEL.create_object (event_filter::bs_type (), true);

    data_saver_.reset (new data_saver_t);
  }

  reservoir::sp_well_t
  reservoir::get_well (const std::string &group_name, const std::string &well_name) const
    {
      return facility_list_->get_well (group_name, well_name);
    }

  reservoir::sp_well_t
  reservoir::get_well (const std::string &well_name) const
  {
    return facility_list_->get_well (well_name);
  }


  void
  reservoir::add_well (BS_SP (well) w)
  {
    facility_list_->add_well (w);
  }


  reservoir::sp_well_controller_t
  reservoir::create_well_controller (const sp_well_t &owner_well)
  {
    BS_ASSERT (owner_well);

    sp_well_controller_t wc = well_controller_factory_->create_controller ();
    BS_ASSERT (wc) (owner_well->get_name ());

    if (wc)
      owner_well->set_controller (wc);

    return wc;
  }

  reservoir::sp_well_limit_operation_t
  reservoir::create_well_limit_operation (const sp_well_t &owner_well, wells::limit_operation_type operation)
  {
    BS_ASSERT (owner_well);

    sp_well_limit_operation_t wlo = well_limit_operation_factory_->create_limit (operation);
    BS_ASSERT (wlo) (operation) (owner_well->get_name ());

    if (wlo)
      owner_well->set_limit_operation (wlo);

    return wlo;
  }


  reservoir::sp_connection_t
  reservoir::create_connection ()
  {
    sp_connection_t con = well_factory_->create_connection ();

    if (!con)
      throw bs_exception ("reservoir::create_connection", "can't create connection");

    return con;
  }

  void
  reservoir::save_data (const sp_storage_t &storage) const
    {
      BS_ASSERT (facility_list_);
      facility_list_->save_data (storage);
    }

  reservoir::item_t
  reservoir::pressure () const
    {
      BS_ASSERT (false && "NOT IMPL YET");
      return 0;
    }

  bool
  reservoir::check_limits (const sp_params_t &/*params*/) const
    {
      BS_ASSERT (false && "NOT IMPL YET");
      return false;
    }

  void
  reservoir::pre_small_step ()
  {
    for_each_facility (*facility_list_, closure <void, facility_t> (&facility_t::pre_small_step));
  }
  void
  reservoir::restart_small_step ()
  {
    for_each_facility (*facility_list_, closure <void, facility_t> (&facility_t::restart_small_step));
  }
  void
  reservoir::pre_newton_step ()
  {
    for_each_facility (*facility_list_, closure <void, facility_t> (&facility_t::pre_newton_step));
  }
  void
  reservoir::restart_newton_step ()
  {
    for_each_facility (*facility_list_, closure <void, facility_t> (&facility_t::restart_newton_step));
  }

  namespace detail
  {
    void
    init_rows (BS_SP (facility_manager) &wells, spv_long &rows, t_long cells)
    {
      rows->init (cells + 1);

      facility_manager::well_const_iterator_t wit = wells->wells_begin ();
      facility_manager::well_const_iterator_t we = wells->wells_end ();
      for (; wit != we; ++wit)
        {
          wit->second->fill_rows (rows);
        }

      t_long *r = rows->data ();
      for (size_t i = 0, cnt = rows->size () - 1; i < cnt; ++i)
        {
          r[i + 1] += r[i];
        }
    }
  }

  // FIXME: remove, obsolete
  void
  reservoir::init_jacobian (const BS_SP (jacobian) & /*jacobian*/, index_t /*n_cells*/)
  {
  }

  void
  reservoir::end_jacobian (BS_SP (jacobian) &jacobian, t_double dt, t_long block_size, t_long cells)
  {
    BS_SP (bcsr_matrix_iface) mx = jacobian->get_matrix ("facility");

    spv_long rows = mx->get_rows_ptr ();
    spv_long cols = mx->get_cols_ind ();
    spv_float values = mx->get_values ();

    blue_sky::detail::init_rows (facility_list_, rows, cells);
    if (rows->empty ())
      return ;

    t_long cols_count = rows->back ();
    if (cols_count == 0)
      {
        bs_throw_exception ("cols count is 0");
      }

    cols->init (cols_count);
	cols->assign(-1);
    values->init (cols_count * block_size * block_size, 0.0);
    markers_.assign (rows->size (), 0);

    mx->set_n_cols ((t_long)cols->size ());

    // FIXME: WTF?? n_rows
    //mx->n_rows = rows->size () - 1;

    facility_manager::well_const_iterator_t wit = facility_list_->wells_begin ();
    facility_manager::well_const_iterator_t we = facility_list_->wells_end ();
    for (; wit != we; ++wit)
      {
        wit->second->fill_jacobian (dt, block_size, rows, cols, values, markers_);
      }

#ifdef _DEBUG
    for (size_t i = 0, cnt = cols->size (); i < cnt; ++i)
      {
        if ((*cols)[i] == -1)
          {
            bs_throw_exception ("Cols i is negative");
          }
      }
#endif
  }

  void
  reservoir::calc_wells (int istart, double dt, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, BS_SP (jacobian) &jacobian)
  {
    facility_manager::well_const_iterator_t wit = facility_list_->wells_begin ();
    facility_manager::well_const_iterator_t we = facility_list_->wells_end ();
    for (; wit != we; ++wit)
      {
        wit->second->process (istart, dt, calc_model, mesh, jacobian);
      }
  }

  void
  reservoir::fill_rhs_wells (double dt, const sp_calc_model_t &calc_model, rhs_item_array_t &rhs, 
                             bool /*update_after_gauss_elimination*/) const
  {
    index_t n_phases = calc_model->n_phases;
    bool is_g = calc_model->is_gas ();
    bool is_o = calc_model->is_oil ();
    bool is_w = calc_model->is_water ();

    facility_manager::well_const_iterator_t wit = facility_list_->wells_begin ();
    facility_manager::well_const_iterator_t we = facility_list_->wells_end ();
    for (; wit != we; ++wit)
      {
        wit->second->fill_rhs (dt, n_phases, is_g, is_o, is_w, rhs);
      }
  }

  reservoir::sp_facility_manager_t
  reservoir::get_facility_list () const
  {
    return facility_list_;
  }

  void
  reservoir::restore_wells_solution (double dt, const spv_double &p_sol, const spv_double &s_sol, index_t block_size)
  {
    facility_manager::well_const_iterator_t wit = facility_list_->wells_begin ();
    facility_manager::well_const_iterator_t we = facility_list_->wells_end ();
    for (; wit != we; ++wit)
      {
        wit->second->restore_solution (dt, p_sol, s_sol, block_size);
      }
  }

  void
  reservoir::pre_large_step (const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh)
  {
    facility_manager_t::well_const_iterator_t wb = facility_list_->wells_begin ();
    facility_manager_t::well_const_iterator_t we = facility_list_->wells_end ();
    for (; wb != we; ++wb)
      {
        wb->second->pre_large_step (calc_model, mesh);
      }
    //for_each_facility (*facility_list_, closure <void, facility_t, const calc_model_t &, const sp_mesh_iface_t &> (&facility_t::pre_large_step, calc_model, mesh));
  }


  void
  reservoir::add_filter_well (const std::string &well_name)
  {
    event_filter_->add_filter_well (well_name);
  }
  const reservoir::sp_event_filter_t &
  reservoir::get_event_filter () const
  {
    return event_filter_;
  }

  void
  reservoir::write_step_to_storage (const sp_calc_model_t &calc_model, 
    const sp_mesh_iface_t & /*mesh*/, 
    const BS_SP (jacobian) &jacobian, 
    size_t large_time_step_num, 
    size_t total_time_step_num, 
    double time)
  {
    BS_ASSERT (data_saver_);

    data_saver_->write_well_results (calc_model, facility_list_->wells_begin (), facility_list_->wells_end (), time);
    data_saver_->write_fip_results  (calc_model);
    data_saver_->write_calc_model_data (calc_model, jacobian, large_time_step_num, total_time_step_num, time);
  }
  void
  reservoir::write_mesh_to_storage (const sp_mesh_iface_t &mesh)
  {
    BS_ASSERT (data_saver_);
    data_saver_->write_mesh (mesh);
  }
  void
  reservoir::write_starting_date_to_storage (const boost::posix_time::ptime &date)
  {
    BS_ASSERT (data_saver_);
    data_saver_->write_starting_date (date);
  }

  void
  reservoir::open_storage (const std::string &name)
  {
    BS_ASSERT (data_saver_) (name);
    data_saver_->open_storage (name);
  }

  reservoir::sp_well_factory_t
  reservoir::get_well_factory () const
  {
    return well_factory_;
  }
  reservoir::sp_well_controller_factory_t
  reservoir::get_well_controller_factory () const
  {
    return well_controller_factory_;
  }
  reservoir::sp_well_limit_operation_factory_t
  reservoir::get_well_limit_operation_factory () const
  {
    return well_limit_operation_factory_;
  }

  void
  reservoir::set_well_factory (const sp_well_factory_t &factory)
  {
    well_factory_ = factory;
  }
  void
  reservoir::set_well_controller_factory (const sp_well_controller_factory_t &factory)
  {
    well_controller_factory_ = factory;
  }
  void
  reservoir::set_well_limit_operation_factory (const sp_well_limit_operation_factory_t &factory)
  {
    well_limit_operation_factory_ = factory;
  }

  //bs stuff
  BLUE_SKY_TYPE_STD_CREATE (reservoir)
  BLUE_SKY_TYPE_STD_COPY (reservoir)
  BLUE_SKY_TYPE_IMPL (reservoir, objbase, "reservior_seq", "reservior_seq", "reservior_seq")

} // namespace blue_sky

