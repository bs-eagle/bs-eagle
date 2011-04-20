/**
 *       \file  reservoir.cpp
 *      \brief  implementation of reservoir class
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  16.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"

#include "reservoir.h"
#include "well_connection.h"
#include "facility_manager.h"
#include "jacobian.h"
#include "calc_model.h"
#include "event_filter.h"

#ifdef _HDF5
#include "bs_hdf5_storage.h"
#include "results_hdf5_writer.h"
#endif

#include "closure.h"
#include "for_each_well.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "scal_3p.h"
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
#ifdef _HDF5
    hdf5                          = BS_KERNEL.create_object (bs_hdf5_storage::bs_type (), true);
#endif

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


  reservoir::sp_well_t
  reservoir::create_well (const std::string &group_name, const std::string &well_name)
  {
    sp_well_t w = well_factory_->create_well (group_name, well_name);
    BS_ASSERT (w) (group_name) (well_name);

    if (w)
      {
        facility_list_->add_well (w);
      }

    return w;
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

  void
  reservoir::init_rows (index_array_t &rows_) const
    {
      for_each_facility (*facility_list_, closure <void, facility_t, index_array_t &> (&facility_t::fill_rows, rows_));

      // correct rows values
      v_long &rows = *rows_;
      for (index_t i = 0, cnt = (index_t)rows.size () - 1; i < cnt; ++i)
        {
          rows[i + 1] += rows[i];
        }
    }

  void
  reservoir::init_jacobian (const sp_jmatrix_t &jmx, index_t n_cells)
  {
    spv_long rows = jmx->get_facility_matrix ()->get_rows_ptr ();

    // FIXME: dereference
    rows->init (n_cells + 1, 0);
    init_rows (rows);
  }

  void
  reservoir::end_jacobian (item_t dt, const sp_calc_model_t &calc_model, sp_jacobian_t &jacobian)
  {
    const sp_jacobian_matrix_t &jmx (jacobian->get_jmatrix ());

    const spv_long &rows    = jmx->get_facility_matrix ()->get_rows_ptr ();
    spv_long cols           = jmx->get_facility_matrix ()->get_cols_ind ();
    spv_double values       = jmx->get_facility_matrix ()->get_values ();

    if (rows->empty ())
      return ;

    index_t block_size = calc_model->n_phases;
    //index_t b_sqr = block_size * block_size;
    //cols->init (rows->back (), -1);
    //values->init (rows->back () * b_sqr, 0);
    //markers_.assign (rows->size (), 0);

    //jmx->get_facility_matrix ()->n_block_size = calc_model->n_phases;
    //jmx->get_facility_matrix ()->n_rows = (index_t) rows->size () - 1;
    //jmx->get_facility_matrix ()->n_cols = (index_t) cols->size ();

    // FIXME:
    t_long cols_count = rows->back ();
    jmx->get_facility_matrix ()->init (rows->size () - 1, cols_count, calc_model->n_phases, 1);

    for_each_facility (*facility_list_, closure <void, facility_t, double, index_t, const spv_long &, spv_long &, spv_double &, stdv_long &> (&facility_t::fill_jacobian,
        dt, block_size, rows, cols, values, markers_));

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
  reservoir::calc_wells (int istart, double dt, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, sp_jmatrix_t &jmatrix)
  {
    for_each_facility (*facility_list_, closure <void, facility_t, bool, double, const sp_calc_model_t &, const sp_mesh_iface_t &, sp_jmatrix_t &> (&facility_t::process,
      istart, dt, calc_model, mesh, jmatrix));
  }

  void
  reservoir::fill_rhs_wells (double dt, const sp_calc_model_t &calc_model, rhs_item_array_t &rhs, bool update_after_gauss_elimination) const
  {
    for_each_facility (*facility_list_, closure <void, facility_t, double, index_t, bool, bool, bool, rhs_item_array_t &> (&facility_t::fill_rhs,
      dt, calc_model->n_phases, calc_model->is_gas (), calc_model->is_oil (), calc_model->is_water (), rhs));
  }

  reservoir::sp_facility_manager_t
  reservoir::get_facility_list () const
  {
    return facility_list_;
  }

  void
  reservoir::restore_wells_solution (double dt, const spv_double &p_sol, const spv_double &s_sol, index_t block_size)
  {
    for_each_facility (*facility_list_, closure <void, facility_t, double, const spv_double &, const spv_double &, index_t> (&facility_t::restore_solution, dt, p_sol, s_sol, block_size));
  }

  void
  reservoir::pre_large_step (const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh)
  {
    typename facility_manager_t::well_const_iterator_t wb = facility_list_->wells_begin ();
    typename facility_manager_t::well_const_iterator_t we = facility_list_->wells_end ();
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

#ifdef _HDF5
  void
  reservoir::write_mesh_to_hdf5 (const smart_ptr <rs_mesh_iface, true> &mesh) const
    {
      smart_ptr <rs_smesh_iface, true> smesh (mesh, bs_dynamic_cast ());
      hdf5::write_mesh_to_hdf5 (hdf5, smesh);
    }

  void
  reservoir::open_hdf5_file (const std::string &filename) const
    {
      hdf5->lock ()->open (filename.c_str ());
    }

  void
  reservoir::close_hdf5_file () const
    {
      hdf5->lock ()->close ();
    }

  void
  reservoir::write_step_to_hdf5 (const sp_calc_model_t &calc_model,
                           const sp_mesh_iface_t &mesh,
                           const sp_jmatrix_t &jmx, 
                           int l_time_step_num, int total_ts_count,
                           item_t time) const
    {
#ifdef _MPI
  int n_elem_global = mesh->mpi_decomp->mpi_global_n_active_elmts;
  int n_elem_local = mesh->mpi_decomp->get_n_local_own (); //local number of active elments without ghost cells
  int mpi_offset = mesh->mpi_decomp->get_elem_own_displ (); //offset in file
  int mpi_start = mesh->mpi_decomp->get_recv_l (); //offset in array-starting with first own element
#else
  int n_elem_global = mesh->get_n_active_elements ();
  int n_elem_local = n_elem_global;
  int mpi_offset = 0;
  int mpi_start = 0;
#endif

  hdf5::write_well_results (*hdf5, calc_model->well_res, calc_model->ts_params->get_bool (fi_params::WRITE_CONN_RESULTS_TO_HDF5));
  hdf5::write_fip_results (*hdf5, calc_model->fip_res);

  int n_phases = calc_model->n_phases;
  // temporary storage for values, need to extract saturations
  // and converting to float if need.
  shared_vector <float> tmp_buf (n_elem_local);
  shared_vector <int> main_var_tmp (n_elem_local);
  //tmp_buf.resize (n_elem_local);

  if ((l_time_step_num % calc_model->ts_params->get_int (fi_params::SAVE_STEP_DATA_PERIOD) == 0)
      || l_time_step_num == total_ts_count
      || l_time_step_num == 1)
    {
      if (calc_model->ts_params->get_bool (fi_params::WRITE_PRESSURE_TO_HDF5))
        {
          copy_to (&calc_model->pressure[mpi_start], &tmp_buf[0], n_elem_local);
          hdf5->add_to_results ("/results/pressure", &tmp_buf[0], n_elem_global, n_elem_local, mpi_offset, time);
        }

      if (calc_model->ts_params->get_bool (fi_params::WRITE_GAS_OIL_RATIO_TO_HDF5))
        {
          if (calc_model->is_gas () && calc_model->is_oil ())
            {
              copy_to (&calc_model->gas_oil_ratio[mpi_start], &tmp_buf[0], n_elem_local);
              hdf5->add_to_results ("/results/gor", &tmp_buf[0], n_elem_global, n_elem_local, mpi_offset, time);
            }
        }

      copy_to (&calc_model->main_variable[mpi_start], &main_var_tmp[0], n_elem_local);
      hdf5->add_to_results ("/results/mainvar", &main_var_tmp[0], n_elem_global, n_elem_local, mpi_offset, time);
      
      if (calc_model->ts_params->get_bool (fi_params::WRITE_SATURATION_TO_HDF5))
        {
          // for 1-phase no saturations
          if (calc_model->n_phases > 1)
            {
              int d_w = calc_model->phase_d[FI_PHASE_WATER];
              int d_g = calc_model->phase_d[FI_PHASE_GAS];

              if (calc_model->is_gas ())
                {
                  copy_to (&calc_model->saturation_3p[n_phases * mpi_start], &tmp_buf[0], n_phases * n_elem_local, n_phases/*stride*/, d_g/*start*/);
                  hdf5->add_to_results ("/results/sgas", &tmp_buf[0], n_elem_global, n_elem_local, mpi_offset, time);
                }
              if (calc_model->is_water ())
                {
                  copy_to (&calc_model->saturation_3p[n_phases * mpi_start], &tmp_buf[0], n_phases * n_elem_local, n_phases /*stride*/, d_w/*start*/);
                  hdf5->add_to_results ("/results/swat", &tmp_buf[0], n_elem_global, n_elem_local, mpi_offset, time);
                }
            }
        }//saturation

#ifdef BS_BOS_CORE_USE_CSR_ILU_CFL_PREC
      if (calc_model->ts_params->get_bool (fi_params::WRITE_CFL_TO_HDF5) && 
          calc_model->ts_params->get_bool (fi_params::USE_CFL))
      {
        if (jmx->get_cfl_vector ().empty ())
          {
            jmx->get_cfl_vector ().assign (n_elem_local, 0);
          }

        copy_to (&jmx->get_cfl_vector ()[0], &tmp_buf[0], n_elem_local);
        hdf5->add_to_results ("/results/cfl", &tmp_buf[0], n_elem_global, n_elem_local, mpi_offset, time);
      }
#endif

#if 0 //TODO: mesh->n_planes
       //if (model->ts_params.get_b_param (FI_PARAMS_B_WRITE_PLANE_FLOW_RATES_TO_HDF5))
        {
          tmp_buf.resize (mesh->n_planes);

          // parse plane_flow_rates, separating data for each direction of plane and phase
          // data layout: {X_oil X_water X_gas Y_oil Y_water Y_gas Z_oil Z_water Z_gas}, so stride is 3*model->n_phases.
          int d_w = calc_model->phase_d[FI_PHASE_WATER];
          int d_g = calc_model->phase_d[FI_PHASE_GAS];
          int d_o = calc_model->phase_d[FI_PHASE_OIL];

          //create group. open-close else
          hid_t group_planes = hdf5->begin_object (hdf5->get_file_id (), "/results/plane_flow_rate");
          H5Gclose (group_planes);

          if (FI_CHK_OIL (calc_model->phases))
            {
              copy_to_float (calc_model->full_step_plane_flow_rate + d_o * mesh->n_planes, &tmp_buf[0], mesh->n_planes);
              hdf5->add_to_results ("/results/plane_flow_rate/oil", &tmp_buf[0], mesh->n_planes, mesh->n_planes, 0, time);
            }
          if (FI_CHK_GAS (calc_model->phases))
            {
              copy_to_float (calc_model->full_step_plane_flow_rate + d_g * mesh->n_planes, &tmp_buf[0], mesh->n_planes);
              hdf5->add_to_results ("/results/plane_flow_rate/gas", &tmp_buf[0], mesh->n_planes, mesh->n_planes, 0, time);
            }
          if (FI_CHK_WATER (calc_model->phases))
            {
              copy_to_float (calc_model->full_step_plane_flow_rate + d_w * mesh->n_planes, &tmp_buf[0], mesh->n_planes);
              hdf5->add_to_results ("/results/plane_flow_rate/wat", &tmp_buf[0], mesh->n_planes, mesh->n_planes, 0, time);
            }
        }
#endif
    }
    }
#endif // #ifdef HDF5

  void
  reservoir::write_step_to_storage (const sp_calc_model_t &calc_model, 
    const sp_mesh_iface_t &mesh, 
    const sp_jmatrix_t &jmx, 
    size_t large_time_step_num, 
    size_t total_time_step_num, 
    double time)
  {
    BS_ASSERT (data_saver_);

    data_saver_->write_well_results (calc_model, facility_list_->wells_begin (), facility_list_->wells_end (), time);
    data_saver_->write_fip_results  (calc_model);
    data_saver_->write_calc_model_data (calc_model, jmx, large_time_step_num, total_time_step_num, time);
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

