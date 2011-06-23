/**
 * */
#ifndef BS_EAGLE_DATA_SAVER_HDF5_H_
#define BS_EAGLE_DATA_SAVER_HDF5_H_

#include "bs_hdf5_storage_v2.h"
#include "well_results_storage.h"

namespace blue_sky {

  struct data_saver::impl
  {
    typedef strategy_t::item_t   item_t;

    typedef calc_model calc_model_t;
    typedef rs_mesh_iface mesh_iface_t;
    typedef jacobian_matrix jacobian_matrix_t; 
    typedef facility_manager facility_manager_t;

    typedef facility_manager_t::well_const_iterator_t      well_iterator_t;

    typedef smart_ptr <calc_model_t>      sp_calc_model_t;
    typedef smart_ptr <mesh_iface_t>      sp_mesh_iface_t;
    typedef smart_ptr <jacobian_matrix_t> sp_jmatrix_t;

    impl (const std::string &name)
    : file_ (name)
    {
    }

    template <size_t N>
    std::string concat (boost::array <std::string, N> const &list)
    {
      std::string s = "";
      for (size_t i = 0; i < N; ++i)
        {
          s += list[i] + ",";
        }

      return s;
    }

    void
    write_well_results (const sp_calc_model_t &calc_model,
      well_iterator_t wb,
      const well_iterator_t &we,
      double time)
    {
      bool write_conn_results_to_hdf5 = calc_model->ts_params->get_bool (fi_params::WRITE_CONN_RESULTS_TO_HDF5);
      hdf5_file &file = file_;

      typedef smart_ptr <well> sp_well_t;
      typedef calc_model_t::sp_connection_t                  sp_connection_t;

      boost::array <std::string, 35> d_params_fields = {
        "prod.oil",
        "prod.water",
        "prod.gas",
        "prod.liquid",
        "inje.oil" ,
        "inje.water",
        "inje.gas",
        "_", "_", "_", "_", "_", "_", "_",
        "bhp",
        "bhp_history",
        "total.prod.oil",
        "total.prod.water",
        "total.prod.gas",
        "total.prod.liquid",
        "total.inje.oil",
        "total.inje.water",
        "total.inje.gas",
        "_", "_", "_", "_", "_", "_",
        "_", "_", "_", "_", "_", "_",
      };

      boost::array <std::string, 2> i_params_fields = {
        "_", "status_old",
      };

      file["/meta/wells"]
        .write ("d_params", concat (d_params_fields))
        .write ("i_params", concat (i_params_fields))
        ;

      //TODO: groups
      for (well_iterator_t w = wb; w != we; ++w)
        {
          sp_well_t ws (w->second, bs_dynamic_cast ());

          file["/wells/" + ws->get_name ()]
            .write ("group", std::string ("FIELD"))
            .write ("d_params", 
              hdf5_pod (-ws->rate ().prod.oil)
                << -ws->rate ().prod.water 
                << -ws->rate ().prod.gas
                << -ws->rate ().prod.liquid
                <<  ws->rate ().inj.oil
                <<  ws->rate ().inj.water
                <<  ws->rate ().inj.gas
                << 0 << 0 << 0 << 0 << 0 << 0 << 0
                << ws->get_well_controller ()->bhp ()
                << ws->get_well_controller ()->bhp_history ()
                << -ws->rate_total ().prod.oil
                << -ws->rate_total ().prod.water
                << -ws->rate_total ().prod.gas
                << -ws->rate_total ().prod.liquid
                <<  ws->rate_total ().inj.oil
                <<  ws->rate_total ().inj.water
                <<  ws->rate_total ().inj.gas
                << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << 0)
            .write ("i_params", 
              hdf5_pod (0) << detail::get_status_old (ws))
            ;

          if (write_conn_results_to_hdf5)
            {
              typedef well well_t;
              typename well_t::connection_iterator_t it = ws->connections_begin (), e = ws->connections_end ();
              for (; it != e; ++it)
                {
                  const sp_connection_t &ci (*it);
                  int cell = ci->n_block ();

                  hdf5_group_v2 group = file[(boost::format ("/wells/%s/connections/%d") % ws->get_name () % cell).str ()];
                  hdf5_struct <item_t> d_params;

                  if (!ci->is_shut ())
                    {
                      d_params.write (
                        hdf5_pod (-ci->rate ().prod.oil)
                        << -ci->rate ().prod.water
                        << -ci->rate ().prod.gas
                        << -ci->rate ().prod.liquid
                        <<  ci->rate ().inj.oil
                        <<  ci->rate ().inj.water
                        <<  ci->rate ().inj.gas
                      );
                    }
                  else
                    {
                      d_params.write (
                        hdf5_pod (0) << 0 << 0 << 0 << 0 << 0 << 0
                      );
                    }

                  group.write ("d_params", 
                    d_params << ci->get_fact () 
                      << ci->get_cur_bhp () 
                      << ci->get_bulkp () 
                      << ci->connection_depth);

                  group.write ("i_params", 
                    hdf5_pod ((size_t)ci->get_status ()) 
                      << 0 
                      << ci->i_coord ()
                      << ci->j_coord () 
                      << ci->k_coord () 
                      << ci->k_coord ());
                }
            }
        }
    }

    void
    write_vector (const char *name, 
      const shared_vector <item_t> &array, 
      double time,
      shared_vector <float> &temp)
    {
      const size_t array_size = array.size ();
      if (array_size)
        {
          temp.resize (array_size);

          temp[0]     = array[0];
          item_t min_ = array[0];
          item_t max_ = array[0];
          for (size_t i = 1, cnt = array_size; i < cnt; ++i)
            {
              temp[i] = array[i];
              min_    = (std::min) (min_, array[i]);
              max_    = (std::max) (max_, array[i]);
            }

          file_[name]
            .write ("value", temp)
            .write ("dates", hdf5_pod (time) << max_ << min_);
        }
    }
    void
    write_vector (const char *name,
      const shared_vector <item_t> &array,
      size_t stride, 
      size_t offset,
      double time,
      shared_vector <float> &temp)
    {
      const size_t array_size = array.size ();
      if (array_size)
        {
          temp.resize (array_size);

          if (!stride)
            stride = 1;

          temp[0]     = array[offset];
          item_t min_ = array[offset];
          item_t max_ = array[offset];
          for (size_t i = offset + stride, j = 1; i < array_size; i += stride, ++j)
            {
              temp[j] = static_cast <item_t> (array[i]);
              min_    = (std::min) (min_, array[i]);
              max_    = (std::max) (max_, array[i]);
            }

          file_[name]
            .write ("values", temp)
            .write ("dates", hdf5_pod (time) << max_ << min_);
        }
    }

    void
    write_calc_model_data (const sp_calc_model_t &calc_model,
      const sp_jmatrix_t &jmx,
      size_t large_time_step_num,
      size_t total_time_step_num,
      double time)
    {
      smart_ptr< fi_params> params = calc_model->ts_params;
      size_t save_step_data_period = params->get_int (fi_params::SAVE_STEP_DATA_PERIOD);

      hdf5_file &file = file_;

      if (large_time_step_num % save_step_data_period == 0
        || large_time_step_num == total_time_step_num
        || large_time_step_num == 1)
        {
          file["/results/mainvar"]
            .write ("values", calc_model->main_variable)
            .write ("dates", hdf5_buffer (time));

          shared_vector <float> temp_data;
          if (params->get_bool (fi_params::WRITE_PRESSURE_TO_HDF5))
            {
              write_vector ("/results/pressure", calc_model->pressure, time, temp_data);
            }

          if (params->get_bool (fi_params::WRITE_GAS_OIL_RATIO_TO_HDF5) 
            && calc_model->is_gas () && calc_model->is_oil ())
            {
              write_vector ("/results/gor", calc_model->gas_oil_ratio, time, temp_data);
            }

          if (params->get_bool (fi_params::WRITE_SATURATION_TO_HDF5) && calc_model->n_phases > 1)
            {
              int d_w = calc_model->phase_d[FI_PHASE_WATER];
              int d_g = calc_model->phase_d[FI_PHASE_GAS];

              if (calc_model->is_gas ())
                {
                  write_vector ("/results/sgas",
                    calc_model->saturation_3p, calc_model->n_phases, d_g,
                    time, temp_data);
                }
              if (calc_model->is_water ())
                {
                  write_vector ("/results/swat",
                    calc_model->saturation_3p, calc_model->n_phases, d_w,
                    time, temp_data);
                }
            }

          if (params->get_bool (fi_params::WRITE_CFL_TO_HDF5) 
            && params->get_bool (fi_params::USE_CFL))
            {
              if (jmx->get_cfl_vector ().empty ())
                {
                  jmx->get_cfl_vector ().assign (calc_model->pressure.size (), 0);
                }

              file["/results/cfl"]
                .write ("values", jmx->get_cfl_vector ())
                .write ("dates", hdf5_buffer (time));
            }
        }
    }

    void
    write_mesh (const sp_mesh_iface_t &mesh)
    {
      hdf5_group_v2 mesh_group = file_["/mesh"];
      hdf5_group_v2 mesh_data_group = file_["/mesh/data"];

      mesh->save_info (mesh_group);
      mesh->save_data (mesh_data_group);
      //smart_ptr <rs_smesh_iface > smesh (mesh, bs_dynamic_cast ());
      //if (!smesh)
      //  {
      //    bs_throw_exception (boost::format ("Can't cast mesh to structired mesh, %s") % mesh->bs_resolve_type ().stype_);
      //  }
      //const typename rs_smesh_iface::index_point3d_t &dims = smesh->get_dimens ();

      //hdf5_group_v2 group = file_["/mesh"];

      //group.write ("initial_data", 
      //  hdf5_pod (dims[0]) << dims[1] << dims[2]
      //    << smesh->get_n_active_elements () << smesh->get_n_active_elements ());

      //group.write ("original_elements", smesh->get_ext_to_int ());
      //group.write ("original_planes",   smesh->get_ext_to_int ());
    }

    void
    write_starting_date (const boost::posix_time::ptime &date)
    {
      boost::gregorian::date start_date = date.date ();
      boost::gregorian::date base_date (1900, 1, 1);
      double starting_date = (start_date - base_date).days () + 2;

      file_["/initial_data"].write ("starting_date", hdf5_pod (starting_date));
    }

    template <size_t N>
    shared_vector <strategy_t::item_t>
    convert (calc_model::data_array_t const &data, boost::array <strategy_t::item_t, N> calc_model_data::*x)
    {
      shared_vector <strategy_t::item_t> new_data;
      for (t_long i = 0, cnt = data.size (); i < cnt; ++i)
        {
          calc_model_data const &d = data[i];
          new_data.insert (new_data.end (), (d.*x).begin (), (d.*x).end ());
        }

      return new_data;
    }
    shared_vector <strategy_t::item_t>
    convert (calc_model::data_array_t const &data, strategy_t::item_t calc_model_data::*x)
    {
      shared_vector <strategy_t::item_t> new_data;
      for (t_long i = 0, cnt = data.size (); i < cnt; ++i)
        {
          calc_model_data const &d = data[i];
          new_data.insert (new_data.end (), d.*x);
        }

      return new_data;
    }

    void
    write_calc_model_data (const char *name, 
                           sp_calc_model_t const &cm, 
                           t_long large_step, 
                           t_long small_step, 
                           t_long iteration, 
                           double time)
    {
      char path[1024] = {0};
      sprintf (path, "/data/%s/%ld/%ld/%ld", name, large_step, small_step, iteration);
      printf ("WRITE TO %s\n", path);

      calc_model::data_array_t const &data = cm->data;

      file_[path]
        .write ("cap_pressure",               convert (data, &calc_model_data::cap_pressure))
        .write ("s_deriv_cap_pressure",       convert (data, &calc_model_data::s_deriv_cap_pressure))
        .write ("relative_perm",              convert (data, &calc_model_data::relative_perm))
        .write ("s_deriv_relative_perm",      convert (data, &calc_model_data::s_deriv_relative_perm))
        .write ("p_deriv_gas_oil_ratio",      convert (data, &calc_model_data::p_deriv_gas_oil_ratio))
        .write ("invers_fvf",                 convert (data, &calc_model_data::invers_fvf))
        .write ("p_deriv_invers_fvf",         convert (data, &calc_model_data::p_deriv_invers_fvf))
        .write ("gor_deriv_invers_fvf",       convert (data, &calc_model_data::gor_deriv_invers_fvf))
        .write ("invers_visc",                convert (data, &calc_model_data::invers_viscosity))
        .write ("p_deriv_invers_visc",        convert (data, &calc_model_data::p_deriv_invers_viscosity))
        .write ("gor_deriv_invers_visc",      convert (data, &calc_model_data::gor_deriv_invers_viscosity))
        .write ("invers_visc_fvf",            convert (data, &calc_model_data::invers_visc_fvf))
        .write ("p_deriv_invers_visc_fvf",    convert (data, &calc_model_data::p_deriv_invers_visc_fvf))
        .write ("gor_deriv_invers_visc_fvf",  convert (data, &calc_model_data::gor_deriv_invers_visc_fvf))
        .write ("density",                    convert (data, &calc_model_data::density))
        .write ("p_deriv_density",            convert (data, &calc_model_data::p_deriv_density))
        .write ("gor_deriv_density",          convert (data, &calc_model_data::gor_deriv_density))
        .write ("poro",                       convert (data, &calc_model_data::porosity))
        .write ("p_deriv_poro",               convert (data, &calc_model_data::p_deriv_porosity))
        .write ("truns_mult",                 convert (data, &calc_model_data::truns_mult))
        .write ("p_deriv_truns_mult",         convert (data, &calc_model_data::p_deriv_truns_mult))
        .write ("mob",                        convert (data, &calc_model_data::mobility))
        .write ("p_deriv_mob",                convert (data, &calc_model_data::p_deriv_mobility))
        .write ("s_deriv_mob",                convert (data, &calc_model_data::s_deriv_mobility))
        .write ("prev_fluid_volume",          convert (data, &calc_model_data::prev_fluid_volume))
        ;

      sprintf (path, "/data/%s/rock_grid/%ld/%ld/%ld", name, large_step, small_step, iteration);
      file_[path]
        .write ("ntg", cm->rock_grid_prop->net_to_gros)
        .write ("poro", cm->rock_grid_prop->porosity_p_ref)
        .write ("perm", cm->rock_grid_prop->permeability)
        .write ("volume", cm->rock_grid_prop->volume)
        .write ("multpv", cm->rock_grid_prop->multpv)
        ;
    }

    hdf5_file file_;
  };


} // namespace blue_sky


#endif // #ifndef BS_EAGLE_DATA_SAVER_HDF5_H_

