#ifndef BS_RESERVOIR_H_
#define BS_RESERVOIR_H_

#include "calc_well.h"
#include "event_filter.h"
#include "jacobian.h"

#ifdef _HDF5
#include "bs_hdf5_storage.h"
#include "rs_smesh_iface.h"
#endif

namespace blue_sky
  {

  template <typename strategy_t>
  class calc_model;

  template <typename strategy_t>
  class reservoir_simulator;

  template <typename strategy_t>
  class facility_manager;

  template <typename strategy_t>
  class well;

  template <typename strategy_t>
  class well_factory;

  class data_storage_interface;

  class fi_params;

  namespace wells
    {
    template <typename strategy_t>
    class well_controller;

    template <typename strategy_t>
    class well_rate_control;

    template <typename strategy_t>
    class connection;

    template <typename strategy_t>
    class well_controller_factory;

    class well_limit_operation_factory;

    class well_limit_operation;
  }

  template <typename strategy_t>
  class BS_API_PLUGIN reservoir : public objbase
    {
    public:

      typedef typename strategy_t::item_t                     item_t;
      typedef typename strategy_t::index_t                    index_t;
      typedef typename strategy_t::item_array_t               item_array_t;
      typedef typename strategy_t::rhs_item_array_t           rhs_item_array_t;
      typedef typename strategy_t::index_array_t              index_array_t;

      typedef facility_base <strategy_t>                      facility_t;
      typedef well <strategy_t>                               well_t;
      typedef wells::well_controller <strategy_t>             well_controller_t;
      typedef wells::well_limit_operation                     well_limit_operation_t;
      typedef wells::well_rate_control <strategy_t>           well_rate_control_t;
      typedef wells::connection <strategy_t>                  connection_t;
      typedef facility_manager <strategy_t>                   facility_manager_t;
      typedef wells::well_controller_factory <strategy_t>     controller_factory_t;
      typedef wells::well_limit_operation_factory             limit_operation_factory_t;
      typedef well_factory <strategy_t>                       well_factory_t;
      typedef calc_model <strategy_t>                         calc_model_t;
      typedef rs_mesh_iface <strategy_t>                      mesh_iface_t;
      typedef jacobian <strategy_t>                           jacobian_t;
      typedef jacobian_matrix <strategy_t>                    jacobian_matrix_t;

      typedef reservoir_simulator <strategy_t>                reservoir_simulator_t;
      typedef rate_data <strategy_t>                          rate_data_t;
      typedef typename rate_data_t::rate_data_inner           rate_data_inner_t;

      typedef smart_ptr <reservoir_simulator_t, true>         sp_top_t;

      typedef smart_ptr <well_t, true>                        sp_well_t;
      typedef smart_ptr <well_controller_t, true>        			sp_well_controller_t;
      typedef smart_ptr <well_limit_operation_t, true>   			sp_well_limit_operation_t;
      typedef smart_ptr <well_rate_control_t, true>						sp_rate_control_t;
      typedef smart_ptr <connection_t, true>              		sp_connection_t;
      typedef smart_ptr <mesh_iface_t, true>                  sp_mesh_iface_t;
      typedef smart_ptr <calc_model_t, true>                  sp_calc_model_t;
      typedef smart_ptr <jacobian_t, true>                    sp_jacobian_t;
      typedef smart_ptr <jacobian_matrix_t, true>             sp_jacobian_matrix_t;
      typedef smart_ptr <jacobian_matrix_t, true>             sp_jmatrix_t;

      typedef smart_ptr <facility_manager_t, true>						sp_facility_manager_t;
      typedef smart_ptr <data_storage_interface, true>        sp_storage_t;

      typedef smart_ptr <controller_factory_t, true>			    sp_well_controller_factory_t;
      typedef smart_ptr <limit_operation_factory_t, true>	    sp_well_limit_operation_factory_t;
      typedef smart_ptr <well_factory_t, true>							  sp_well_factory_t;

      typedef smart_ptr <fi_params, true>                     sp_params_t;
      typedef smart_ptr <event_filter, true>                  sp_event_filter_t;

#ifdef _HDF5
      typedef smart_ptr <bs_hdf5_storage, true>               sp_bs_hdf5_storage;
#endif

    public:

      BLUE_SKY_TYPE_DECL_T (reservoir <strategy_t>);
      ~reservoir ();

    public:

      sp_well_t                 get_well (const std::string &group_name, const std::string &well_name) const;
      sp_well_t                 get_well (const std::string &well_name) const;

      sp_well_t                 create_well (const std::string &group_name, const std::string &well_name);
      sp_well_controller_t      create_well_controller (const sp_well_t &owner_well);
      sp_well_limit_operation_t create_well_limit_operation (const sp_well_t &owner_well, wells::limit_operation_type operation);
      sp_connection_t           create_connection ();

      sp_rate_control_t					create_bhp_control (bool is_prod, const sp_calc_model_t &calc_model);
      sp_rate_control_t					create_rate_control (bool is_prod, wells::rate_control_type control_type, const sp_calc_model_t &calc_model);

      void											save_data (const sp_storage_t &storage) const;

      item_t										pressure () const;

      bool                      check_limits (const sp_params_t &params) const;

      void                      pre_large_step (const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh);
      void                      pre_small_step ();
      void                      pre_newton_step ();
      void                      restart_small_step ();
      void                      restart_newton_step ();

      void                      init_jacobian (const sp_jmatrix_t &jmx, index_t n_cells);
      void                      end_jacobian (item_t dt, const sp_calc_model_t &calc_model, sp_jacobian_t &jacobian);

      void                      restore_wells_solution (double dt, const item_array_t &p_sol, const item_array_t &s_sol, index_t block_size);

      void                      calc_wells (int istart, double dt, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, sp_jmatrix_t &jmatrix);
      void                      fill_rhs_wells (double dt, const sp_calc_model_t &calc_model, rhs_item_array_t &rhs, bool update_after_gauss_elimination) const;

      sp_facility_manager_t			get_facility_list () const;
      size_t                    get_connections_count () const;

      void                      add_filter_well (const std::string &well_name);
      const sp_event_filter_t   &get_event_filter () const;

      sp_well_factory_t                 get_well_factory () const;
      sp_well_controller_factory_t      get_well_controller_factory () const;
      sp_well_limit_operation_factory_t get_well_limit_operation_factory () const;

      void                      set_well_factory (const sp_well_factory_t &factory);
      void                      set_well_controller_factory (const sp_well_controller_factory_t &factory);
      void                      set_well_limit_operation_factory (const sp_well_limit_operation_factory_t &factory);

#ifdef _HDF5
      void                      open_hdf5_file (const std::string &filename) const;
      void                      close_hdf5_file () const;
      void                      write_step_to_hdf5 (const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, const sp_jmatrix_t &jmx, int, int, item_t time) const;
      void                      write_mesh_to_hdf5 (const smart_ptr <rs_mesh_iface<strategy_t>, true> &mesh) const;
      const smart_ptr<bs_hdf5_storage, true> get_hdf5_file () const {return hdf5;}
#endif

      const rate_data_t &
      rate () const
      {
        return rate_;
      }

    private:
      void                      init_rows (index_array_t &rows) const;

    public:

      rate_data_t                         rate_;
      rate_data_t                         rate_rc_;
      rate_data_t                         rate_wefac_;
      rate_data_t                         rate_rc_wefac_;
      rate_data_t                         rate_initial_;
      rate_data_t                         rate_total_;

    private:

      sp_facility_manager_t								facility_list_;											//!< manager for facilities (wells, unnamed wells and other)
      sp_well_factory_t										well_factory_;											//!< factory for wells and connections
      sp_well_controller_factory_t				well_controller_factory_;						//!< factory for well controllers
      sp_well_limit_operation_factory_t		well_limit_operation_factory_;			//!< factory for well limit operations

      sp_event_filter_t                   event_filter_;

      index_array_t                       markers_;

#ifdef _HDF5
      sp_bs_hdf5_storage                  hdf5;
#endif
    };


}	// namespace blue_sky

#endif	// #ifndef BS_RESERVOIR_H_

