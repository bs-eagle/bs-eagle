/**
 * \file calc_well.h
 * \brief calculated wells
 * \author Sergey Miryanov
 * \date 23.06.2008
 * */
#ifndef BS_CALC_WELL_H_
#define BS_CALC_WELL_H_

#include "facility_base.h"
#include "well_controller.h"
#include "well_limit_operation.h"
#include "fi_params.h"
#include "well_type_helper.h"

// TODO: BUG
#include "calc_well_pressure.h"
#include "calc_rho.h"
#include "calc_perf_bhp_base.h"
#include "calc_perf_density_base.h"
#include "well_rate_control.h"
#include "array_ext.h"

namespace blue_sky
  {

  ///////////////////////////////////////////////////////////////////////////
  // fwd declarations
  class physical_constants;

  template <typename strategy_t>
  class calc_model;

  namespace wells
    {
    template <typename strategy_t>
    class connection;

    namespace compute_factors
      {
      template <typename strategy_t>
      struct peaceman_model;

      template <typename strategy_t>
      struct baby_odeh_model;
    }
  }
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  enum well_state_type
  {
    well_open,
    well_shut,

    well_state_total,
  };

  well_state_type
  well_state_cast (const std::string &str);
  ///////////////////////////////////////////////////////////////////////////

  template <typename strategy_t>
  struct well_state
    {
      auto_value <well_state_type, well_open>   state;
      auto_value <bool, false>                  is_work;
    };

  template <typename strategy_t>
  class BS_API_PLUGIN well : public facility_base<strategy_t>
    {
    public:

      typedef well <strategy_t>                         well_t;
      typedef typename strategy_t::item_array_t         item_array_t;
      typedef typename strategy_t::rhs_item_array_t     rhs_item_array_t;
      typedef typename strategy_t::index_array_t        index_array_t;
      typedef typename strategy_t::index_t              index_t;
      typedef typename strategy_t::item_t               item_t;
      typedef typename strategy_t::rhs_item_t           rhs_item_t;

      typedef rate_data <strategy_t>                    rate_data_t;
      typedef typename rate_data_t::rate_data_inner     rate_data_inner_t;
      typedef well_state <strategy_t>                   well_state_t;

      typedef wells::well_controller <strategy_t>       well_controller_t;
      typedef wells::well_rate_control <strategy_t>     well_rate_control_t;
      typedef calc_model <strategy_t>                   calc_model_t;
      typedef calc_model_data <strategy_t>              calc_model_data_t;
      typedef wells::connection <strategy_t>            connection_t;
      typedef rs_mesh_iface <strategy_t>                mesh_iface_t;
      typedef jacobian_matrix <strategy_t>              jacobian_matrix_t;

      typedef pvt_oil <strategy_t>                      pvt_oil_t;
      typedef pvt_dead_oil <strategy_t>                 pvt_dead_oil_t;
      typedef pvt_gas <strategy_t>                      pvt_gas_t;
      typedef pvt_water <strategy_t>                    pvt_water_t;

      typedef calc_well_pressure_base <strategy_t>      calc_well_pressure_t;
      typedef calc_rho_base <strategy_t>                calc_rho_base_t;
      typedef calc_perf_density_base <strategy_t>       calc_perf_density_t;
      typedef calc_perf_bhp_base <strategy_t>           calc_perf_bhp_t;

      typedef nc_ptr <calc_model_t>                     nc_calc_model_t;
      typedef nc_ptr <mesh_iface_t>                     nc_mesh_iface_t;
      typedef nc_ptr <jacobian_matrix_t>                nc_jmatrix_t;

      typedef smart_ptr <well_controller_t, true>       sp_well_controller_t;
      typedef smart_ptr <wells::well_limit_operation, true>     sp_well_limit_operation_t;
      typedef smart_ptr <calc_model_t, true>            sp_calc_model_t;
      typedef smart_ptr <mesh_iface_t, true>            sp_mesh_iface_t;
      typedef smart_ptr <well_t, true>                  sp_well_t;

      typedef smart_ptr <pvt_oil_t, true>               sp_pvt_oil_t;
      typedef smart_ptr <pvt_dead_oil_t, true>          sp_pvt_dead_oil_t;
      typedef smart_ptr <pvt_gas_t, true>               sp_pvt_gas_t;
      typedef smart_ptr <pvt_water_t, true>             sp_pvt_water_t;

      typedef smart_ptr <jacobian_matrix_t, true>       sp_jmatrix_t;
      typedef smart_ptr <connection_t, true>            sp_connection_t;

      typedef smart_ptr <calc_well_pressure_t, true>    sp_calc_well_pressure_t;
      typedef smart_ptr <calc_rho_base_t, true>         sp_calc_rho_t;
      typedef smart_ptr <calc_perf_density_t, true>     sp_calc_perf_density_t;
      typedef smart_ptr <calc_perf_bhp_t, true>         sp_calc_perf_bhp_t;

      typedef smart_ptr <fi_params, true>               sp_params_t;

    public:
      virtual sp_connection_t 
      add_connection (index_t i_coord, index_t j_coord, index_t k_coord, index_t n_block)
      {
        bs_throw_exception ("PURE CALL");
      }
      virtual sp_connection_t
      get_connection (index_t idx) const
      {
        bs_throw_exception ("PURE CALL");
      }
      virtual sp_connection_t
      get_connection_map (index_t n_block) const
      {
        bs_throw_exception ("PURE CALL");
      }
      virtual size_t 
      get_connections_count () const
      {
        bs_throw_exception ("PURE CALL");
      }

      void set_coord (index_t i_coord, index_t j_coord);
      void set_bhp_depth (item_t bhp_depth);
      void set_state (well_state_type well_state, const sp_calc_model_t &calc_model);
      void set_exploitation_factor (item_t exploitation_factor);

      void set_controller (sp_well_controller_t controller);
      void set_limit_operation (sp_well_limit_operation_t limit_operation);

      sp_well_controller_t get_controller () const;
      sp_well_limit_operation_t get_limit_operation () const;

      const std::string &get_name () const;
      void set_name (const std::string &name);

      item_t get_bhp_depth () const
        {
          return bhp_depth_;
        }

      item_t bhp () const
        {
          return bhp_;
        }

      const rate_data_t &rate () const
        {
          return rate_;
        }
      const rate_data_t &rate_total () const
        {
          return rate_total_;
        }

      const rate_data_inner_t &
      rate_prod () const
      {
        return rate_.prod;
      }
      const rate_data_inner_t &
      rate_inj () const
      {
        return rate_.inj;
      }

      bool is_open () const
        {
          if (well_state_.state == well_open)
            return true;
          return false;
        }

      const std::string &name () const;

      void set_bhp (item_t bhp);

      bool fi_check_limits () const;

      item_t get_reference_depth (const sp_mesh_iface_t &mesh) const;
      item_t get_reference_pressure () const;

      bool check_connections_bhp (const item_array_t &pressure) const;

      bool is_bhp () const;
      bool is_rate () const;
      bool is_shut () const;

      const sp_well_controller_t &get_well_controller () const;

      item_t get_input_rate () const;

      bool check_shut (const sp_calc_model_t &calc_model);

      void reset_init_approx ();

      virtual array_ext <item_t> get_ww_value ();
      virtual array_ext <item_t> get_bw_value ();

      virtual void eliminate      (rhs_item_t *array, index_t rw_index, index_t wr_index, double dt, index_t block_size) const;
      virtual void process        (bool is_start, double dt, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, sp_jmatrix_t &jmatrix);
      virtual void clear_data     ();

      virtual void fill_rows (index_array_t &rows) const;
      virtual void fill_jacobian (double dt, index_t block_size, const index_array_t &rows, index_array_t &cols, rhs_item_array_t &values, index_array_t &markers) const;
      virtual void fill_rhs (double dt, index_t n_phases, bool is_g, bool is_o, bool is_w, rhs_item_array_t &rhs) const;
      virtual void restore_solution (double dt, const item_array_t &p_sol, const item_array_t &s_sol, index_t block_size);

      virtual void custom_init (const sp_calc_model_t &mdl);

      virtual void pre_large_step (const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh);
      virtual void pre_small_step ();
      virtual void pre_newton_step ();
      virtual void restart_small_step ();
      virtual void restart_newton_step ();

      virtual ~well ();

      well (const std::string &well_name);

    public:
      BLUE_SKY_TYPE_DECL_T (well <strategy_t>);

    protected:

      void shut_well (const sp_calc_model_t &well);
      void reset_rhs_block ();

      void compute_connection_factor (const physical_constants &internal_constants,
                                      const sp_params_t &params,
                                      const sp_mesh_iface_t &mesh,
                                      const item_array_t &perm,
                                      const item_array_t &ntg,
                                      bool ro_calc_flag);

    public:

      std::string                 name_;

      auto_value <index_t, -1>    i_coord_;
      auto_value <index_t, -1>    j_coord_;

      sp_well_controller_t        well_controller_;

//  private:
    public:
      auto_value <item_t, -1>     bhp_depth_;
      auto_value <item_t, 1>      exploitation_factor_;

      well_state_t                well_state_;
      well_state_t                saved_well_state_;
      well_state_t                saved_niter_well_state_;

      sp_well_limit_operation_t   well_limit_operation_;

      rate_data_t                 rate_;
      rate_data_t                 rate_total_;
      rate_data_t                 rate_rc_;

      item_t                      bhp_;
      item_t                      gor_;

    protected:

      auto_value <bool, false>    init_approx_is_calc_;
      auto_value <item_t>         input_reference_depth_;

      sp_calc_well_pressure_t     calc_well_pressure_;
      sp_calc_rho_t               calc_rho_;
      sp_calc_perf_density_t      calc_perf_density_;
      sp_calc_perf_bhp_t          calc_perf_bhp_;

    protected:
      index_array_t               open_connections_;
    };

  template <typename strategy_t>
  class BS_API_PLUGIN well_factory : public objbase
    {
    public:

      typedef well <strategy_t>                 well_t;
      typedef wells::connection <strategy_t>    connection_t;
      typedef smart_ptr <well_t>                sp_well_t;
      typedef smart_ptr <connection_t, true>    sp_connection_t;

    public:

      virtual ~well_factory () {}

      virtual sp_well_t             create_well (const std::string &group_name, const std::string &well_name) const;
      virtual sp_connection_t create_connection () const;

      BLUE_SKY_TYPE_DECL_T (well_factory <strategy_t>);
    };

  bool
  calc_well_register_types (const blue_sky::plugin_descriptor &pd);

  bool
  well_factory_register_type (const blue_sky::plugin_descriptor &pd);

} // namespace blue_sky


#endif  // #ifndef BS_CALC_WELL_H_
