/**
 *       \file  calc_well.h
 *      \brief  Base class for wells
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  23.06.2008
 *  \copyright  This source code is released under the terms of
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_CALC_WELL_H_
#define BS_CALC_WELL_H_

#include "facility_base.h"
#include "well_controller.h"
#include "well_limit_operation.h"
#include "fi_params.h"
#include "well_type_helper.h"
#include "connection_iterator.h"
#include "well_facility.h"
#include "shared_vector.h"

// TODO: BUG
#include "calc_well_pressure.h"
#include "calc_rho.h"
#include "calc_perf_bhp_base.h"
#include "calc_perf_density_base.h"

#include "kernel_signals.h"

namespace blue_sky
  {

  ///////////////////////////////////////////////////////////////////////////
  // fwd declarations
  class physical_constants;
  class calc_model;

  namespace wells
    {
    class connection;

    namespace compute_factors
      {
      struct peaceman_model;
      struct baby_odeh_model;
    }
  }
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  /**
   * \enum  well_state_type
   * \brief Type of well status
   * */
  enum well_state_type
  {
    well_open,          //!< Well is open
    well_shut,          //!< Well is shut

    well_state_total,
  };

  well_state_type
  well_state_cast (const std::string &str);
  ///////////////////////////////////////////////////////////////////////////

  /**
   * \class well_state
   * \brief Incapsulates state of well
   * */
  struct well_state
    {
      auto_value <well_state_type, well_open>   state;    //!< Status of well directed by model
      auto_value <bool, false>                  is_work;  //!< Is well work by inner state
    };

  /**
   * \class well
   * \brief Base class for wells
   * */
  class BS_API_PLUGIN well : public facility_base
    {
    public:

      typedef facility_base base_t;
      typedef well well_t;
      typedef strategy_t::item_array_t         item_array_t;
      typedef strategy_t::rhs_item_array_t     rhs_item_array_t;
      typedef strategy_t::index_array_t        index_array_t;
      typedef strategy_t::index_t              index_t;
      typedef strategy_t::item_t               item_t;
      typedef strategy_t::rhs_item_t           rhs_item_t;

      typedef rate_data rate_data_t;
      typedef rate_data_t::rate_data_inner     rate_data_inner_t;
      typedef well_state well_state_t;

      typedef wells::well_controller well_controller_t;
      typedef wells::well_rate_control well_rate_control_t;
      typedef calc_model calc_model_t;
      typedef calc_model_data calc_model_data_t;
      typedef wells::connection connection_t;
      typedef rs_mesh_iface mesh_iface_t;
      typedef jacobian_matrix jacobian_matrix_t;

      typedef pvt_oil pvt_oil_t;
      typedef pvt_dead_oil pvt_dead_oil_t;
      typedef pvt_gas pvt_gas_t;
      typedef pvt_water pvt_water_t;

      typedef calc_well_pressure_base calc_well_pressure_t;
      typedef calc_rho_base calc_rho_base_t;
      typedef calc_perf_density_base calc_perf_density_t;
      typedef calc_perf_bhp_base calc_perf_bhp_t;

      typedef connection_iterator connection_iterator_t;
      typedef wells::well_facility_iface well_facility_t;

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

      typedef smart_ptr <well_facility_t>               sp_well_facility_t;
      typedef shared_vector <sp_well_facility_t>        well_facility_list_t;

    public:
      /**
       * \brief  Adds primary connection (perforation) to well and return it
       * \param  i_coord i coordinate of perforation
       * \param  j_coord j coordinate of perforation
       * \param  k_coord k coordinate of perforation
       * \param  n_block Index of block (cell) in mesh for (i, j, k) coordinates
       * \return Created connection
       * */
      virtual sp_connection_t
      add_primary_connection (index_t i_coord, index_t j_coord, index_t k_coord, index_t n_block)
      {
        bs_throw_exception ("PURE CALL");
      }
      /**
       * \brief  Adds secondary connection (perforation) to well and return it
       * \param  i_coord i coordinate of perforation
       * \param  j_coord j coordinate of perforation
       * \param  k_coord k coordinate of perforation
       * \param  n_block Index of block (cell) in mesh for (i, j, k) coordinates
       * \return Created connection
       * */
      virtual sp_connection_t
      add_secondary_connection (index_t i_coord, index_t j_coord, index_t k_coord, index_t n_block)
      {
        bs_throw_exception ("PURE CALL");
      }
      /**
       * \brief  Returns connection (perforation) with n_block
       * \param  n_block Value of block index
       * \return connection instance on success otherwise null pointer
       * */
      virtual sp_connection_t
      get_connection_map (index_t n_block) const
      {
        bs_throw_exception ("PURE CALL");
      }

      /**
       * \brief  Returns iterator for begin of primary and
       *         secondary connections
       * \return Begin iterator
       * */
      virtual connection_iterator_t
      connections_begin () const
      {
        bs_throw_exception ("PURE CALL");
      }

      /**
       * \brief  Returns iterator for end of primary and
       *         secondary connections
       * \return End iterator
       * */
      virtual connection_iterator_t
      connections_end () const
      {
        bs_throw_exception ("PURE CALL");
      }

      /**
       * \brief  Returns true if no any connections
       * \return True if no any connections
       * */
      virtual bool
      is_no_connections () const
      {
        bs_throw_exception ("PURE CALL");
      }

      /**
       * \brief  Returns true if no primary connections
       * \return True if no primary connections
       * */
      virtual bool
      is_no_primary_connections () const
      {
        bs_throw_exception ("PURE CALL");
      }

      /**
       * \brief  Returns first primary connection
       * \return First primary connection or throw exception
       *         is_no_primary_connections == true
       * */
      virtual sp_connection_t
      get_first_connection () const
      {
        bs_throw_exception ("PURE CALL");
      }

      /**
       * \brief  Returns true if connection n_block is
       *         in primary_connection list
       * \param  n_block Number of connection block
       * \return True if connection is in primary_connection list
       * */
      virtual bool
      is_primary_connection (index_t n_block) const
      {
        bs_throw_exception ("PURE CALL");
      }

      /**
       * \brief  Returns true if connection it is in
       *         primary_connection list
       * \param  it Connection iterator
       * \return True if connection is in primary_connection list
       * */
      virtual bool
      is_primary_connection (const connection_iterator_t &it) const
      {
        bs_throw_exception ("PURE CALL");
      }


      /**
       * \brief  Sets coordinates of heel
       * \param  i_coord i coordinate of heel
       * \param  j_coord j coordinate of heel
       * */
      void
      set_coord (index_t i_coord, index_t j_coord);

      /**
       * \brief  Sets BHP reference depth
       * \param  bhp_depth BHP reference depth
       * */
      void
      set_bhp_depth (item_t bhp_depth);

      /**
       * \brief  Sets initial state of well
       * \param  well_state initial well state
       * \param  calc_model
       * */
      void
      set_state (well_state_type well_state, const sp_calc_model_t &calc_model);

      /**
       * \brief  Sets exploitation factor (wefac)
       * \param  exploitation_factor
       * \todo   Should be renamed to wefac
       * */
      void
      set_exploitation_factor (item_t exploitation_factor);

      /**
       * \brief  Sets well controller
       * \param  controller Pointer to well_control instance
       * */
      void
      set_controller (sp_well_controller_t controller);

      /**
       * \brief  Sets well limit operation - operation that
       *         should be performed if well limits reached
       * \param  limit_operation Pointer to well_limit_operation instance
       * */
      void
      set_limit_operation (sp_well_limit_operation_t limit_operation);

      /**
       * \brief  Returns well_controler instance
       * \return well_controller instance
       * */
      sp_well_controller_t
      get_controller () const;

      /**
       * \brief  Returns well_limit_operation instance
       * \return well_limit_operation instance
       * */
      sp_well_limit_operation_t
      get_limit_operation () const;

      /**
       * \brief  Returns name of the well
       * \return Name of the well
       * */
      const std::string &
      get_name () const;

      /**
       * \brief  Sets name of the well
       * \param  name Name of the well
       * */
      void
      set_name (const std::string &name);

      /**
       * \brief  Returns BHP reference depth
       * \return Value of BHP reference depth
       * */
      item_t
      get_bhp_depth () const
        {
          return bhp_depth_;
        }

      /**
       * \brief  Returns BHP value
       * \return Value of BHP
       * */
      item_t
      bhp () const
        {
          return bhp_;
        }

      /**
       * \brief  Returns rate data
       * \return Rate data
       * */
      const rate_data_t &
      rate () const
        {
          return rate_;
        }
      /**
       * \brief  Returns total rate data
       * \return Total rate data
       * */
      const rate_data_t &
      rate_total () const
        {
          return rate_total_;
        }

      /**
       * \brief  Returns production part of rate data
       * \return Production part of rate data
       * */
      const rate_data_inner_t &
      rate_prod () const
      {
        return rate_.prod;
      }
      /**
       * \brief  Returns injection part of rate data
       * \return Injection part of rate data
       * */
      const rate_data_inner_t &
      rate_inj () const
      {
        return rate_.inj;
      }

      /**
       * \brief  Checks if well is open
       * \return True if well is open otherwise false
       * */
      bool
      is_open () const
        {
          if (well_state_.state == well_open)
            return true;
          return false;
        }

      /**
       * \brief  Returns name of the well
       * \return Name of the well
       * */
      const std::string &
      name () const;

      /**
       * \brief  Sets value of BHP
       * \param  bhp Value of BHP
       * */
      void
      set_bhp (item_t bhp);

      /**
       * \todo May be this method is obsolete
       * */
      bool
      fi_check_limits () const;

      /**
       * \brief  Returns calculated reference depth
       * \param  mesh
       * \return Value of calculated reference depth
       * */
      item_t
      get_reference_depth (const sp_mesh_iface_t &mesh) const;

      /**
       * \brief  Returns BHP value
       * \return BHP value
       * */
      item_t
      get_reference_pressure () const;

      /**
       * \brief  Checks connections BHP
       * \param  pressure Array of pressure
       * \return True if any connection have valid BHP value otherwise false
       * */
      bool
      check_connections_bhp (const item_array_t &pressure) const;

      /**
       * \brief  Returns true if well controlled by BHP
       * \return True if well controlled by BHP
       * */
      bool
      is_bhp () const;

      /**
       * \brief  Returns true if well controlled by rate
       * \return True if well controlled by rate
       * */
      bool
      is_rate () const;

      /**
       * \brief  Returns true is well is shuted
       * \return True if well state is shut
       * */
      bool
      is_shut () const;

      /**
       * \brief  Returns well_controller instance
       * \return Pointer to well_controller instance
       * */
      const sp_well_controller_t &
      get_well_controller () const;

      /**
       * \brief  Returns input rate (from well_controller)
       * \return Input rate
       * */
      item_t
      get_input_rate () const;

      /**
       * \brief  Checks well on shut if not shut fills
       *         open_connections_ array
       * \return True if well is shut
       * */
      virtual bool
      check_shut ()
      {
        bs_throw_exception ("PURE CALL");
      }

      /**
       * \brief  Resets init_approx_is_calc_ flag
       * */
      void
      reset_init_approx ();

      /**
       * \brief  Calculates rate and deriv values for well and
       *         well perforations (connections)
       * \param  is_start
       * \param  dt
       * \param  calc_model
       * \param  mesh
       * \param  jmatrix
       * */
      virtual void
      process (bool is_start, double dt, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, sp_jmatrix_t &jmatrix);

      /**
       * \brief  Clears well and well perforations data
       * */
      virtual void
      clear_data ();

      /**
       * \brief  Fills Jacobian rows
       * \param  rows Array of Jacobian Rows
       * */
      virtual void
      fill_rows (index_array_t &rows) const
      {
        bs_throw_exception ("PURE CALL");
      }

      /**
       * \brief  Fills Jacobian colls and values, uses eliminate for
       *         fill values
       * \param  dt
       * \param  block_size
       * \param  rows
       * \param  cols
       * \param  values
       * \param  markers
       * */
      virtual void
      fill_jacobian (double dt, index_t block_size, const index_array_t &rows, index_array_t &cols, rhs_item_array_t &values, index_array_t &markers) const
      {
        bs_throw_exception ("PURE CALL");
      }

      /**
       * \brief  Fills rhs array with rate values
       * \param  dt
       * \param  n_phases
       * \param  is_g
       * \param  is_o
       * \param  is_w
       * \param  rhs Array of rhs values
       * */
      virtual void
      fill_rhs (double dt, index_t n_phases, bool is_g, bool is_o, bool is_w, rhs_item_array_t &rhs) const
      {
        bs_throw_exception ("PURE CALL");
      }

      /**
       * \brief  Restores solution
       * \param  dt
       * \param  p_sol primary solution vector
       * \param  s_sol secondary solution vector
       * \param  block_size size of one block in vectors
       * */
      virtual void
      restore_solution (double dt, const item_array_t &p_sol, const item_array_t &s_sol, index_t block_size);

      /**
       * \brief  Custom init
       * \param  mdl Pointer to calc_model instance
       * \todo   Specify more details
       * */
      virtual void
      custom_init (const sp_calc_model_t &mdl);

      /**
       * \brief  Performs actions before start of each large step
       * \param  calc_model
       * \param  mesh
       * */
      virtual void
      pre_large_step (const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh);

      /**
       * \brief  Performs actions before start of each small step
       * */
      virtual void
      pre_small_step ();

      /**
       * \brief  Performs actions before start of each newton step
       * */
      virtual void
      pre_newton_step ();

      /**
       * \brief  Restores 'internal-state' of well if small step failed
       * */
      virtual void
      restart_small_step ();

      /**
       * \brief  Restores 'internal-state' of well if newton step failed
       * */
      virtual void
      restart_newton_step ();

      /**
       * \brief  Calculates rate and deriv values for well and
       *         well perforations (connections). Process function
       *         should call
       *         pre_process_facilities -> process -> post_process_facilities
       *         To guarantee this we implement process function in process_impl.
       * \param  is_start
       * \param  dt
       * \param  calc_model
       * \param  mesh
       * \param  jmatrix
       * */
      virtual void
      process_impl (bool is_start, double dt, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh, sp_jmatrix_t &jmatrix);

      /**
       * \brief  add well facility to facility list
       * */
      void
      add_well_facility (const sp_well_facility_t &facility);

      /**
       * \brief  delete well facility from facility list by iterator
       * */
      void
      delete_well_facility (typename well_facility_list_t::iterator iter);

      /**
       * \brief  return well facility list
       * */
      well_facility_list_t &
      get_facility_list ()
        {
          return well_facility_list_;
        }

      /**
       * \brief  well dtor
       * */
      virtual ~well ();

      /**
       * \brief  well ctor
       * \param  name of new well
       * \todo   Obsolete
       * */
      well (const std::string &well_name);

    public:
      //! blue-sky type declaration
      BLUE_SKY_TYPE_DECL_T (well);

    protected:

      /**
       * \brief  Shuts well, sets bulkp and bhp values
       *         for each perforation (connection) to
       *         value of calc_model pressure
       * \param  calc_model
       * */
      void
      shut_well (const sp_calc_model_t &well);

      /**
       * \brief  Computes perforations (connections) factor
       * \param  internal_constants
       * \param  params
       * \param  mesh
       * \param  perm
       * \param  ntg
       * \param  ro_calc_flag
       * */
      void
      compute_connection_factor (const physical_constants &internal_constants,
                                 const sp_params_t &params,
                                 const sp_mesh_iface_t &mesh,
                                 const item_array_t &perm,
                                 const item_array_t &ntg,
                                 bool ro_calc_flag);


      /**
       * \brief  Performs facilities actions before process
       * */
      void
      pre_process_facilities ();

      /**
       * \brief  Performs facilities actions after process
       * */
      void
      post_process_facilities ();

    public:
      DECLARE_EVENT_LIST_V2 (well,
        ((connection_change, (const sp_well_t &, const index_t &), 2))
      );

    public:

      std::string                 name_;                      //!< Name of the well

      auto_value <index_t, -1>    i_coord_;                   //!< i coordinate of hell (well head)
      auto_value <index_t, -1>    j_coord_;                   //!< j coordinate of hell (well head)

      sp_well_controller_t        well_controller_;           //!< well_controller instance

//  private:
    public:
      auto_value <item_t, -1>     bhp_depth_;                 //!< BHP reference depth
      auto_value <item_t, 1>      exploitation_factor_;       //!< Exploitation factor (wefac)

      well_state_t                well_state_;                //!< State of well
      well_state_t                saved_well_state_;          //!< State of well on begin of small (or large) step
      well_state_t                saved_niter_well_state_;    //!< State of well on begin of newton step

      sp_well_limit_operation_t   well_limit_operation_;      //!< well_limit_operation instance

      rate_data_t                 rate_;                      //!< Rate data
      rate_data_t                 rate_total_;                //!< Total rate data (from begin of simulation)
      rate_data_t                 rate_rc_;                   //!< Rate data in reservoir conditions

      item_t                      bhp_;                       //!< BHP value
      item_t                      gor_;                       //!< gas_oil_ratio value

    protected:

      auto_value <bool, false>    init_approx_is_calc_;       //!< Flag that specified should initial approximation calculated
      auto_value <item_t>         input_reference_depth_;     //!< \todo Should be removed

      sp_calc_well_pressure_t     calc_well_pressure_;        //!< Instance of object that calculates well pressure
      sp_calc_rho_t               calc_rho_;                  //!< \todo Should be removed
      sp_calc_perf_density_t      calc_perf_density_;         //!< Instance of object that calculates well perforations density
      sp_calc_perf_bhp_t          calc_perf_bhp_;             //!< Instance of object that calculates well perforations BHP

      well_facility_list_t        well_facility_list_;        //!< Well facilities list
    };

  /**
   * \class well_factory
   * \brief Creates wells and perforations
   * \todo  create_connection should be removed because
   *        perforations now created in add_connection
   * */
  class BS_API_PLUGIN well_factory : public objbase
    {
    public:

      typedef well well_t;
      typedef wells::connection connection_t;
      typedef smart_ptr <well_t>                sp_well_t;
      typedef smart_ptr <connection_t, true>    sp_connection_t;

    public:

      /**
       * \brief  well_factory dtor
       * */
      virtual ~well_factory () {}

      /**
       * \brief  Creates well with name well_name and in group group_name
       * \param  group_name Name of group to which the well belongs
       * \param  well_name Name of well
       * \return Instance of well
       * */
      virtual sp_well_t
      create_well (const std::string &group_name, const std::string &well_name) const;

      /**
       * \todo   Obsolete, should be removed
       * */
      virtual sp_connection_t
      create_connection () const;

      //! blue-sky type declaration
      BLUE_SKY_TYPE_DECL_T (well_factory);
    };

  /**
   * \brief  Registers well types in blue-sky kernel
   * \param  pd plugin_desriptor
   * \return True if all types registered successfully
   * */
  bool
  calc_well_register_types (const blue_sky::plugin_descriptor &pd);

  /**
   * \brief  Registers well factory types in blue-sky kernel
   * \param  pd plugin_desriptor
   * \return True if all types registered successfully
   * */
  bool
  well_factory_register_type (const blue_sky::plugin_descriptor &pd);

} // namespace blue_sky


#endif  // #ifndef BS_CALC_WELL_H_
