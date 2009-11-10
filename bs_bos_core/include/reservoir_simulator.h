/**
	\file reservoir_simulator.h
	\brief main reservoir simulator of bos core
	\author Nikonov Max
*/

#ifndef RESERVOIR_SIMULATOR_H
#define RESERVOIR_SIMULATOR_H

#include "simulator_events.h"
#include "data_manager.h"
#include "event_manager.h"
#include "jacobian.h"

namespace blue_sky
{

  template <typename strategy_t>
  class reservoir;

  template <typename strategy_t>
  class calc_model;

  template <typename strategy_t>
  class reservoir;

  template <typename strategy_t>
  struct BS_API_PLUGIN trans_multipliers_calc;

  template <typename strategy_t>
  struct main_loop_calc_base;

  class data_storage_interface;

  template <typename strategy_t>
  class keyword_manager;

  /**
  	\class reservoir_simulator
  	\brief main bos class for data and process manipulate
  */
  template <class strategy_t>
  class BS_API_PLUGIN reservoir_simulator : public bs_node
    {
      //private:
      //struct mstatus_traits;                                      //!< tree node sorting traits

    public:
      // typedefs
      typedef reservoir_simulator < strategy_t >  this_t;         //!< this type
      typedef smart_ptr < this_t, true >          sp_this_t;      //!< smart pointer to this_t

      typedef data_manager < strategy_t >         dm_t;           //!< data_manager type
      typedef smart_ptr < dm_t, true >            sp_dm_t;        //!< smart_ptr to data_manager type

      typedef rs_mesh_iface < strategy_t >              mesh_iface_t;         //!< rs_mesh_iface type
      typedef smart_ptr < mesh_iface_t, true >          sp_mesh_iface_t;      //!< smart_ptr to rs_mesh_iface type

      typedef linear_solver_base < strategy_t >   solver_t;       //!< linear_solver type
      typedef smart_ptr < solver_t, true>         sp_solver_t;    //!< smart_ptr to solver_t type

      typedef event_manager <strategy_t>          em_t;           //!< event_manager type
      typedef smart_ptr < em_t, true >            sp_em_t;        //!< smart_ptr to event_manager type
      typedef typename em_t::sp_event_base_list   sp_event_base_list_t;

      typedef calc_model < strategy_t >           calc_model_t;   //!< calc_model type
      typedef smart_ptr < calc_model_t, true >    sp_calc_model_t;//!< smart_ptr to calc_model type

      typedef reservoir <strategy_t>							reservoir_t;		//!< short name
      typedef smart_ptr <reservoir_t, true>       sp_reservoir_t; //!< short name

      typedef data_storage_interface              facility_storage_t;
      typedef smart_ptr <facility_storage_t, true>  sp_facility_storage_t;

      typedef jacobian <strategy_t>							  jacobian_t;		//!< short name
      typedef smart_ptr <jacobian_t, true>        sp_jacobian_t; //!< short name

      typedef typename strategy_t::item_t         item_t;         //!< item type
      typedef typename strategy_t::index_t        index_t;        //!< index type
      typedef typename strategy_t::item_array_t   item_array_t;   //!< array of items type
      typedef typename strategy_t::index_array_t  index_array_t;  //!< array of indexes type

      typedef trans_multipliers_calc <strategy_t> trans_multipliers_calc_t;
      typedef main_loop_calc_base <strategy_t>    main_loop_calc_t;

      typedef keyword_manager <strategy_t>             keyword_manager_t;
      typedef smart_ptr <keyword_manager_t, true >	    sp_keyword_manager_t;

      typedef jacobian_matrix <strategy_t>        jmatrix_t;
      typedef smart_ptr <jmatrix_t, true>         sp_jmatrix_t;

      // methods
      void                  set_mesh (const sp_mesh_iface_t&);
      const sp_dm_t         &get_data_manager () const;
      const sp_em_t         &get_event_manager () const;
      const sp_calc_model_t &get_calc_model () const;
      const sp_mesh_iface_t &get_mesh () const;
      const sp_jacobian_t   &get_jacobian () const;
      const sp_reservoir_t  &get_reservoir () const;

      void init();

      // main loop
      void main_loop ();

      main_loop_calc_t *
      get_main_loop ();

      //void main_loop_iteration (const typename event_manager::event_map::iterator & it);

      //! dtor
      virtual ~reservoir_simulator();

      void set_subnode_in_tree (const std::string &name, sp_obj obj)
      {
        bs_node::erase (name);
        bs_node::insert (obj, name, false);
      }

      void simulate (const std::string &path);
      void read_keyword_file_and_init (const std::string &path);

      void 
      pre_large_step (const sp_event_base_list_t &event_list);

      std::string
      model_filename () const;

      //! blue_sky class declarations
      BLUE_SKY_TYPE_DECL_T(reservoir_simulator)

      typedef bs_array <std::string> signal_params_t;

      DECLARE_EVENT_LIST (reservoir_simulator,
                          signal_params_t,
                          12,
                          ((begin, (clock_t), 1)
                           , (newton_iter_fail, (), 0)
                           , (newton_iter_success, (), 0)
                           , (before_fi_operator, (), 0)
                           , (before_jacobian_setup, (), 0)
                           , (before_jacobian_solve, (), 0)
                           , (before_restore_solution, (), 0)
                           , (simulation_start, (), 0)
                           , (large_step_start, (), 0)
                           , (small_step_start, (), 0)
                           , (end, (clock_t), 1)
                           , (post_read, (), 0)
                          ));


      //private:
      sp_dm_t                 dm;
      sp_em_t                 em;
      sp_calc_model_t         cm;
      sp_mesh_iface_t         mesh;
      sp_reservoir_t          reservoir_;
      sp_facility_storage_t   facility_storage_;
      sp_jacobian_t           jacobian_;
      sp_keyword_manager_t    keyword_manager_;
      std::string             model_filename_;

      smart_ptr <main_loop_calc_t, false> mloop;
    };

  bool
  reservoir_simulator_register_types (const plugin_descriptor &pd);
}

#endif // RESERVOIR_SIMULATOR_H
