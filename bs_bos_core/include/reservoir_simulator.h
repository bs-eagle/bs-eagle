/**
 *       \file  reservoir_simulator.h
 *      \brief  main file of reservoir simulator
 *     \author  Nikonov Max
 *       \date  10.11.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef RESERVOIR_SIMULATOR_H
#define RESERVOIR_SIMULATOR_H

#include "simulator_events.h"
#include "hydrodynamic_model.h"
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
   * \class reservoir_simulator
   * \brief main class of reservoir_simulator implementation,
   *        implements bs_node
   * */
  template <class strategy_t>
  class BS_API_PLUGIN reservoir_simulator : public bs_node
    {
    public:
      // typedefs
      typedef bs_node                               base_t;
      typedef reservoir_simulator < strategy_t >    this_t;               //!< shortname for this type
      typedef smart_ptr < this_t, true >            sp_this_t;            //!< smart pointer to this_t

      typedef hydrodynamic_model < strategy_t >           dm_t;                 //!< hydrodynamic_model type
      typedef smart_ptr < dm_t, true >              sp_hdm_t;              //!< smart_ptr to hydrodynamic_model type

      typedef rs_mesh_iface < strategy_t >          mesh_iface_t;         //!< rs_mesh_iface type
      typedef smart_ptr < mesh_iface_t, true >      sp_mesh_iface_t;      //!< smart_ptr to rs_mesh_iface type

      typedef linear_solver_base < strategy_t >     solver_t;             //!< linear_solver type
      typedef smart_ptr < solver_t, true>           sp_solver_t;          //!< smart_ptr to solver_t type

      typedef event_manager <strategy_t>            em_t;                 //!< event_manager type
      typedef smart_ptr < em_t, true >              sp_em_t;              //!< smart_ptr to event_manager type
      typedef typename em_t::sp_event_base_list     sp_event_base_list_t; //!< list of scheduled events

      typedef calc_model < strategy_t >             calc_model_t;         //!< calc_model type
      typedef smart_ptr < calc_model_t, true >      sp_calc_model_t;      //!< smart_ptr to calc_model type

      typedef reservoir <strategy_t>                reservoir_t;          //!< reservoir type
      typedef smart_ptr <reservoir_t, true>         sp_reservoir_t;       //!< smart_ptr to reservoir type

      typedef data_storage_interface                facility_storage_t;   //!< facility_storage type
      typedef smart_ptr <facility_storage_t, true>  sp_facility_storage_t;//!< smart_ptr to facility_storage type

      typedef jacobian <strategy_t>                 jacobian_t;           //!< jacobian type
      typedef smart_ptr <jacobian_t, true>          sp_jacobian_t;        //!< smart_ptr to jacobian type

      typedef typename strategy_t::item_t           item_t;               //!< item type
      typedef typename strategy_t::index_t          index_t;              //!< index type
      typedef typename strategy_t::item_array_t     item_array_t;         //!< array of items type
      typedef typename strategy_t::index_array_t    index_array_t;        //!< array of indexes type

      //! transmissibility multipliers calculator type
      typedef trans_multipliers_calc <strategy_t>   trans_multipliers_calc_t; 
      typedef main_loop_calc_base <strategy_t>      main_loop_calc_t;     //!< base type of main_loop_calc type

      typedef keyword_manager <strategy_t>          keyword_manager_t;    //!< keyword_manager type
      typedef smart_ptr <keyword_manager_t, true >  sp_keyword_manager_t; //!< smart_ptr to keyword_manager type

      typedef jacobian_matrix <strategy_t>          jmatrix_t;            //!< jacobian_matrix type
      typedef smart_ptr <jmatrix_t, true>           sp_jmatrix_t;         //!< smart_ptr to jacobian_matrix type

    public:
      /**
       * \brief  setups pointer to mesh
       * \param  mesh pointer to mesh
       * */
      void                  
      set_mesh (const sp_mesh_iface_t&);

      /**
       * \brief  returns pointer to hydrodynamic_model instance
       * */
      const sp_hdm_t &
      get_hydrodynamic_model () const;

      /**
       * \brief  returns pointer to event_manager instance
       * */
      const sp_em_t &
      get_event_manager () const;

      /**
       * \brief  returns pointer to calc_model instance
       * */
      const sp_calc_model_t &
      get_calc_model () const;

      /**
       * \brief  returns pointer to mesh instance
       * */
      const sp_mesh_iface_t &
      get_mesh () const;

      /**
       * \brief  returns pointer to jacobian instance
       * */
      const sp_jacobian_t &
      get_jacobian () const;

      /**
       * \brief  returns pointer to reservoir instance
       * */
      const sp_reservoir_t &
      get_reservoir () const;

      /**
       * \brief  inits reservoir simulator
       * */
      void 
      init();

      /**
       * \brief  creates (if needed) main_loop_calc instance 
       *         and starts main simulation loop
       * \return may throw exception
       * */
      void 
      main_loop ();

      /**
       * \brief  creates (if needed) main_loop_calc instance 
       *         and returns main_loop_calc
       * \return may throw exception
       * */
      main_loop_calc_t *
      get_main_loop ();

      /**
       * \brief  reservoir_simulator dtor
       * */
      virtual ~reservoir_simulator();

      /**
       * \brief  sets subnode in tree by subnode name
       * \param  name name of subnode
       * \param  obj pointer to new subnode object
       * */
      void 
      set_subnode_in_tree (const std::string &name, sp_obj obj)
      {
        bs_node::erase (name);
        bs_node::insert (obj, name, false);
      }

      /**
       * \brief  load model and starts main simulation loop
       * \param  path path to the model file
       * \return may throw exception
       * */
      void 
      simulate (const std::string &path);

      /**
       * \brief  load model and init it
       * \param  path path to the model file
       * \return may throw exception
       * */
      void 
      read_keyword_file_and_init (const std::string &path);

      /**
       * \brief  default actions that should be executed
       *         before start of the following large simulation
       *         step
       * \param  event_list list of scheduled events that should be
       *         executed before start of the step
       * */
      void 
      pre_large_step (const sp_event_base_list_t &event_list);

      /**
       * \brief  returns name of model file
       * */
      std::string
      model_filename () const;

      BLUE_SKY_TYPE_DECL_T(reservoir_simulator)

      //! storage for event params
      typedef bs_array <std::string> signal_params_t;

      //! declaration of reservoir_simulator events
      DECLARE_EVENT_LIST_V2 (reservoir_simulator,
        ((begin, (clock_t), 1))
        ((newton_iter_fail, (), 0))
        ((newton_iter_success, (), 0))
        ((before_fi_operator, (), 0))
        ((before_jacobian_setup, (), 0))
        ((before_jacobian_solve, (), 0))
        ((before_restore_solution, (), 0))
        ((simulation_start, (), 0))
        ((large_step_start, (), 0))
        ((small_step_start, (), 0))
        ((end, (clock_t), 1))
        ((post_read, (), 0))
        ((pre_read, (sp_this_t), 1))
      );

      sp_hdm_t                 hdm;                 //!< pointer to hydrodynamic_model instance
      sp_em_t                 em;                 //!< pointer to event_manager instance
      sp_calc_model_t         cm;                 //!< pointer to calc_model instance
      sp_mesh_iface_t         mesh;               //!< pointer to mesh instance
      sp_reservoir_t          reservoir_;         //!< pointer to reservoir instance
      sp_facility_storage_t   facility_storage_;  //!< pointer to facility_storage instance
      sp_jacobian_t           jacobian_;          //!< pointer to jacobian instance
      sp_keyword_manager_t    keyword_manager_;   //!< pointer to keyword manager instance
      std::string             model_filename_;    //!< name of model file

      smart_ptr <main_loop_calc_t, false> mloop;  //!< pointer to main_loop_calc instance
    };

  /**
   * \brief  registers reservoir_simulator types in blue-sky kernel
   * \param  pd plugin descriptor
   * \return true if all types registered successfully
   * */
  bool
  reservoir_simulator_register_types (const plugin_descriptor &pd);
}

#endif // RESERVOIR_SIMULATOR_H
