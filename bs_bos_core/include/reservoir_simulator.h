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
#include "hdm.h"
#include "event_manager.h"
#include "jacobian.h"

namespace blue_sky
{
  class reservoir;
  class calc_model;
  struct main_loop_calc_base;
  class data_storage_interface;

  /**
   * \class reservoir_simulator
   * \brief main class of reservoir_simulator implementation,
   *        implements bs_node
   * */
  class BS_API_PLUGIN reservoir_simulator : public bs_node
    {
    public:
      // typedefs
      typedef bs_node                                 base_t;
      typedef smart_ptr <reservoir_simulator, true >  sp_this_t;

      typedef smart_ptr <hdm, true>                   sp_hdm_t;              //!< smart_ptr to hdm type
      typedef smart_ptr <rs_mesh_iface, true >        sp_mesh_iface_t;      //!< smart_ptr to rs_mesh_iface type
      typedef smart_ptr <lsolver_iface, true>         sp_solver_t;          //!< smart_ptr to solver_t type
      typedef smart_ptr <event_manager, true >        sp_em_t;              //!< smart_ptr to event_manager type
      typedef event_manager::sp_event_base_list       sp_event_base_list_t; //!< list of scheduled events
      typedef smart_ptr < calc_model, true >          sp_calc_model_t;      //!< smart_ptr to calc_model type
      typedef smart_ptr <reservoir, true>             sp_reservoir_t;       //!< smart_ptr to reservoir type
      typedef data_storage_interface                  facility_storage_t;   //!< facility_storage type
      typedef smart_ptr <facility_storage_t, true>    sp_facility_storage_t;//!< smart_ptr to facility_storage type
      typedef smart_ptr <jacobian, true>              sp_jacobian_t;        //!< smart_ptr to jacobian type

    public:
      /**
       * \brief  setups pointer to mesh
       * \param  mesh pointer to mesh
       * */
      void                  
      set_mesh (const sp_mesh_iface_t&);

      /**
       * \brief  returns pointer to hdm instance
       * */
      sp_hdm_t 
      get_hdm () const;

      /**
       * \brief  returns pointer to event_manager instance
       * */
      sp_em_t 
      get_event_manager () const;

      /**
       * \brief  returns pointer to calc_model instance
       * */
      sp_calc_model_t 
      get_calc_model () const;

      /**
       * \brief  returns pointer to mesh instance
       * */
      sp_mesh_iface_t
      get_mesh () const;

      /**
       * \brief  returns pointer to jacobian instance
       * */
      sp_jacobian_t 
      get_jacobian () const;

      /**
       * \brief  returns pointer to reservoir instance
       * */
      sp_reservoir_t 
      get_reservoir () const;

      /**
       * \brief returns pointer to facility_storage
       * */
      smart_ptr <facility_storage_t> 
      get_facility_storage () const;

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
      main_loop_calc_base *
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

      BLUE_SKY_TYPE_DECL (reservoir_simulator)


    protected:
      sp_hdm_t                hdm_;               //!< pointer to hdm instance
      sp_em_t                 em;                 //!< pointer to event_manager instance
      sp_calc_model_t         cm;                 //!< pointer to calc_model instance
      sp_reservoir_t          reservoir_;         //!< pointer to reservoir instance
      sp_facility_storage_t   facility_storage_;  //!< pointer to facility_storage instance
      sp_jacobian_t           jacobian_;          //!< pointer to jacobian instance
      std::string             model_filename_;    //!< name of model file

      smart_ptr <main_loop_calc_base, false> mloop;  //!< pointer to main_loop_calc instance

    public:

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
        ((before_ready, (), 0))
        ((after_ready, (), 0))
      );
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
