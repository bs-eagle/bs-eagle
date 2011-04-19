/**
 *       \file  reservoir.h
 *      \brief  reservoir is a storage for and manager of facilities. also
 *              stores reservoir rates.
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  16.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_RESERVOIR_H_
#define BS_RESERVOIR_H_

#include "calc_well.h"
#include "event_filter.h"
#include "jacobian.h"

#ifdef _HDF5
#include "bs_hdf5_storage.h"
#include "rs_smesh_iface.h"
#endif

#include <boost/shared_ptr.hpp>

namespace blue_sky
  {

  class calc_model;
  class reservoir_simulator;
  class facility_manager;
  class well;
  class well_factory;
  class data_storage_interface;
  class fi_params;

  namespace wells
  {
    class well_controller;
    class well_rate_control;
    class connection;
    class well_controller_factory;
    class well_limit_operation_factory;
    class well_limit_operation;
  }

  struct data_saver;

  /**
   * \class reservoir
   * \brief storage for and manager of facilities. also
   *        stores reservoir rates.
   * */
  class BS_API_PLUGIN reservoir : public objbase
    {
    public:

      typedef t_double                                        item_t;                             //!< item value (floating point)
      typedef t_long                                          index_t;                            //!< index value (integral type)
      typedef spv_double                                      item_array_t;                       //!< type for array of item_t values
      typedef spv_float                                       rhs_item_array_t;                   //!< type for array of rhs_item_t values
      typedef spv_long                                        index_array_t;                      //!< type for array of index_t values

      typedef facility_base                                   facility_t;                         //!< base type for facilities
      typedef well                                            well_t;                             //!< base type for wells
      typedef wells::well_controller                          well_controller_t;                  //!< well_controller type
      typedef wells::well_limit_operation                     well_limit_operation_t;             //!< well_llimit_operation type
      typedef wells::well_rate_control                        well_rate_control_t;                //!< well_rate_control type
      typedef wells::connection                               connection_t;                       //!< base type for well connections (perforations)
      typedef facility_manager                                facility_manager_t;                 //!< facility_manager type
      typedef wells::well_controller_factory                  controller_factory_t;               //!< type for well_controller factory
      typedef wells::well_limit_operation_factory             limit_operation_factory_t;          //!< type for well_limit_operation factory
      typedef well_factory                                    well_factory_t;                     //!< type for well factory
      typedef calc_model                                      calc_model_t;                       //!< calc_model type
      typedef rs_mesh_iface                                   mesh_iface_t;                       //!< rs_mesh_iface type
      typedef jacobian                                        jacobian_t;                         //!< jacobian type
      typedef jac_matrix_iface                                jacobian_matrix_t;                  //!< jacobian_matrix type

      typedef reservoir_simulator                             reservoir_simulator_t;              //!< reservoir_simulator type
      typedef rate_data                                       rate_data_t;                        //!< type for rate data holder
      typedef rate_data_t::rate_data_inner                    rate_data_inner_t;                  //!< type for internal data of rate data holder

      typedef data_saver                                      data_saver_t;

      typedef smart_ptr <reservoir_simulator_t, true>         sp_top_t;                           //!< smart_ptr to reservoir_simulator type

      typedef smart_ptr <well_t, true>                        sp_well_t;                          //!< smart_ptr to well type
      typedef smart_ptr <well_controller_t, true>             sp_well_controller_t;               //!< smart_ptr to well_controller type
      typedef smart_ptr <well_limit_operation_t, true>        sp_well_limit_operation_t;          //!< smart_ptr to well_limit_operation type
      typedef smart_ptr <well_rate_control_t, true>           sp_rate_control_t;                  //!< smart_ptr to rate_control type
      typedef smart_ptr <connection_t, true>                  sp_connection_t;                    //!< smart_ptr to connection type
      typedef smart_ptr <mesh_iface_t, true>                  sp_mesh_iface_t;                    //!< smart_ptr to rs_mesh_iface type
      typedef smart_ptr <calc_model_t, true>                  sp_calc_model_t;                    //!< smart_ptr to calc_model type
      typedef smart_ptr <jacobian_t, true>                    sp_jacobian_t;                      //!< smart_ptr to jacobian type
      typedef smart_ptr <jacobian_matrix_t, true>             sp_jacobian_matrix_t;               //!< smart_ptr to jacobian_matrix type
      typedef smart_ptr <jacobian_matrix_t, true>             sp_jmatrix_t;                       //!< smart_ptr to jacobian_matrix type

      typedef smart_ptr <facility_manager_t, true>            sp_facility_manager_t;              //!< smart_ptr to facility_manager type
      typedef smart_ptr <data_storage_interface, true>        sp_storage_t;                       //!< smart_ptr to data_storage_interface type

      typedef smart_ptr <controller_factory_t, true>          sp_well_controller_factory_t;       //!< smart_ptr to well_controller factory type
      typedef smart_ptr <limit_operation_factory_t, true>     sp_well_limit_operation_factory_t;  //!< smart_ptr to well_limit_operation factory type
      typedef smart_ptr <well_factory_t, true>                sp_well_factory_t;                  //!< smart_ptr to well factory type

      typedef smart_ptr <fi_params, true>                     sp_params_t;                        //!< smart_ptr to fi_params type
      typedef smart_ptr <event_filter, true>                  sp_event_filter_t;                  //!< smart_ptr to event_filter type

#ifdef _HDF5
      typedef smart_ptr <bs_hdf5_storage, true>               sp_bs_hdf5_storage;                 //!< smart_ptr to hdf5_storage type
#endif

    public:

      //! blue-sky type declaration
      BLUE_SKY_TYPE_DECL (reservoir);

      /**
       * \brief  reservoir dtor
       * */
      ~reservoir ();

    public:

      /**
       * \brief  returns well with name well_name from group group_name from facility_list_
       * \param  group_name name of well group
       * \param  well_name name of well
       * \return well instance if well exists otherwise null pointer
       * */
      sp_well_t                 
      get_well (const std::string &group_name, const std::string &well_name) const;

      /**
       * \brief  returns well with name well_name from default group from facility_list_
       * \param  well_name name of well
       * \return well instance if well exists otherwise null pointer
       * */
      sp_well_t                 
      get_well (const std::string &well_name) const;

      /**
       * \brief  creates well with name well_name in group group_name
       * \param  group_name name of well group
       * \param  well_name name of well
       * \return well instance if created successfully otherwise null pointer
       * */
      sp_well_t                 
      create_well (const std::string &group_name, const std::string &well_name);

      /**
       * \brief  creates well_controller for given well
       * \param  owner_well well instance that will own created well_controller
       * \return well_controller instance if created successfully otherwise null pointer
       * */
      sp_well_controller_t      
      create_well_controller (const sp_well_t &owner_well);

      /**
       * \brief  creates well_limit_operation for given well
       * \param  owner_well well instance that will own created well_limit_operation
       * \param  operation type of well_limit_operation
       * \return well_limit_operation instance if created successfully otherwise null pointer
       * */
      sp_well_limit_operation_t 
      create_well_limit_operation (const sp_well_t &owner_well, wells::limit_operation_type operation);

      /**
       * \brief  creates well connection
       * \return connection instance if created successfully otherwise null pointer
       * */
      sp_connection_t
      create_connection ();

      /**
       * \brief  saves facilities data in data storage
       * \param  storage data storage
       * */
      void                      
      save_data (const sp_storage_t &storage) const;

      /**
       * \todo obsolete, should be removed
       * */
      item_t                    
      pressure () const;

      /**
       * \todo obsolete, should be removed
       * */
      bool                      
      check_limits (const sp_params_t &params) const;

      /**
       * \brief  for each facility calls pre_large_step
       * \param  calc_model
       * \param  mesh
       * */
      void                      
      pre_large_step (const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh);

      /**
       * \brief  for each facility calls pre_small_step
       * */
      void                      
      pre_small_step ();

      /**
       * \brief  for each facility calls pre_newton_step
       * */
      void                      
      pre_newton_step ();

      /**
       * \brief  for each facility calls restart_small_step
       * */
      void                      
      restart_small_step ();

      /**
       * \brief  for each facility calls restart_newton_step
       * */
      void                      
      restart_newton_step ();

      /**
       * \brief  inits jacobian
       * \param  jmx jacobian_matrix that should be inited
       * \param  n_cells number of cells in mesh
       * */
      void                      
      init_jacobian (const sp_jmatrix_t &jmx, index_t n_cells);

      /**
       * \brief  ends building of jacobian
       * \param  dt current dt
       * \param  calc_model
       * \param  jacobian holds jacobian_matrix that builded
       * \return 
       * */
      void                      
      end_jacobian (BS_SP (jacobian) &jacobian, t_double dt, t_long block_size, t_long cells);

      /**
       * \brief  for each facility calls restore_solution
       * \param  dt
       * \param  p_sol primary solution vector
       * \param  s_sol secondary solution vector
       * \param  block_size size of one block in vectors
       * */
      void                      
      restore_wells_solution (double dt, const spv_double &p_sol, const spv_double &s_sol, index_t block_size);

      /**
       * \brief  for each facility calls process
       * \param  istart
       * \param  dt
       * \param  calc_model
       * \param  mesh
       * \param  jmatrix
       * \return 
       * */
      void                      
      calc_wells (int istart, 
        double dt, 
        const BS_SP (calc_model) &calc_model, 
        const BS_SP (rs_mesh_iface) &mesh, 
        BS_SP (jacobian) &jacobian);

      /**
       * \brief  for each facility calls fill_rhs
       * \param  dt
       * \param  calc_model
       * \param  rhs
       * \param  update_after_gauss_elimination
       * */
      void                      
      fill_rhs_wells (double dt, const sp_calc_model_t &calc_model, rhs_item_array_t &rhs, bool update_after_gauss_elimination) const;

      /**
       * \brief  returns facility list
       * \return facility_list
       * */
      sp_facility_manager_t     
      get_facility_list () const;

      /**
       * \brief  specifies what well should be filtered be event_filter
       * \param  well_name name of well
       * */
      void                      
      add_filter_well (const std::string &well_name);

      /**
       * \brief  returns pointer to event_filter instance
       * \return pointer to event_filter instance
       * */
      const sp_event_filter_t &
      get_event_filter () const;

      /**
       * \brief  returns pointer to well factory instance
       * \return pointer to well factory instance
       * */
      sp_well_factory_t                 
      get_well_factory () const;

      /**
       * \brief  returns pointer to well_controller factory instance
       * \return pointer to well_controller factory instance
       * */
      sp_well_controller_factory_t      
      get_well_controller_factory () const;

      /**
       * \brief  returns pointer to well_limit_operation factory
       * \return pointer to well_limit_operation factory
       * */
      sp_well_limit_operation_factory_t 
      get_well_limit_operation_factory () const;

      /**
       * \brief  sets well factory
       * \param  factory pointer to well factory instance
       * */
      void                      
      set_well_factory (const sp_well_factory_t &factory);

      /**
       * \brief  sets well_controller factory
       * \param  pointer to well_controller factory instance
       * */
      void                      
      set_well_controller_factory (const sp_well_controller_factory_t &factory);

      /**
       * \brief  sets well_limit_operation factory
       * \param  pointer to well_limit_operation factory instance
       * */
      void                      
      set_well_limit_operation_factory (const sp_well_limit_operation_factory_t &factory);

#ifdef _HDF5
      /**
       * \todo describe
       * */
      void                      
      open_hdf5_file (const std::string &filename) const;

      /**
       * \todo describe
       * */
      void                      
      close_hdf5_file () const;

      /**
       * \todo describe
       * */
      void                      
      write_step_to_hdf5 (const sp_calc_model_t &calc_model, 
        const sp_mesh_iface_t &mesh, 
        const BS_SP (jacobian) &jmx, 
        int, 
        int, 
        item_t time) const;

      /**
       * \todo describe
       * */
      void                      
      write_mesh_to_hdf5 (const smart_ptr <rs_mesh_iface, true> &mesh) const;

      /**
       * \todo describe
       * */
      const smart_ptr<bs_hdf5_storage, true> 
      get_hdf5_file () const {return hdf5;}
#endif

      /**
       * \brief  Writes data of current time-step to storage (HDF5 for example)
       * \param  calc_model
       * \param  mesh
       * \param jmx
       * \param 
       * \param
       * \param time Current time-step
       * */
      void
      write_step_to_storage (const sp_calc_model_t &calc_model, 
        const sp_mesh_iface_t &mesh, 
        const BS_SP (jacobian) &jacobian, 
        size_t large_time_step_num, 
        size_t total_time_step_num, 
        double time);

      /**
       * \brief  Writes initial mesh status to storage
       * \param  mesh
       * */
      void
      write_mesh_to_storage (const sp_mesh_iface_t &mesh);

      /**
       * \brief  Writes starting date to storage
       * \param  date
       * */
      void
      write_starting_date_to_storage (const boost::posix_time::ptime &date);

      /**
       * \brief  Opens data storage
       * \param  name Name of the data storage
       * */
      void
      open_storage (const std::string &name);

      /**
       * \brief  returns rate data
       * \return rate data
       * */
      const rate_data_t &
      rate () const
      {
        return rate_;
      }

    private:
      /**
       * \brief  inits rows for builded jacobian
       * \param  rows arrays of row indexes that should be initialized
       * */
      void                      
      init_rows (index_array_t &rows) const;

    public:

      rate_data_t                         rate_;                              //!< reservoir rate
      rate_data_t                         rate_rc_;                           //!< reservoir rate in reservoir conditions
      rate_data_t                         rate_wefac_;                        //!< rate_ * wefac
      rate_data_t                         rate_rc_wefac_;                     //!< rate_rc_ * wefac
      rate_data_t                         rate_initial_;                      //!< inital rate on current step
      rate_data_t                         rate_total_;                        //!< total rate

    private:

      sp_facility_manager_t               facility_list_;                     //!< manager for facilities (wells, unnamed wells and other)
      sp_well_factory_t                   well_factory_;                      //!< factory for wells and connections
      sp_well_controller_factory_t        well_controller_factory_;           //!< factory for well controllers
      sp_well_limit_operation_factory_t   well_limit_operation_factory_;      //!< factory for well limit operations

      sp_event_filter_t                   event_filter_;                      //!< events filter

      stdv_long                           markers_;                           //!< markers, used to build jacobian

#ifdef _HDF5
      sp_bs_hdf5_storage                  hdf5;                               //!< pointer to hdf5_storage instance
#endif

      boost::shared_ptr <data_saver_t>    data_saver_;
    };


} // namespace blue_sky

#endif  // #ifndef BS_RESERVOIR_H_

