/**
 *       \file  calc_model.h
 *      \brief  calc_model declaration
 *     \author  Max Nikonov
 *       \date  07.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef CALC_MODEL_H
#define CALC_MODEL_H

#include BS_FORCE_PLUGIN_IMPORT ()
#include "convert_units.h"
#include "constants.h"
#include "prvd_table.h"
#include "arrays.h"
#include "rocktab_table.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include "fi_params.h"
#include "well_type_helper.h"
#include "norm_calc.h"
#include "shared_vector.h"
#include "jacobian.h"

#include "calc_model_data.h"

namespace blue_sky
  {

  /**
   * \enum  restore_solution_return_type
   * \brief type of return value for restore_solution
   * \todo  describe enum values
   * */
  enum restore_solution_return_type
  {
    SMALL_TIME_STEP_FAIL = -1,      //!< solution restoring failed, newton process should be terminated
    SMALL_TIME_STEP_OK = 0,         //!< solution restored successfully
    SMALL_TIME_STEP_CHOP = 256,     //!< solution restoring failed, newton process should be restarted
  };

  template <typename strategy_t>
  class calc_model;

  template <typename strategy_t>
  class well;

  template <typename strategy_t>
  class reservoir;

  namespace wells
  {
    template <typename strategy_t>
    class connection;
  }

  class well_results_storage;
  class fip_results_storage;
  ///////////////////////////////////////////////////////////////////////////


  /**
   * \class calc_model_data_tmp_holder 
   * \brief holds temporary data 
   *        (shored on new newton iteration, for example)
   * */
  template <typename strategy_t>
  struct BS_API_PLUGIN calc_model_data_tmp_holder
    {
      typedef typename strategy_t::item_t         item_t;         //!< item type (floating point)
      typedef typename strategy_t::item_array_t   item_array_t;   //!< type for array of item_t values
      //! type for array of main_var_type values
      typedef typename strategy_t::template vec <main_var_type>::type main_var_array_t;
      typedef calc_model <strategy_t>             calc_model_t;   //!< calc_model type
      typedef smart_ptr <calc_model_t, true>      sp_calc_model_t;//!< smart_ptr to calc_model type

    public:

      /**
       * \brief  stores data from calc_model in holder
       * \param  calc_model pointer to calc_model instance
       * */
      void 
      save (const sp_calc_model_t &calc_model);

      /**
       * \brief  restores data from holder to calc_model
       * \param  calc_model pointer to calc_model instance
       * */
      void 
      restore (sp_calc_model_t &calc_model);

    public:
      item_array_t              pressure;       //!< pressure array
      item_array_t              saturation_3p;  //!< saturation array (3 phase)
      item_array_t              gas_oil_ratio;  //!< gas_oil_ratio array
      main_var_array_t          main_var;       //!< main variables
    };

  ///////////////////////////////////////////////////////////////////////////
  /**
   * \class calc_model
   * \brief holds calculated data and implements some useful functions
   * */
  template <class strategy_t>
  class BS_API_PLUGIN calc_model : public bs_node
    {
    public:

      typedef strategy_t                                strategy_type;            //!< strategy_type

      typedef calc_model<strategy_t>                    this_t;                   //!< shortname for this type
      typedef smart_ptr<this_t, true>                   sp_this_t;                //!< smart_ptr to this_t type

      typedef jacobian<strategy_t>                      jacobian_t;               //!< jacobian type
      typedef smart_ptr <jacobian_t, true>              sp_jacobian_t;            //!< smart_ptr to jacobian type
      typedef jacobian_matrix <strategy_t>              jacobian_matrix_t;        //!< jacobian_matrix type
      typedef smart_ptr <jacobian_matrix_t, true>       sp_jacobian_matrix_t;     //!< smart_ptr to jacobian_matrix type

      typedef idata                                     idata_t;                  //!< idata type
      typedef smart_ptr<idata_t, true>                  sp_idata_t;               //!< smart_ptr to idata type

      typedef rs_mesh_iface <strategy_t>                mesh_iface_t;             //!< rs_mesh_iface type
      typedef smart_ptr <mesh_iface_t, true>            sp_mesh_iface_t;          //!< smart_ptr to rs_mesh_iface
      
      typedef rs_smesh_iface <strategy_t>               smesh_iface_t;            //!< rs_smesh_iface type
      typedef smart_ptr <smesh_iface_t, true>           sp_smesh_iface_t;         //!< smart_ptr to rs_smesh_iface

      typedef reservoir <strategy_t>                    reservoir_t;              //!< reservoir type
      typedef smart_ptr <reservoir_t, true>             sp_reservoir_t;           //!< smart_ptr to reservoir type

      typedef typename strategy_t::item_array_t         item_array_t;             //!< type for array of item_t values
      typedef typename strategy_t::index_array_t        index_array_t;            //!< type for array of index_t values
      typedef typename strategy_t::item_t               item_t;                   //!< item value (floating point)
      typedef typename strategy_t::index_t              index_t;                  //!< index value (integral type)
      
      typedef typename strategy_t::csr_matrix_t         csr_matrix_t;             //!< shortname for csr_matrix type
      typedef smart_ptr <csr_matrix_t, true>            sp_csr_matrix_t;          //!< smart_ptr to csr_matrix type

      typedef calc_model_data <strategy_t>              data_t;                   //!< calc_model data, each instance for one mesh cell
      typedef typename strategy_t::template vec <data_t>::type data_array_t;      //!< array of calc_model_data values, each value for one mesh cell

      typedef scal_3p<strategy_t>                       scal_3p_t;                //!< scal_3p type
      typedef scale_array_holder <strategy_t>           scale_array_holder_t;     //!< type of holder of scale arrays

      typedef wells::connection <strategy_t>            connection_t;             //!< base type for well connections
      typedef well <strategy_t>                         well_t;                   //!< base type for wells

      typedef smart_ptr< scal_3p_t, true>               sp_scal3p;                //!< smart_ptr to scal_3p type
      typedef smart_ptr <scale_array_holder_t, true>    sp_scale_array_holder_t;  //!< smart_ptr to scale_array_holder type

      typedef smart_ptr< rock_grid< strategy_t >, true> sp_rock_grid;             //!< smart_ptr to rock_grid type
      typedef smart_ptr< fi_params, true>               sp_fi_params;             //!< smart_ptr to fi_params type
      typedef smart_ptr <connection_t, true>            sp_connection_t;          //!< smart_ptr to connection type
      typedef smart_ptr <well_t, true>                  sp_well_t;                //!< smart_ptr to well type

      typedef smart_ptr<well_results_storage, true>     sp_well_results_storage;  //!< smart_ptr to well_results_storage type
      typedef smart_ptr<fip_results_storage, true>      sp_fip_results_storage;   //!< smart_ptr to fip_results_storage type

      typedef pvt_base< strategy_t >                    pvt_base_t;               //!< type of base pvt class
      typedef pvt_oil < strategy_t >                    pvt_oil_t;                //!< pvt_oil type
      typedef pvt_dead_oil< strategy_t >                pvt_dead_oil_t;           //!< pvt_dead_oil type
      typedef pvt_gas< strategy_t >                     pvt_gas_t;                //!< pvt_gas type
      typedef pvt_water< strategy_t >                   pvt_water_t;              //!< pvt_water type

      typedef smart_ptr <pvt_base_t, true>              sp_pvt_t;                 //!< smart_ptr to pvt_base type
      typedef smart_ptr <pvt_dead_oil_t, true>          sp_pvt_oil;               //!< smart_ptr to pvt_dead_oil type
      typedef smart_ptr <pvt_dead_oil_t, true>          sp_pvt_dead_oil;          //!< smart_ptr to pvt_dead_oil type
      typedef smart_ptr <pvt_gas_t, true>               sp_pvt_gas;               //!< smart_ptr to pvt_gas type
      typedef smart_ptr <pvt_water_t, true>             sp_pvt_water;             //!< smart_ptr to pvt_water type

      typedef std::vector< sp_pvt_t >                   sp_pvt_array_t;           //!< type for array of pvt_base objects
      typedef std::vector< sp_pvt_oil >                 sp_pvt_oil_array_t;       //!< type for array of pvt_dead_oil objects
      typedef std::vector< sp_pvt_dead_oil >            sp_pvt_dead_oil_array_t;  //!< type for array of pvt_dead_oil objects
      typedef std::vector< sp_pvt_gas >                 sp_pvt_gas_array_t;       //!< type for array of pvt_gas objects
      typedef std::vector< sp_pvt_water >               sp_pvt_water_array_t;     //!< type for array of pvt_water objects

      typedef std::vector< int >                        vec_i;                    //!< vector of int values

      typedef boost::array <index_t, FI_PHASE_TOT>      phase_d_t;                //!< type for array of shifts for each phase
      typedef boost::array <index_t, FI_PHASE_TOT>      sat_d_t;                  //!< type for array of shifts for each phase, for saturation
      typedef norms_storage <strategy_t>                norms_storage_t;          //!< norms_storage type
      typedef boost::array <item_t, FI_PHASE_TOT>       invers_fvf_avgerage_t;    //!< type for store invers_fvf_average value (by phase)

      //! type for array of main_var_type values
      typedef typename strategy_t::template vec <main_var_type>::type  main_var_array_t;

      //! temporary holder for calc_model data
      typedef calc_model_data_tmp_holder <strategy_t>   calc_model_data_tmp_holder_t; 

    public:

      /**
       * \brief  calc_model dtor
       * */
      ~calc_model();

      /**
       * \brief  inits calc_model
       * */
      void 
      init();

      /**
       * \brief  assignment operator
       * \param  src reference to calc_model instance
       * \return this
       * \todo   should be removed
       * */
      const this_t &
      operator=(const this_t &src);

      /**
       * \brief  inits main arrays
       * \param  input_data pointer to idata instance
       * \param  mesh pointer to mesh instance
       * \return 0 on success otherwise negative integer value
       * \todo   remove return values, throw exceptions instead
       * */
      int 
      init_main_arrays (const sp_idata_t &input_data, const sp_mesh_iface_t &mesh);

      /**
       * \brief  inits arrays for calculation process
       * \param  input_data pointer to idata instance
       * \param  mesh pointer to mesh instance
       * \return 0 on success otherwise negative integer value
       * \todo   remove return values, throw exceptions instead
       * */
      int 
      init_calcul_arrays (const sp_idata_t &input_data, const sp_mesh_iface_t &mesh);

      /**
       * \brief  inits initial conditions
       * \param  input_data pointer to idata instance
       * \param  mesh pointer to mesh instance
       * \return 0 on success otherwise negative integer value
       * \todo   remove return values, throw exceptions instead
       * */
      int 
      set_initial_data (const sp_idata_t &input_data, const sp_mesh_iface_t &mesh);

      /**
       * \brief  calculates equilibrium
       * \param  input_data pointer to idata instance
       * \param  mesh pointer to mesh instance
       * \return 0 on success otherwise negative integer value
       * \todo   remove return values, throw exceptions instead
       * */
      int 
      calc_equil (const sp_idata_t &input_data, const sp_mesh_iface_t &mesh);

      /**
       * \brief  calculates pressure for equil
       * \return 0 on success otherwise negative integer value
       * \todo   remove return values, throw exceptions instead
       * \param[in]   prev_press preview pressure value
       * \param[in]   cur_d
       * \param[in]   h
       * \param[in]   phase
       * \param[in]   i_pvt
       * \param[in]   rs_type
       * \param[in]   depth_goc
       * \param[in]   rs_dat
       * \param[in]   rsvd
       * \param[in]   pbvd
       * \param[out]  p
       * \param[out]  rs
       * */
      int 
      equil_calc_pressure (item_t prev_press, item_t cur_d, item_t h, index_t phase, index_t i_pvt,
                           double rs_type, item_t depth_goc, item_t rs_dat,
                           val_vs_depth *rsvd, val_vs_depth *pbvd,
                           item_t &p, item_t *rs = 0);

      /**
       * \brief  inits pressure array
       * \param  input_data pointer to idata instance
       * \param  mesh pointer to mesh instance
       * \return 0 on success otherwise negative integer value
       * \todo   remove return values, throw exceptions instead
       * */
      int 
      init_pressure (const sp_idata_t &input_data, const sp_mesh_iface_t &mesh);

      /**
       * \brief  inits saturation array
       * \param  input_data pointer to idata instance
       * \param  mesh pointer to mesh instance
       * \return 0 on success otherwise negative integer value
       * \todo   remove return values, throw exceptions instead
       * */
      int 
      init_saturation (const sp_idata_t &input_data, const sp_mesh_iface_t &mesh);

      /**
       * \brief  inits phases variables (rs) and selest number of phases
       * \param  input_data pointer to idata instance
       * \param  mesh pointer to mesh instance
       * \return 0 on success otherwise negative integer value
       * \todo   remove return values, throw exceptions instead
       * */
      int 
      init_rs (const sp_idata_t &input_data, const sp_mesh_iface_t &mesh);


      /**
       * \brief  inits scal 
       * */
      void 
      init_scal ();

      /**
       * \brief  inits pvt arrays
       * \param  pvto_ array of pvt_oil objects
       * \param  pvtg_ array of pvt_gas objects
       * \param  pvtw_ array of pvt_water objects
       * \param  idata_ pointer to idata instance that 
       *                holds input data for pvt
       * */
      void 
      init_pvt_arrays (sp_pvt_oil_array_t &pvto_,
                       sp_pvt_gas_array_t &pvtg_,
                       sp_pvt_water_array_t &pvtw_,
                       const sp_idata_t &idata_);

      /**
       * \brief  returns true if water phase present in model
       * \param  phases phases enum
       * \return true if water phase present in model
       * \todo   obsolete, should be removed
       * */
      static bool 
      is_water_phase (int phases)
      {
        return FI_CHK_WATER (phases);
      }
      /**
       * \brief  returns true if gas phase present in model
       * \param  phases phases enum
       * \return true if gas phase present in model
       * \todo   obsolete, should be removed
       * */
      static bool 
      is_gas_phase (int phases)
      {
        return FI_CHK_GAS (phases);
      }
      /**
       * \brief  returns true if oil phase present in model
       * \param  phases phases enum
       * \return true if oil phase present in model
       * \todo   obsolete, should be removed
       * */
      static bool 
      is_oil_phase (int phases)
      {
        return FI_CHK_OIL (phases);
      }

      /**
       * \brief  return calc_model_data for n_block cell
       * \param  n_block index of cell in mesh
       * \return calc_model_data for n_block cell
       * */
      const data_t &
      get_data (index_t n_block) const
        {
          return data[n_block];
        }

      /**
       * \todo obsolete, should be removed  
       * */
      item_t  
      get_initial_rho (item_t height) const;

      /**
       * \todo obsolete, should be removed  
       * */
      void    
      update_min_pressure_range (item_t min_range);

      /**
       * \todo obsolete, should be removed  
       * */
      void    
      update_max_pressure_range (item_t max_range);

      /**
       * \brief  calculates fluid volume for previous step
       * \param  istart is first call on current step
       * \param  mesh pointer to mesh instance
       * */
      void    
      calc_prev_fluid_volume (bool istart, const sp_mesh_iface_t &mesh);

      /**
       * \brief  restores solution from jacobian_matrix
       * \param  mesh pointer to mesh instance
       * \param  jacobian_mx pointer to jacobian_matrix instance
       * \return see restore_solution_return_type for more details
       * */
      restore_solution_return_type
      restore_solution (const sp_mesh_iface_t &mesh, const sp_jacobian_matrix_t &jacobian_mx);

      /**
       * \brief  applies newton correction and multiplies it by given mult (default 1.0)
       * \param  mult multiplier
       * \param  istart_line_search
       * \param  mesh pointer to mesh instance
       * \param  jacobian_mx pointer to jacobian_matrix instance
       * \return see restore_solution_return_type for more details
       * */
      restore_solution_return_type
      apply_newton_correction (item_t mult, index_t istart_line_search, const sp_mesh_iface_t &mesh, const sp_jacobian_matrix_t &jacobian_mx);

      /**
       * \brief  calculates multiplier (m) for newton correction vector (J(x0) * w = -F(x0), x1 = x0 + m * w)
       * \param  mesh pointer to mesh instance
       * \param  jacobian_mx pointer to jacobian_matrix instance
       * \return multiplier value
       * */
      item_t  
      new_simple_get_cell_solution_mult (const sp_mesh_iface_t &mesh, const sp_jacobian_matrix_t &jacobian_mx);

      /**
       * \brief  returns multiplier for apply_newton_correction function
       * \param  mesh pointer to mesh instance
       * \param  jacobian_mx pointer to jacobian_matrix instance
       * \return multiplier value
       * */
      item_t  
      new_simple_get_cell_solution_mult_2 (const sp_mesh_iface_t &mesh, const sp_jacobian_matrix_t &jmatrix) const;

      /**
       * \brief  restores solution from Jacobian (jacobian_matrix), 
       *         using mult if istart_linear_search != 0
       * \param  mult x^n = x^(n-1) + mult * dx
       * \param  istart_linear_search use mult if != 0
       * \return 0 on success
       * */
      int     
      new_simple_get_cell_solution (const double mult, int istart_linear_search, const sp_mesh_iface_t &msh, const sp_jacobian_matrix_t &jacobian_mx);

      /**
       * \brief  calculates approximal value of So, Sg, Ro from 
       *         full mass of gas (mg_in) and oil (mo_in)
       * \return 0
       * \param[in]  mo_in
       * \param[in]  mg_in
       * \param[in]  poro
       * \param[in]  ifvf_o
       * \param[in]  ifvf_g
       * \param[in]  max_ro
       * \param[out] so
       * \param[out] sg
       * \param[out] ro
       * \param[out] m_var
       * */
      int     
      calc_approx_so_sg_ro (const item_t mo_in, const item_t mg_in, const item_t poro,
                            const item_t ifvf_o, const item_t ifvf_g, const item_t max_ro,
                            // results
                            item_t &so, item_t &sg, item_t &ro,
                            main_var_type &m_var);

      /**
       * \brief  inits jacobian
       * \param  input_data pointer to idata instance
       * \param  mesh pointer to mesh instance
       * */
      void
      init_jacobian (const sp_jacobian_t &input_data, const sp_mesh_iface_t &mesh);


      /**
       * \brief  returns true if water phase present in model
       * \param  phases phases enum
       * \return true if water phase present in model
       * */
      bool 
      is_water () const;

      /**
       * \brief  returns true if gas phase present in model
       * \param  phases phases enum
       * \return true if gas phase present in model
       * */
      bool 
      is_gas () const;

      /**
       * \brief  returns true if oil phase present in model
       * \param  phases phases enum
       * \return true if oil phase present in model
       * */
      bool 
      is_oil () const;

      /**
       * \brief  returns shift value for water phase 
       * \return phase_d[FI_PHASE_WATER]
       * */
      index_t 
      water_shift () const;

      /**
       * \brief  returns shift value for gas phase 
       * \return phase_d[FI_PHASE_GAS]
       * */
      index_t 
      gas_shift () const;

      /**
       * \brief  returns shift value for oil phase 
       * \return phase_d[FI_PHASE_OIL]
       * */
      index_t 
      oil_shift () const;


      //! blue-sky type declaration
      BLUE_SKY_TYPE_DECL_T (calc_model);

    private:

      /**
       * \todo obsolete, should be removed
       * */
      void
      init_boundary_connections (const sp_idata_t &input_data, const sp_mesh_iface_t &mesh);


    public:
      // Variables
      main_var_array_t                                        main_variable;                  //!< (n_elements) main variables (Sg or Ro) per each grid block
      auto_value <well_model_var_type, WELL_MODEL_3VAR>       well_model_var_;                //!< well model var type
      auto_value <well_model_type, BLACK_OIL>                 well_model_type_;               //!< well model type
      norms_storage_t                                         norm;                           //!< norm storage

      sp_pvt_dead_oil_array_t                                 pvt_oil_array;                  //!< array of pvt_oil objects, length == n_pvt_regions
      sp_pvt_water_array_t                                    pvt_water_array;                //!< array of pvt_water objects, length == n_pvt_regions
      sp_pvt_gas_array_t                                      pvt_gas_array;                  //!< array of pvt_gas objects, length == n_pvt_regions

      data_array_t                                            data;                           //!< array of calc_model_data, length == n_elements 

      item_array_t                                            workspace;                      //!< obsolete, n_elements * (n_phase)

      calc_model_data_tmp_holder_t                            old_data_;                      //!< calc_model data stored on previous step
      calc_model_data_tmp_holder_t                            prev_niter_data_;               //!< calc_model data stored on previous newton iteration

      //! pressure on the boundary
      item_array_t                                            bconn_pressure;                 //!< obsolete
      item_array_t                                            bconn_saturation;               //!< obsolete
      item_array_t                                            bconn_gor;                      //!< obsolete
      index_array_t                                           bconn_mainvar;                  //!< obsolete

      auto_value <item_t>                                     linear_search_mult;
      auto_value <bool>                                       lsearch_force_newton_step;      //!< obsolete, force to make newton step in any case

      auto_value <int>                                        b_w_w;                          //!< obsolete 
      auto_value <int>                                        b_w_g;                          //!< obsolete
      auto_value <int>                                        b_w_p;                          //!< obsolete
      auto_value <int>                                        b_g_w;                          //!< obsolete
      auto_value <int>                                        b_g_g;                          //!< obsolete
      auto_value <int>                                        b_g_p;                          //!< obsolete
      auto_value <int>                                        b_o_w;                          //!< obsolete
      auto_value <int>                                        b_o_g;                          //!< obsolete
      auto_value <int>                                        b_o_p;                          //!< obsolete
      index_array_t                                           iwksp;                          //!< obsolete, n_elements
      auto_value <int>                                        multi_well_in_cell_flag;        //!< obsolete, if != 0 use sorting

      auto_value <double>                                     ave_volume;                     //!< average volume

      index_array_t                                           max_norm_counter;               //!< obsolete

      //! base parameters
    public:

      auto_value <int>                                        n_comps;                        //!< number of components
      auto_value <int>                                        n_phases;                       //!< number of phases
      auto_value <int>                                        n_sec_vars;                     //!< number if secondary variables
      auto_value <FI_PHASE_ENUM>                              phases;                         //!< sizeof (int) bit fields (1 -- phase present, 0 -- do not present)
      auto_value <int>                                        n_HCcomps;                      //!< number of hydrocarbon components (n_comps = n_HCcomps + 1)
      auto_value <int>                                        n_HCphases;                     //!< number of hydrocarbon phases (n_phases = n_HCphases + 1)
      auto_value <int>                                        n_pri;                          //!< number of primary variables
      auto_value <int>                                        n_sec;                          //!< number of secondary variables
      auto_value <int>                                        n_vars;                         //!< full number of variables

      phase_d_t                                               phase_d;                        //!< displacement of phase in arrays
      sat_d_t                                                 sat_d;                          //!< displacement of phase in saturation array

      auto_value <int>                                        n_pvt_regions;                  //!< number of pvt regions
      auto_value <int>                                        n_sat_regions;                  //!< number of sat regions
      auto_value <int>                                        n_fip_regions;                  //!< number of fip regions

      physical_constants                                      internal_constants;             //!< physical constants in internal units

      index_array_t                                           pvt_regions;                    //!< (n_elements) index of PVT table for given block zero base
      index_array_t                                           sat_regions;                    //!< (n_elements) index of SAT table for given block zero base
      index_array_t                                           fip_regions;                    //!< (n_elements) index of FIP region for cell
      index_array_t                                           rock_regions;                   //!< (n_elements) index of ROCK regions

      auto_value <RPO_MODEL_ENUM, RPO_DEFAULT_MODEL>          rpo_model;                      //!< 3-ph oil relative permeability model: flag 0, 1 or 2 (stone model)
      sp_scal3p                                               scal_prop;                      //!< scal properties
      sp_rock_grid                                            rock_grid_prop;                 //!< rock and grid properties

      std::vector<rocktab_table <base_strategy_fi> >          rocktab;                        //!< (rocktab table)


      auto_value <double>                                     last_c_norm;                    //!< obsolete
      auto_value <int>                                        approx_flag;                    //!< flag of initial approximation

      sp_fi_params                                            ts_params;                      //!< input model params
      invers_fvf_avgerage_t                                   invers_fvf_average;             //!< (n_phases) 1. / (formation volume factor) for all phases average

      item_array_t                                            plane_flow_rate;                //!< (n_planes * n_phases) flow rates for all planes on current time step
      item_array_t                                            full_step_plane_flow_rate;      //!< (n_planes * n_phases) total flow rates on time step

      item_array_t                                            pressure;                       //!< pressure (n_elements)
      item_array_t                                            saturation_3p;                  //!< (n_phases * n_elements)
      item_array_t                                            gas_oil_ratio;                  //!< gas_oil_ratio (n_elements)

      sp_csr_matrix_t                                         mat;                            //!< obsolete
      sp_well_results_storage                                 well_res;                       //!< storage for well results
      sp_fip_results_storage                                  fip_res;                        //!< storage for fip results
    };
}

#endif // CALC_MODEL_H
