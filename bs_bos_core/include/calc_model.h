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
#include "array_ext.h"
#include "jacobian.h"

#include "calc_model_data.h"

namespace blue_sky
  {

  enum restore_solution_return_type
  {
    SMALL_TIME_STEP_FAIL = -1,
    SMALL_TIME_STEP_OK = 0,
    SMALL_TIME_STEP_CHOP = 256,     //#define YS_SMALL_TIME_STEP_CHOP                 (256)
  };

  template <typename strategy_t>
  class calc_model;

  template <typename strategy_t>
  struct well_mobility_calc;

  template <typename strategy_t>
  class well;

  template <typename strategy_t>
  class reservoir;

  namespace wells
  {
    template <typename strategy_t>
    class connection;
  }

  /*
  template <typename strategy_t>
  class BS_API_PLUGIN rs_mesh_iface;

  template <typename strategy_t>
  class BS_API_PLUGIN rs_mesh_struct_iface;
  */

  class well_results_storage;
  class fip_results_storage;
  ///////////////////////////////////////////////////////////////////////////


  template <typename strategy_t>
  struct BS_API_PLUGIN calc_model_data_tmp_holder
    {
      typedef typename strategy_t::item_t         item_t;
      typedef typename strategy_t::item_array_t   item_array_t;
      typedef typename strategy_t::template vec <main_var_type>::type main_var_array_t;
      typedef calc_model <strategy_t>             calc_model_t;
      typedef smart_ptr <calc_model_t, true>      sp_calc_model_t;

public:

      void save (const sp_calc_model_t &calc_model);
      void restore (sp_calc_model_t &calc_model);

public:
      item_array_t              pressure;
      item_array_t              saturation_3p;
      item_array_t              gas_oil_ratio;
      main_var_array_t          main_var;
    };

  ///////////////////////////////////////////////////////////////////////////
  template <class strategy_t>
  class BS_API_PLUGIN calc_model : public bs_node
    {
    public:

      typedef strategy_t                                strategy_type;

      typedef calc_model<strategy_t>						        this_t;
      typedef smart_ptr<this_t, true>                   sp_this_t;

      typedef jacobian<strategy_t>                      jacobian_t;
      typedef smart_ptr <jacobian_t, true>              sp_jacobian_t;
      typedef jacobian_matrix <strategy_t>              jacobian_matrix_t;
      typedef smart_ptr <jacobian_matrix_t, true>       sp_jacobian_matrix_t;

      typedef idata                                     idata_t;
      typedef smart_ptr<idata_t, true>                  sp_idata_t;

      typedef rs_mesh_iface <strategy_t> 								mesh_iface_t;
      typedef smart_ptr <mesh_iface_t, true>            sp_mesh_iface_t;
      
      typedef rs_smesh_iface <strategy_t> 			  smesh_iface_t;
      typedef smart_ptr <smesh_iface_t, true>     sp_smesh_iface_t;

      typedef reservoir <strategy_t>						      	reservoir_t;		//!< short name
      typedef smart_ptr <reservoir_t, true>             sp_reservoir_t; //!< short name

      typedef typename strategy_t::item_array_t         item_array_t;
      typedef typename strategy_t::index_array_t        index_array_t;
      typedef typename strategy_t::item_t               item_t;
      typedef typename strategy_t::index_t              index_t;
      
      typedef typename strategy_t::csr_matrix_t         csr_matrix_t;
      typedef smart_ptr <csr_matrix_t, true>            sp_csr_matrix_t;

      typedef calc_model_data <strategy_t>              data_t;
      typedef typename strategy_t::template vec <data_t>::type data_array_t;

      typedef scal_3p<strategy_t>                       scal_3p_t;
      typedef scale_array_holder <strategy_t>						scale_array_holder_t;

      typedef wells::connection <strategy_t>            connection_t;
      typedef well <strategy_t>                         well_t;

      typedef smart_ptr< scal_3p_t, true>               sp_scal3p;
      typedef smart_ptr <scale_array_holder_t, true>		sp_scale_array_holder_t;

      typedef smart_ptr< rock_grid< strategy_t >, true> sp_rock_grid;
      typedef smart_ptr< fi_params, true>               sp_fi_params;
      typedef smart_ptr <connection_t, true>            sp_connection_t;
      typedef smart_ptr <well_t, true>                  sp_well_t;

      typedef smart_ptr<well_results_storage, true>     sp_well_results_storage;
      typedef smart_ptr<fip_results_storage, true>      sp_fip_results_storage;

      typedef pvt_base< strategy_t >                    pvt_base_t;
      typedef pvt_oil < strategy_t >										pvt_oil_t;
      typedef pvt_dead_oil< strategy_t >                pvt_dead_oil_t;
      typedef pvt_gas< strategy_t >                     pvt_gas_t;
      typedef pvt_water< strategy_t >                   pvt_water_t;

      typedef smart_ptr <pvt_base_t, true>              sp_pvt_t;
      typedef smart_ptr <pvt_dead_oil_t, true>          sp_pvt_oil;
      typedef smart_ptr <pvt_dead_oil_t, true>          sp_pvt_dead_oil;
      typedef smart_ptr <pvt_gas_t, true>               sp_pvt_gas;
      typedef smart_ptr <pvt_water_t, true>             sp_pvt_water;

      typedef std::vector< sp_pvt_t >                   sp_pvt_array_t;
      typedef std::vector< sp_pvt_oil >                 sp_pvt_oil_array_t;
      typedef std::vector< sp_pvt_dead_oil >            sp_pvt_dead_oil_array_t;
      typedef std::vector< sp_pvt_gas >                 sp_pvt_gas_array_t;
      typedef std::vector< sp_pvt_water >               sp_pvt_water_array_t;

      typedef std::vector< int >                        vec_i;

      typedef boost::array <index_t, FI_PHASE_TOT>			phase_d_t;
      typedef boost::array <index_t, FI_PHASE_TOT>			sat_d_t;

      typedef typename strategy_t::template vec <main_var_type>::type  main_var_array_t;

      typedef norms_storage <strategy_t>                norms_storage_t;

      typedef calc_model_data_tmp_holder <strategy_t>   calc_model_data_tmp_holder_t;
      typedef well_mobility_calc <strategy_t>           well_mobility_calc_t;

      typedef wells::type_helper <strategy_t>           helper_t;

      typedef typename helper_t::item_q_rate_t          item_q_rate_t;
      typedef typename helper_t::item_q_rate_inflow_t   item_q_rate_inflow_t;
      typedef typename helper_t::item_gas_rate_t        item_gas_rate_t;
      typedef typename helper_t::invers_fvf_avgerage_t  invers_fvf_avgerage_t;

      typedef boost::array <index_t, FI_PHASE_TOT>      up_cell_array_t;

    public:

      ~calc_model();

      void init();

      const this_t &operator=(const this_t &src);

      // initialize main arrays
      int init_main_arrays (const sp_idata_t &input_data, const sp_mesh_iface_t &mesh);

      // initialize arrays for calculation process
      int init_calcul_arrays (const sp_idata_t &input_data, const sp_mesh_iface_t &mesh);

      // init initial conditions
      int set_initial_data (const sp_idata_t &input_data, const sp_mesh_iface_t &mesh);
      // calc equilibrium
      int calc_equil (const sp_idata_t &input_data, const sp_mesh_iface_t &mesh);
      // calc pressure for equil
      int equil_calc_pressure (item_t prev_press, item_t cur_d, item_t h, index_t phase, index_t i_pvt,
                               double rs_type, item_t depth_goc, item_t rs_dat,
                               val_vs_depth *rsvd, val_vs_depth *pbvd,
                               item_t &p, item_t *rs = 0);
      // initialize saturation array
      int init_pressure (const sp_idata_t &input_data, const sp_mesh_iface_t &mesh);

      // initialize saturation array
      int init_saturation (const sp_idata_t &input_data, const sp_mesh_iface_t &mesh);

      // initialize phases variables (rs) and select number of phases
      int init_rs (const sp_idata_t &input_data, const sp_mesh_iface_t &mesh);

      //int initialize_gas_oil_fluid_volume (const sp_idata_t &input_data);

      // prepare arrays for rst file
      //int prepare_data_for_rst (const sp_idata_t &input_data, const sp_mesh &mesh);

      // make initialization for gas-oil model (2ph or 3phase) by solving equation of mass
      //int fi_initialization_gas_oil_model (int istart, const sp_idata_t &input_data);

      void init_scale_arrays (const sp_scale_array_holder_t &array_,
                              const sp_idata_t &idata_,
                              ARRAY_NAME socr,
                              ARRAY_NAME scr,
                              ARRAY_NAME su,
                              ARRAY_NAME sl,
                              ARRAY_NAME pc);
      void init_scal ();

      void init_pvt_arrays (sp_pvt_oil_array_t &pvto_,
                            sp_pvt_gas_array_t &pvtg_,
                            sp_pvt_water_array_t &pvtw_,
                            const sp_idata_t &idata_);

      static bool is_water_phase (int phases)
      {
        return FI_CHK_WATER (phases);
      }
      static bool is_gas_phase (int phases)
      {
        return FI_CHK_GAS (phases);
      }
      static bool is_oil_phase (int phases)
      {
        return FI_CHK_OIL (phases);
      }

      const data_t &get_data (index_t n_block) const
        {
          return data[n_block];
        }

      item_t  get_initial_rho (item_t height) const;
      void    update_min_pressure_range (item_t min_range);
      void    update_max_pressure_range (item_t max_range);
      void    calc_prev_fluid_volume (bool istart, const sp_mesh_iface_t &mesh);

      restore_solution_return_type
      restore_solution (const sp_mesh_iface_t &mesh, const sp_jacobian_matrix_t &jacobian_mx);

      restore_solution_return_type
      apply_newton_correction (item_t mult, index_t istart_line_search, const sp_mesh_iface_t &mesh, const sp_jacobian_matrix_t &jacobian_mx);

      item_t  new_simple_get_cell_solution_mult (const sp_mesh_iface_t &mesh, const sp_jacobian_matrix_t &jacobian_mx);
      item_t  new_simple_get_cell_solution_mult_2 (const sp_mesh_iface_t &mesh, const sp_jacobian_matrix_t &jmatrix) const;
      int     new_simple_get_cell_solution (const double mult, int istart_linear_search, const sp_mesh_iface_t &msh, const sp_jacobian_matrix_t &jacobian_mx);
      int     calc_approx_so_sg_ro (const item_t mo_in, const item_t mg_in, const item_t poro,
                                    const item_t ifvf_o, const item_t ifvf_g, const item_t max_ro,
                                    // results
                                    item_t &so, item_t &sg, item_t &ro,
                                    main_var_type &m_var);

      void
      init_jacobian (const sp_jacobian_t &input_data, const sp_mesh_iface_t &mesh);


      bool is_water () const;
      bool is_gas () const;
      bool is_oil () const;

      index_t water_shift () const;
      index_t gas_shift () const;
      index_t oil_shift () const;


      BLUE_SKY_TYPE_DECL_T (calc_model);

    private:

      void
      init_boundary_connections (const sp_idata_t &input_data, const sp_mesh_iface_t &mesh);


    public:
      // Variables
      main_var_array_t                                        main_variable;                  //!< (n_elements) main variables (Sg or Ro) per each grid block
      auto_value <well_model_var_type, WELL_MODEL_3VAR>       well_model_var_;								//!< well model
      auto_value <well_model_type, BLACK_OIL>                 well_model_type_;
      norms_storage_t                                         norm;

      sp_pvt_dead_oil_array_t                                 pvt_oil_array;                  //!< (n_pvt_regions)
      sp_pvt_water_array_t                                    pvt_water_array;
      sp_pvt_gas_array_t                                      pvt_gas_array;

      data_array_t                                            data;                           //!< n_elements size

      // workspace
      item_array_t              workspace;																										//!< n_elements * (n_phase)

      calc_model_data_tmp_holder_t    old_data_;
      calc_model_data_tmp_holder_t    prev_niter_data_;

      //! pressure on the boundary
      item_array_t              bconn_pressure;
      item_array_t              bconn_saturation;
      item_array_t              bconn_gor;
      index_array_t             bconn_mainvar;

      auto_value <item_t>       linear_search_mult;
      auto_value <bool>         lsearch_force_newton_step;		// flag : force to make newton step in any case

      auto_value <int>          b_w_w;
      auto_value <int>          b_w_g;
      auto_value <int>          b_w_p;
      auto_value <int>          b_g_w;
      auto_value <int>          b_g_g;
      auto_value <int>          b_g_p;
      auto_value <int>          b_o_w;
      auto_value <int>          b_o_g;
      auto_value <int>          b_o_p;
      index_array_t             iwksp;												//!< n_elements
      auto_value <int>          multi_well_in_cell_flag;			//!< if != 0 use sorting

      auto_value <double>       ave_volume;

      index_array_t             max_norm_counter;

      //! base parameters
    public:

      auto_value <int>          n_comps;                            //!< number of components
      auto_value <int>          n_phases;                           //!< number of phases
      auto_value <int>          n_sec_vars;                         //!< number if secondary variables
      auto_value <FI_PHASE_ENUM> phases;                            //!< sizeof (int) bit fields (1 -- phase present, 0 -- do not present)
      auto_value <int>          n_HCcomps;                          //!< number of hydrocarbon components (n_comps = n_HCcomps + 1)
      auto_value <int>          n_HCphases;                         //!< number of hydrocarbon phases (n_phases = n_HCphases + 1)
      auto_value <int>          n_pri;                              //!< number of primary variables
      auto_value <int>          n_sec;                              //!< number of secondary variables
      auto_value <int>          n_vars;                             //!< full number of variables

      phase_d_t                 phase_d;											//!< displacement of phase in arrays
      sat_d_t		                sat_d;												//!< displacement of phase in saturation array

      auto_value <int>          n_pvt_regions;                      //!< number of pvt regions
      auto_value <int>          n_sat_regions;                      //!< number of sat regions
      auto_value <int>          n_fip_regions;                      //!< number of fip regions

      physical_constants        internal_constants;				//!< physical constants in internal units

      index_array_t             pvt_regions;              //!< (n_elements) index of PVT table for given block zero base
      index_array_t             sat_regions;              //!< (n_elements) index of SAT table for given block zero base
      index_array_t             fip_regions;              //!< (n_elements) index of FIP region for cell
      index_array_t             rock_regions;             //!< (n_elements) index of ROCK regions

      auto_value <RPO_MODEL_ENUM, RPO_DEFAULT_MODEL>
      rpo_model;                 //!< 3-ph oil relative permeability model: flag 0, 1 or 2 (stone model)
      sp_scal3p                 scal_prop;                 //!< scal properties
      sp_rock_grid              rock_grid_prop;            //!< rock and grid properties

      std::vector<rocktab_table <base_strategy_fi> > rocktab;     //!< (rocktab table)


      auto_value <double>       last_c_norm;
      auto_value <int>          approx_flag;                        // flag of initial approximation

      sp_fi_params              ts_params;
      invers_fvf_avgerage_t     invers_fvf_average;    //!< (n_phases) 1. / (formation volume factor) for all phases average

      item_array_t              plane_flow_rate;                  //!< (n_planes * n_phases) flow rates for all planes on current time step
      item_array_t              full_step_plane_flow_rate;        //!< (n_planes * n_phases) total flow rates on time step

      item_array_t              pressure;
      item_array_t              saturation_3p;                    //!< (n_phases * n_elements)
      item_array_t              gas_oil_ratio;

      sp_csr_matrix_t           mat;
      sp_well_results_storage   well_res;
      sp_fip_results_storage    fip_res;
    };
}

#endif // CALC_MODEL_H
