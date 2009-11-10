/*! \file fi_params.h
		\brief Comtains fi_params - named_pbase inherited class
		\author Nikonov Max
*/

#ifndef FI_PARAMS_H
#define FI_PARAMS_H

namespace blue_sky
  {
  class BS_API_PLUGIN fi_params : public named_pbase
    {
    public:
      //! blue-sky class declaration
      BLUE_SKY_TYPE_DECL(fi_params);
    public:
      //! enumeration definitions
      PROP_BASE_IDX_DECL_BEGIN(fi_params,named_pbase)
      //! indexes for floating point parameters
      FI_PARAMS_START,
      PVT_PRESSURE_RANGE_MAX,                     //!< PVT tables maximum pressure
      PVT_PRESSURE_RANGE_MIN,                     //!< PVT tables minimal pressure
      PVTO_RS_RANGE_MAX,                          //!< PVTO maximum RS for slop table
      PVTO_RS_RANGE_MIN,                          //!< PVTO minimum RS for slop table
      MAX_TS,                                     //!< maximum allowed time step length
      MIN_TS,                                     //!< minimum allowed time step length
      FIRST_TS,                                   //!< first simulation time step length
      LIN_SOLV_RESIDUAL,                          //!< tolerance for linear solver
      NEWTON_RESIDUAL,                            //!< tolerance for newton process
      TS_INC_MULT,                                //!< multiplier for incrementing time step length
      TS_DEC_MULT,                                //!< multiplier for decrementing time step length
      OVERDRAFT,                                  //!< overdraft factor (any time step could be multiplied by this factor to achieve end of report time step)
      MAX_P_CORRECTION,                           //!< maximum allowed newton correction for pressure
      MAX_S_CORRECTION,                           //!< maximum allowed newton correction for saturation
      MAX_RS_CORRECTION,                          //!< maximum allowed newton correction for gas oil ratio
      MAX_P_WELL_CORRECTION,                      //!< maximum allowed newton correction for well pressure
      WAT_ROW_MULT,                               //!< multiplier for water equation in Jacobian
      GAS_ROW_MULT,                               //!< multiplier for gas equation in Jacobian
      OIL_ROW_MULT,                               //!< multiplier for oil equation in Jacobian
      PRESSURE_COL_MULT,                          //!< multiplier for pressure derivates column in Jacobian
      TS_OMEGA,                                   //!< omega coef in time step controling
      TS_DP,                                      //!< pressure change for time step controling
      TS_DS,                                      //!< saturation change for time step controling
      DP_MIN_CHOP,                                //!< minimum pressure chop
      DS_MIN_CHOP,                                //!< minimum saturation chop
      MAX_ALLOWED_LIN_SOLV_RESID,                 //!< maximum allowed residual
      AMG_RESID,                                  //!< tolerance for AMG preconditioner
      TS_DRS,                                     //!< Rs change for time step controling
      MAX_NORM_ON_TS,                             //!< if norm on time step is greater than this value restart occur
      DRSDT,                                      //!< Maximum rate of increase of solution GOR
      GAS_NORM_MULT,                              //!< multiplier for gas norm
      // compositional
      //COMP_PARAMS_D_MASS_BALANCE,                   //!< Scaled mass balance equation residuals. Equations are scaled by the accumulation terms (component mass dividing timestep size).
      //COMP_PARAMS_D_PHASE_EQUIL,                    //!< Scaled phase equilibrium relation residuals. These are scaled by the component fugacities in the gas phase
      //COMP_PARAMS_D_WELL_RATE,                      //!< Scaled well equation residuals, scaled by well input rate
      //COMP_PARAMS_D_MAX_DP,                         //!< The maximum scaled pressure change (absolute pressure change divided by the average reservoir pressure)
      //COMP_PARAMS_D_MAX_DS,                         //!< The maximum absolute saturation change
      //COMP_PARAMS_D_MAX_DXCP,                       //!< The maximum absolute component mole fraction change
      //COMP_PARAMS_D_TS_DP,                          //!< time step controling pressure change
      //COMP_PARAMS_D_TS_DS,                          //!< time step controling saturation change
      //COMP_PARAMS_D_TS_DXCP,                        //!< time step controling mole fraction change
      //COMP_PARAMS_D_MAX_P_CORRECTION,               //!< maximum allowed newton correction for pressure
      //COMP_PARAMS_D_MAX_S_CORRECTION,               //!< maximum allowed newton correction for saturation
      //COMP_PARAMS_D_MAX_XCP_CORRECTION,             //!< maximum allowed newton correction for mole fraction
      // black oil
      MASS_BALANS_ERROR,                          //!< maximum allowed mass balans error
      MAX_NORM_ON_FIRST_N,                        //!< maximum allowed norm on first newton iteration
      P_DIMENSION_LESS_SCALE_FACTOR,              //!< scale factor for dimension less pressure
      APPL_CHOP_PERC,                             //!< typically (0.25-0.4)
      APPL_CHOP_BND,                              //!< typically (0.1-0.25)
      LIN_SOLV_MATBAL_RESIDUAL,                   //!< mat. balans tolerance for linear solver
      S_RHS_NORM,                                 //!< max value for secondary rhs
      D_TOTAL,                                    //!< all parameters should be before this item

      //! indexes for integer parameters
      PVT_INTERP_POINTS,                          //!< number of interpolation points
      NEWTON_ITERS_NUM,                           //!< maximum allowed newton iterations number
      LIN_ITERS_NUM,                              //!< maximum allowed linear iterations number
      NEWTON_ITERS_INC,                           //!< number of newton iterations to increment step length
      NEWTON_ITERS_AMG,                           //!< number of newton iterations to build full amg setup
      APPROX_STEPS,                               //!< number of approximation steps
      SELECT_SOL_STEPS,                           //!< number of steps to select newton correction force term
      LIN_SOLVER_TYPE,                            //!< type of linear solver
      GMRES_ORTONORM_VLEN,                        //!< number of vectors in GMRES to ortonorm
      AMG_LIN_ITERS_NUM,                          //!< maximum allowed AMG solver iterations
      WELL_NEWTON_ITERS_NUM,                      //!< maximum number of well internal newton iterations
      PREC_TYPE,                                  //!< type of the preconditioner for linear solver
      PREC_TYPE_ONE_PHASE,                        //!< type of the preconditioner for linear solver for 1-phase problem
      MIN_CELLS_IN_REGION,                        //!< minimum allowed cells in region for ACTNUM = 1
      SAVE_STEP_DATA_PERIOD,                      //!< save step data every n step
      I_TOTAL,                                      //!< all parameters should be before this item

      //! indexes for boolean parameters
      PRINT_PVTO_TABLE,                           //!< if this item is set print PVTO to LOG file
      PRINT_PVTW_TABLE,                           //!< print pvtw interpolation table
      PRINT_PVTG_TABLE,                           //!< print pvtg interpolation table
      STORE_PANE_FLOW_RATES,
      NEW_TS_SELECTION_ALGORITHM,                 //!< use new algorithm for ts selection
      SAVE_STEP_MAPS,                             //!< save maps on time step in the rst file
      SAVE_INITIAL_DATA,                          //!< save initial data in the rst file
      SAVE_RESTART_DATA,                          //!< save restart information
      NEWTRANS,                                   //!< correct transmissibility calculation between cells
      CHECKS_FOR_MULTI_CONN,                      //!< check for multi connections in one cell
      USE_TIMESTEP_FILE,                          //!< save time steps into *tsl file
      SIMPLE_GET_CELL_SOLUTION,                   //!< enable simple solution restore algorithm
      USE_LOW_SKIN_TRANS_MULT,                    //!< if this item is set use calculation of trans mults for low skin
      SAVE_MESH_DATA,                             //!< save initial data in the rst file
      SAVE_CALC_DATA,                             //!< save solution vector at each large time step into .dmp file
      LOAD_CALC_DATA,                             //!< load initial guesses for Newton iterations from .dmp file
      SAVE_NORM_DATA,                             //!< save cells norms to .cel file
      DENS_JAC_SCALE,                             //!< enable densities based multiplier's for gas and oil equation's in Jacobian
      NEWTON_TUNING,
      SAVE_CROSSFLOW,
      G_FRACTURE,
      CLAMP_PRESSURE,                             //!< clamp pressure values at each Newton iteration between possible min and max
      G_FRACTURE_FINAL,
      SAVE_CONN_DATA,                             //!< save connection rates to rst file
      DISABLE_FRACTURE_CHECK,
      COMPRESS_IRR,                               //!< compress irregular matrix
      REMOVE_ISOLATED_REGIONS,                    //!< make isolated regions with out wells inactive
      SAVE_WRATES_TO_ASCII_FILE,                  //!< save well rates to ascii file
      CREATE_HDF5_FILE,                           //!< create hdf5 file
      WRITE_PRESSURE_TO_HDF5,                     //!< write pressure to hdf5
      WRITE_SATURATION_TO_HDF5,                   //!< write saturation to hdf5
      FIX_SOIL_BUG,
      DISABLE_CROSSFLOW,
      WRITE_GAS_OIL_RATIO_TO_HDF5,                //!< write gas-oil ratio values to hdf5 file
      WRITE_FIP_RESULTS_TO_HDF5,                  //!< write data by fip regions to hdf5 file
      WRITE_WELL_RESULTS_TO_HDF5,                 //!< write well data to hdf5 file
      WRITE_CONN_RESULTS_TO_HDF5,                 //!< write connection data to hdf5 file
      WRITE_PLANE_FLOW_RATES_TO_HDF5,             //!< write plane flow rates to hdf5 file
      USE_CFL,                                    //!< compute CFL in fi_operator and use csr_ilu_cfl_prec, used only if BS_BOS_CORE_USE_CSR_ILU_CFL_PREC defined
      WRITE_CFL_TO_HDF5,                          //!< write CFL vector to hdf5 file
      B_TOTAL,                                      //!< all parameters should be before this item
      FI_PARAMS_TOTAL,
      PROP_BASE_IDX_DECL_END

      //! virtual dtor
      virtual ~fi_params ();

      //! method for access property table
      PBASE_ACCESS_MS(fi_params)

    public:
      typedef smart_ptr<FRead, true> sp_reader_t;

      //! Assignment operator
      const fi_params &operator=(const fi_params&);
      //! Update operator
      const fi_params &operator+=(const fi_params&);

      //! set default values for properties
      void set_default_values ();

      // read params from keyword data file
      void read_from_keyword_data_file (sp_reader_t r);

      //! write data to result file
      //int write_data_to_res_file (res_frame_ctrl *parent_frame, int id = YS_ID_FI_PARAMS);

      //! read data from result file
      //int read_data_from_res_file (res_frame_ctrl *parent_frame);

      //! returns name by index
      virtual const std::string &get_params_name (idx_type idx);
    };

#ifdef BSPY_EXPORTING_PLUGIN
  void
  inline py_export_fi_params ()
  {
    using namespace boost::python;
    class_ <fi_params, bases <named_pbase>, boost::noncopyable> ("fi_params", no_init)
      ;

    register_ptr_to_python <smart_ptr <fi_params, true> > ();
    implicitly_convertible <smart_ptr <fi_params, true>, smart_ptr <named_pbase, true> > ();
  }
#endif

}

#endif // FI_PARAMS_H
