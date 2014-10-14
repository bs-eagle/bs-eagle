/**
 *       \file  fi_params.cpp
 *      \brief  fi_params class implementation
 *     \author  Nikonov Max
 *       \date  27.06.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "fi_params.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "bos_reader_iface.h"
#include "convert_units.h"
#include "constants.h"
#include "err_num_def.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {

  void fi_params::set_default_values ()
  {
    physical_constants ph_const;

    this->params_names.resize(fi_params::FI_PARAMS_TOTAL+1);
    ereg <fi_params> (PVT_PRESSURE_RANGE_MAX, "PVT_PRESSURE_MAX", "PVT tables maximum pressure");
    ereg <fi_params> (PVT_PRESSURE_RANGE_MIN, "PVT_PRESSURE_MIN", "PVT tables minimal pressure");
    ereg <fi_params> (PVTO_RS_RANGE_MAX,      "PVTO_RS_MAX", "PVTO maximum RS for slop table");
    ereg <fi_params> (PVTO_RS_RANGE_MIN,      "PVTO_RS_MIN", "PVTO minimum RS for slop table");
    ereg <fi_params> (MAX_TS,                 "TS_MAX", "maximum allowed time step length");
    ereg <fi_params> (MIN_TS,                 "TS_MIN", "minimum allowed time step length");
    ereg <fi_params> (FIRST_TS,               "TS_FIRST", "first simulation time step length");
    ereg <fi_params> (LIN_SOLV_RESIDUAL,      "LITERS_RESID", "tolerance for linear solver");
    ereg <fi_params> (NEWTON_RESIDUAL,        "NITERS_RESID", "tolerance for newton process");
    ereg <fi_params> (TS_INC_MULT,            "TS_INC_MULT", "multiplier for incrementing time step length");
    ereg <fi_params> (TS_DEC_MULT,            "TS_DEC_MULT", "multiplier for decrementing time step length");
    ereg <fi_params> (OVERDRAFT,              "TS_OVERDRAFT", "overdraft factor (any time step could be multiplied by this factor to achieve end of report time step)");
    ereg <fi_params> (MAX_P_CORRECTION,       "P_CORR_MAX", "maximum allowed newton correction for pressure");
    ereg <fi_params> (MAX_S_CORRECTION,       "S_CORR_MAX", "maximum allowed newton correction for saturation");
    ereg <fi_params> (MAX_RS_CORRECTION,      "RS_CORR_MAX", "maximum allowed newton correction for gas oil ratio");
    ereg <fi_params> (MAX_P_WELL_CORRECTION,  "WELL_P_CORR__MAX", "maximum allowed newton correction for well pressure");
    EREG(fi_params,WAT_ROW_MULT,"multiplier for water equation in Jacobian");
    EREG(fi_params,GAS_ROW_MULT,"multiplier for gas equation in Jacobian");
    EREG(fi_params,OIL_ROW_MULT,"multiplier for oil equation in Jacobian");
    ereg <fi_params> (PRESSURE_COL_MULT,      "P_COL_MULT", "multiplier for pressure derivates column in Jacobian");
    EREG(fi_params,TS_OMEGA,"omega coef in time step controling");
    EREG(fi_params,TS_DP,"pressure change for time step controling");
    EREG(fi_params,TS_DS,"saturation change for time step controling");
    EREG(fi_params,DP_MIN_CHOP,"minimum pressure chop");
    EREG(fi_params,DS_MIN_CHOP,"minimum saturation chop");
    ereg <fi_params> (MAX_ALLOWED_LIN_SOLV_RESID, "LITERS_MAX_RESID", "maximum allowed residual");
    EREG(fi_params,AMG_RESID,"tolerance for AMG preconditioner");
    EREG(fi_params,TS_DRS,"Rs change for time step controling");
    EREG(fi_params,MAX_NORM_ON_TS,"if norm on time step is greater than this value restart occur");
    EREG(fi_params,DRSDT,"Maximum rate of increase of solution GOR");
    EREG(fi_params,GAS_NORM_MULT,"multiplier for gas norm");
    // compositional
    /*EREG(fi_params,COMP_PARAMS_D_MASS_BALANCE,"Scaled mass balance equation residuals. Equations are scaled by the accumulation terms (component mass dividing timestep size).");
    EREG(fi_params,COMP_PARAMS_D_PHASE_EQUIL,"Scaled phase equilibrium relation residuals. These are scaled by the component fugacities in the gas phase");
    EREG(fi_params,COMP_PARAMS_D_WELL_RATE,"Scaled well equation residuals, scaled by well input rate");
    EREG(fi_params,COMP_PARAMS_D_MAX_DP,"The maximum scaled pressure change (absolute pressure change divided by the average reservoir pressure)");
    EREG(fi_params,COMP_PARAMS_D_MAX_DS,"The maximum absolute saturation change");
    EREG(fi_params,COMP_PARAMS_D_MAX_DXCP,"The maximum absolute component mole fraction change");
    EREG(fi_params,COMP_PARAMS_D_TS_DP,"time step controling pressure change");
    EREG(fi_params,COMP_PARAMS_D_TS_DS,"time step controling saturation change");
    EREG(fi_params,COMP_PARAMS_D_TS_DXCP,"time step controling mole fraction change");
    EREG(fi_params,COMP_PARAMS_D_MAX_P_CORRECTION,"maximum allowed newton correction for pressure");
    EREG(fi_params,COMP_PARAMS_D_MAX_S_CORRECTION,"maximum allowed newton correction for saturation");
    EREG(fi_params,COMP_PARAMS_D_MAX_XCP_CORRECTION,"maximum allowed newton correction for mole fraction");*/
    // black oil
    EREG(fi_params,MASS_BALANS_ERROR,"maximum allowed mass balans error");
    EREG(fi_params,MAX_NORM_ON_FIRST_N,"maximum allowed norm on first newton iteration");
    EREG(fi_params,P_DIMENSION_LESS_SCALE_FACTOR,"scale factor for dimension less pressure");
    EREG(fi_params,APPL_CHOP_PERC,"APPL_CHOP_PERC");
    EREG(fi_params,APPL_CHOP_BND,"APPL_CHOP_BND");
    EREG(fi_params,LIN_SOLV_MATBAL_RESIDUAL,"mat. balans tolerance for linear solver");
    EREG(fi_params,D_TOTAL,"all parameters should be before this item");

    //! indexes for integer parameters
    EREG(fi_params,PVT_INTERP_POINTS,"number of interpolation points");
    ereg <fi_params> (NEWTON_ITERS_NUM,       "NITERS_NUM", "maximum allowed newton iterations number");
    ereg <fi_params> (LIN_ITERS_NUM,          "LITERS_NUM", "maximum allowed linear iterations number");
    ereg <fi_params> (NEWTON_ITERS_INC,       "NITERS_INC", "number of newton iterations to increment step length");
    ereg <fi_params> (NEWTON_ITERS_AMG,       "NITERS_AMG", "number of newton iterations to build full amg setup");
    EREG(fi_params,APPROX_STEPS,"number of approximation steps");
    EREG(fi_params,SELECT_SOL_STEPS,"number of steps to select newton correction force term");
    ereg <fi_params> (LIN_SOLVER_TYPE,        "LSOLV_TYPE", "type of linear solver");
    EREG(fi_params,GMRES_ORTONORM_VLEN,"number of vectors in GMRES to ortonorm");
    ereg <fi_params> (AMG_LIN_ITERS_NUM,      "AMG_LITERS_NUM", "maximum allowed AMG solver iterations");
    ereg <fi_params> (WELL_NEWTON_ITERS_NUM,  "WELL_NITERS_NUM", "maximum number of well internal newton iterations");
    EREG(fi_params,PREC_TYPE,"type of the preconditioner for linear solver");
    ereg <fi_params> (PREC_TYPE_ONE_PHASE,    "1PHASE_PREC_TYPE", "type of the preconditioner for linear solver for 1 phase problem");
    EREG(fi_params,MIN_CELLS_IN_REGION,"minimum allowed cells in region for ACTNUM = 1");
    EREG(fi_params,SAVE_STEP_DATA_PERIOD,"save step data every n step");
    EREG(fi_params,I_TOTAL,"all parameters should be before this item");

    //! indexes for boolian parameters
    EREG(fi_params,PRINT_PVTO_TABLE,"if this item is set print PVTO to LOG file");
    EREG(fi_params,PRINT_PVTW_TABLE,"print pvtw interpolation table");
    EREG(fi_params,PRINT_PVTG_TABLE,"print pvtg interpolation table");
    EREG(fi_params,STORE_PANE_FLOW_RATES,"boolean store pane flow rates");
    ereg <fi_params> (NEW_TS_SELECTION_ALGORITHM, "NEW_TS_SELECTION", "use new algorithm for ts selection");
    EREG(fi_params,SAVE_STEP_MAPS,"save maps on time step in the rst file");
    EREG(fi_params,SAVE_INITIAL_DATA,"save initial data in the rst file");
    EREG(fi_params,SAVE_RESTART_DATA,"save restart information");
    EREG(fi_params,NEWTRANS,"correct transmissibility calculation between cells");
    ereg <fi_params> (CHECKS_FOR_MULTI_CONN,      "CHECK_FOR_MULTY_CONN", "check for multi connections in one cell");
    EREG(fi_params,USE_TIMESTEP_FILE,"save time steps into *tsl file");
    EREG(fi_params,SIMPLE_GET_CELL_SOLUTION,"enable simple solution restore algorithm");
    EREG(fi_params,USE_LOW_SKIN_TRANS_MULT,"if this item is set use calculation of trans mults for low skin");
    EREG(fi_params,SAVE_MESH_DATA,"save initial data in the rst file");
    EREG(fi_params,SAVE_CALC_DATA,"save solution vector at each large time step into .dmp file");
    EREG(fi_params,LOAD_CALC_DATA,"load initial guesses for Newton iterations from .dmp file");
    EREG(fi_params,SAVE_NORM_DATA,"save cells norms to .cel file");
    EREG(fi_params,DENS_JAC_SCALE,"enable densities based multiplier's for gas and oil equation's in Jacobian ");
    EREG(fi_params,NEWTON_TUNING,"newton tuning");
    EREG(fi_params,SAVE_CROSSFLOW,"save crossflow");
    EREG(fi_params,G_FRACTURE,"g fracture");
    EREG(fi_params,CLAMP_PRESSURE,"clamp pressure values at each Newton iteration between possible min and max");
    EREG(fi_params,G_FRACTURE_FINAL,"g fracture final");
    EREG(fi_params,SAVE_CONN_DATA,"save connection rates to rst file");
    EREG(fi_params,DISABLE_FRACTURE_CHECK,"disable fracture check");
    EREG(fi_params,COMPRESS_IRR,"compress irregular matrix");
    EREG(fi_params,REMOVE_ISOLATED_REGIONS,"make isolated regions with out wells inactive");
    EREG(fi_params,SAVE_WRATES_TO_ASCII_FILE,"save well rates to ascii file");
    EREG(fi_params,CREATE_HDF5_FILE,"create hdf5 file");
    EREG(fi_params,WRITE_PRESSURE_TO_HDF5,"write pressure to hdf5");
    EREG(fi_params,WRITE_SATURATION_TO_HDF5,"write saturation to hdf5");
    EREG(fi_params,FIX_SOIL_BUG,"fix soil bug");
    EREG(fi_params,DISABLE_CROSSFLOW,"disable crossflow");
    EREG(fi_params, WRITE_PRESSURE_TO_HDF5,"write pressure values to hdf5 file");
    EREG(fi_params, WRITE_SATURATION_TO_HDF5,"write saturation values to hdf5 file");
    EREG(fi_params, WRITE_GAS_OIL_RATIO_TO_HDF5,"write gas-oil ratio values to hdf5 file");
    EREG(fi_params, WRITE_FIP_RESULTS_TO_HDF5,"write data by fip regions to hdf5 file");
    EREG(fi_params, WRITE_WELL_RESULTS_TO_HDF5,"write well data to hdf5 file");
    EREG(fi_params, WRITE_CONN_RESULTS_TO_HDF5,"write connection data to hdf5 file");
    EREG(fi_params, USE_CFL, "compute CFL in fi_operator and use csr_ilu_cfl_prec, used only if BS_BOS_CORE_USE_CSR_ILU_CFL_PREC defined");
    EREG(fi_params, WRITE_CFL_TO_HDF5, "write CFL vector to HDF5, used only if BS_BOS_CORE_USE_CSR_ILU_CFL_PREC defined and USE_CFL turned on");
    EREG(fi_params,B_TOTAL,"all parameters should be before this item");
    EREG(fi_params,fi_params::FI_PARAMS_TOTAL,"size of parameters array");

    resize(fi_params::FI_PARAMS_TOTAL); // to resize names table and table of properties

    set_bool (PRINT_PVTO_TABLE, 1);                    // print PVTO table by default
    set_bool (PRINT_PVTW_TABLE, 1);                    // print PVTW table by default
    set_bool (PRINT_PVTG_TABLE, 1);                    // print PVTW table by default
    set_bool (STORE_PANE_FLOW_RATES, 0);               // do not store flow rates by default
    set_bool (NEW_TS_SELECTION_ALGORITHM, 1);          // use old algorithm
    set_bool (SAVE_STEP_MAPS, 1);                      // save step maps on time step in the rst file
    set_bool (SAVE_RESTART_DATA, 1);                   // save restart data
    set_bool (SAVE_INITIAL_DATA, 1);                   // save initial data
    set_bool (NEWTRANS, 1);                            // use NEWTRAN option
    set_bool (CHECKS_FOR_MULTI_CONN, 1);               // check for
    set_bool (SIMPLE_GET_CELL_SOLUTION, 1);            // read all time steps from *tsl file
#ifndef COMPOSITIONAL
    set_bool (USE_TIMESTEP_FILE, 0);                   // read all time steps from *tsl file
#else
    set_bool (USE_TIMESTEP_FILE, 0);                   // read all time steps from *tsl file
#endif
    set_bool (USE_LOW_SKIN_TRANS_MULT, 1);             // use calculation of multipliers for low skin-factor
    set_bool (SAVE_MESH_DATA, 1);                      // save mesh data
    set_bool (SAVE_CALC_DATA, 0);                      // save calc data
    set_bool (LOAD_CALC_DATA, 0);                      // load calc data
    set_bool (SAVE_NORM_DATA, 0);                      // save cells norms
    set_bool (DENS_JAC_SCALE, 1);                      // scale Jacobian equations
    set_bool (SAVE_CROSSFLOW, 0);
    set_bool (NEWTON_TUNING, 0);
    set_bool (G_FRACTURE, 0);
    set_bool (CLAMP_PRESSURE, 0);
    set_bool (G_FRACTURE_FINAL, 0);
    set_bool (SAVE_CONN_DATA, 1);                      // save connection data
    set_bool (DISABLE_FRACTURE_CHECK, 0);
    set_bool (COMPRESS_IRR, 0);
    set_bool (REMOVE_ISOLATED_REGIONS, 0);             // remove regions by default
    set_bool (SAVE_WRATES_TO_ASCII_FILE, 1);           // save well rates to ascii file by default
    set_bool (CREATE_HDF5_FILE, 1);                    // create hdf5 file
    set_bool (WRITE_PRESSURE_TO_HDF5, 1);              // write pressure to hdf5
    set_bool (WRITE_SATURATION_TO_HDF5, 0);           // write saturation to hdf5
    set_bool (FIX_SOIL_BUG, 1);//0);
    set_bool (DISABLE_CROSSFLOW, 0);
    set_bool (WRITE_PRESSURE_TO_HDF5, 1);
    set_bool (WRITE_SATURATION_TO_HDF5, 1);
    set_bool (WRITE_GAS_OIL_RATIO_TO_HDF5, 1);
    set_bool (WRITE_PLANE_FLOW_RATES_TO_HDF5, 0);
    set_bool (WRITE_FIP_RESULTS_TO_HDF5, 1);
    set_bool (WRITE_WELL_RESULTS_TO_HDF5, 1);
    set_bool (WRITE_CONN_RESULTS_TO_HDF5, 0);
    set_bool (USE_CFL, 0);
    set_bool (WRITE_CFL_TO_HDF5, 0);

    // default value for integer params
    set_int (PVT_INTERP_POINTS, 20);                   // number of interpolation points for PVT tables
    set_int (LIN_ITERS_NUM, 30);                       // max number of linear iteration
    set_int (NEWTON_ITERS_NUM, 12);                    // max number of newtonian iteration
    set_int (NEWTON_ITERS_INC, 6);                     // if niters < maxniters - niters_inc increase ts
    set_int (NEWTON_ITERS_AMG, 20);                    // if niters < niters_amg build full amg setup
    set_int (APPROX_STEPS, 0);                         // number of approximation steps
    set_int (SELECT_SOL_STEPS, 2);                     // number of steps produced to select force term in newton correction (0 -- no selection)
    set_int (LIN_SOLVER_TYPE, FI_LIN_SOLVER_GMRES);    // linear solver type
    set_int (GMRES_ORTONORM_VLEN, 10);                 // number of ortonorm vectors in GMRES
    set_int (AMG_LIN_ITERS_NUM, 10);                   // maximum number of AMG iterations
    set_int (WELL_NEWTON_ITERS_NUM, 1/*20*/);               // maximum number of well internal newton iterations
    set_int (PREC_TYPE, FI_LIN_PREC_CPR);              // preconditioner type for linear solver
    set_int (PREC_TYPE_ONE_PHASE, FI_LIN_PREC_ILU);              // preconditioner type for linear solver (1-phase problem)
    set_int (MIN_CELLS_IN_REGION, 0);
    set_int (SAVE_STEP_DATA_PERIOD, 1);                // by default save step data at all steps

    // default value for floating points params
    set_float (PVT_PRESSURE_RANGE_MAX, (ph_const.atmospheric_pressure * 1000));  // default max pressure in reservoir (atm.)
    set_float (PVT_PRESSURE_RANGE_MIN, (ph_const.atmospheric_pressure * 1));     // default min pressure in reservoir (atm.)
    set_float (PVTO_RS_RANGE_MIN, 1);                  // default min RS
    set_float (PVTO_RS_RANGE_MAX, 400);                // default max RS
    set_float (NEWTON_RESIDUAL, (1.0e-2));             // residual value
    set_float (LIN_SOLV_RESIDUAL, (1.0e-4));           // residual value
    set_float (LIN_SOLV_MATBAL_RESIDUAL, (-1));        // matbal residual value
    set_float (MAX_ALLOWED_LIN_SOLV_RESID, 0.05);      // maximum allowed residual value
    set_float (AMG_RESID, (1.0e-1));                   // AMG preconditioner solver tolerance
    set_float (MAX_TS, 100);                           // maximum allowed time step
    set_float (MIN_TS, (5e-6));                        // minimum allowed time step
    set_float (FIRST_TS, 1);
    set_float (TS_INC_MULT, 2);                        // increment multiplier
    set_float (TS_DEC_MULT, 0.5);                      // decrement multiplier
    set_float (OVERDRAFT, 1.4);                        // overdraft multiplier
    set_float (MAX_P_CORRECTION, 300);//00);                 // pressure correction on newton iteration should be <
    set_float (MAX_S_CORRECTION, 1);                 // saturation correction on newton iteration should be <
    set_float (MAX_RS_CORRECTION, 1000);                 // Rs correction on newton iteration should be <
    set_float (MAX_P_WELL_CORRECTION, 3000);            // well pressure correction on newton iteration should be <
    set_float (PRESSURE_COL_MULT, 1.0);                // mult for pressue column in Jacobian
    set_float (WAT_ROW_MULT, 1);                       // mult for water row in Jacobian (do not change)
    set_float (GAS_ROW_MULT, 1);                       // mult for gas row in Jacobian
    set_float (OIL_ROW_MULT, 1);                       // mult for oil row in Jacobian
    set_float (TS_OMEGA, 0.5);                         // time step controling omega param
    set_float (TS_DP, 60);                             // time step controling pressure change
    set_float (TS_DS, 0.05);                           // time step controling saturation change
    set_float (TS_DRS, 10);                            // time step controling Rs change
    set_float (DP_MIN_CHOP, (1e-5));                   // minimum chop in pressure change
    set_float (DS_MIN_CHOP, (1e-7));                   // minimum chop in saturation change
    set_float (MAX_NORM_ON_TS, 10);
    set_float (DRSDT, (-1));                           // Maximum rate of increase of solution GOR
    set_float (GAS_NORM_MULT, 1);                      // default gas norm multiplier

    // default parameters for compositional model
    // newton converge parameters
    //set_param(COMP_PARAMS_D_MASS_BALANCE, 0.1);                   // Scaled mass balance equation residuals. Equations are scaled by the accumulation terms (component mass dividing timestep size).
    //set_param(COMP_PARAMS_D_PHASE_EQUIL, 0.02);                   // Scaled phase equilibrium relation residuals. These are scaled by the component fugacities in the gas phase
    //set_param(COMP_PARAMS_D_WELL_RATE, 0.001);                    // Scaled well equation residuals, scaled by well input rate
    //set_param(COMP_PARAMS_D_MAX_DP, 0.05);                        // The maximum scaled pressure change (absolute pressure change divided by the average reservoir pressure)
    //set_param(COMP_PARAMS_D_MAX_DS, 0.005);                       // The maximum absolute saturation change
    //set_param(COMP_PARAMS_D_MAX_DXCP, 0.001);                     // The maximum absolute component mole fraction change
    //set_param(COMP_PARAMS_D_TS_DP, 20);                           // time step controling pressure change
    //set_param(COMP_PARAMS_D_TS_DS, 0.2);                          // time step controling saturation change
    //set_param(COMP_PARAMS_D_TS_DXCP, 0.02);                       // time step controling mole fraction change
    //set_param(COMP_PARAMS_D_MAX_P_CORRECTION, 200);               // pressure correction on newton iteration should be <
    //set_param(COMP_PARAMS_D_MAX_S_CORRECTION, 0.5);               // saturation correction on newton iteration should be <
    //set_param(COMP_PARAMS_D_MAX_XCP_CORRECTION, 0.2);             // mole fraction correction on newton iteration should be <

    set_float (MASS_BALANS_ERROR, (1.0e-5));           // maximum allowed mass balans error
    set_float (MAX_NORM_ON_FIRST_N, 1);
    set_float (P_DIMENSION_LESS_SCALE_FACTOR, 0.001);  // scale factor for making pressure dimension less

    // FIXME: set properly value
    set_float (S_RHS_NORM, 0.01);
  }

  const std::string & fi_params::get_params_name (idx_type idx)
  {
    return params_names[idx].name;
  }

  // And special blue-sky objbase class implementations

  fi_params::fi_params (bs_type_ctor_param /*param*/)
  : bs_refcounter(), named_pbase()
  {
    set_default_values ();
  }

  fi_params::fi_params (const fi_params& prop)
  : bs_refcounter(prop), named_pbase(prop)
  {
    if (&prop != this)
      *this = prop;

    set_default_values ();
  }

  fi_params::~fi_params ()
  {
    set_default_values ();
  }

  const fi_params &fi_params::operator=(const fi_params &src)
  {
    if (&src == this)
      return *this;

    named_pbase::operator=(src);
    return *this;
  }

  const fi_params &fi_params::operator+=(const fi_params &src)
  {
    if (&src == this)
      return *this;

    named_pbase::operator+=(src);
    return *this;
  }

  /*!
  * \brief read params from keyword file in format
  *        -------------------------------------------------------
  *        D_PARAM_NAME  D_VALUE
  *        ...
  *        I_PARAM_NAME  I_VALUE
  *        ...
  *        B_PARAM_NAME  B_VALUE
  *        ...
  *        /
  *        where D_PARAM_NAME specified in d_param_name_table, I_PARAM_NAME specified in i_param_name_table,
  *        B_PARAM_NAME specified in b_param_name_table, D_VALUE -- floating point value,
  *        I_VALUE -- integer value, B_VALUE -- should be word "ON" or word "OFF"
  *        list of params should be finished by "/"
  * \param r -- pointer to file master class
  */
  void fi_params::read_from_keyword_data_file (sp_reader_t r)
  {
    char key[CHAR_BUF_LEN];
    char key1[CHAR_BUF_LEN];
    char buf[CHAR_BUF_LEN];
    int len, i;
    char *strt, *end_ptr;
    t_double d_tmp;
    t_long i_tmp;
    int flag;
    std::ostringstream out_s;
    const sp_reader_t &lreader (r);


    //BSOUT << priority (sn::read, lev::med) << "Process params:" << bs_end;
    for (;;)
      {
        len = lreader->read_line (buf, CHAR_BUF_LEN);
        if (len == YS_EOF)
          return;
        if (len < 0)              // break if EOF
          {
            switch (len)
              {
              case YS_NOTHING_TO_READ:
                break;
              default:
                return;
              }
            break;
          }
        /* End of list  */
        if (buf[0] == '/')
          break;

        // Read connection info
        lreader->unwrap (buf, 2);
        //if ((current_error_num = lreader->) != YS_SUCCESS)
        //  {
        //    out_s << "Error in " << lreader->get_prefix() << ": bad string " << buf;
        //    BS_ASSERT(false) (out_s.str());
        //    throw bs_exception("fi_params class",out_s.str().c_str());
        //    BOSERR (section::read_data, level::error) << out_s << bs_end;
        //  }
        strt = buf;
        lreader->scanf_s (strt, &end_ptr, key);
        //if (lreader-> != YS_SUCCESS)
        //  {
        //    out_s << "Error in " << lreader->get_prefix() << ": can't read parameter name from " << strt;
        //    BS_ASSERT(false) (out_s.str());
        //    throw bs_exception("fi_params class",out_s.str().c_str());
        //    BOSERR (section::read_data, level::error) << out_s << bs_end;
        //  }
        strt = end_ptr;
        flag = 1;

        for (i = FI_PARAMS_START+1; i < D_TOTAL && flag; ++i)
          {
            if (!strcmp (key, this->get_name((named_pbase::idx_type)i).c_str()))
              {
                // read value
                lreader->scanf_fp (strt, &end_ptr, &d_tmp);
                //if (lreader-> != YS_SUCCESS)
                //  {
                //    out_s << "Error in " <<lreader->get_prefix() << ": can't read value for parameter " << key
                //    << " from " << strt;
                //    BS_ASSERT(false) (out_s.str());
                //    throw bs_exception("fi_params class",out_s.str().c_str());
                //    BOSERR (section::read_data, level::error) << out_s << bs_end;
                //  }
                BOSOUT (section::read_data, level::medium) << key << " = " << d_tmp << bs_end;

                this->set_float ((idx_type)i,  d_tmp);
                flag = 0;
              }

          }
        // loop through integer params
        for (i = D_TOTAL+1; i < I_TOTAL && flag; ++i)
          {
            if (!strcmp (key, this->get_name((named_pbase::idx_type)i).c_str()))
              {
                // read value
                lreader->scanf_int (strt, &end_ptr, &i_tmp);
                //if (lreader-> != YS_SUCCESS)
                //  {
                //    out_s << "Error in " <<lreader->get_prefix() << ": can't read value for parameter " << key
                //    << " from " << strt;
                //    BS_ASSERT(false) (out_s.str());
                //    throw bs_exception("fi_params class",out_s.str().c_str());
                //    BOSERR (section::read_data, level::error) << out_s.str() << bs_end;
                //  }
                BOSOUT (section::read_data, level::medium) << key << " = " << i_tmp << bs_end;

                this->set_int ((idx_type)i,  i_tmp);
                flag = 0;
              }
          }
        // loop through boolian params
        for (i = I_TOTAL; i < B_TOTAL && flag; ++i)
          {
            if (!strcmp (key, this->get_name((named_pbase::idx_type)i).c_str()))
              {
                // read value
                lreader->scanf_s (strt, &end_ptr, key1);
                //if (lreader-> != YS_SUCCESS)
                //  {
                //    out_s << "Error in " << lreader->get_prefix() << ": can't read value for parameter " << key
                //    << " from " << strt << " value should be \'ON\' or \'OFF\' or \'YES\' or \'NO\'";
                //    BS_ASSERT(false) (out_s.str());
                //    throw bs_exception("fi_params class",out_s.str().c_str());
                //    BOSERR (section::read_data, level::error) << out_s.str() << bs_end;
                //  }
                if (!strcmp (key1, "ON") || !strcmp (key1, "YES") || !strcmp (key1, "1"))
                  {
                    BOSOUT (section::read_data, level::medium) << key << " = " << key1 << bs_end;
                    this->set_param ((idx_type)i, 1);
                  }
                else if (!strcmp (key1, "OFF") || !strcmp (key1, "NO") || !strcmp (key1, "0"))
                  {
                    BOSOUT (section::read_data, level::medium) << key << " = " << key1 << bs_end;
                    this->set_param ((idx_type)i, 0);
                  }
                else
                  {
                    out_s << "Error in " << lreader->get_prefix() << ": can't read value for parameter " << key
                    << " from " << strt << " value should be \'ON\' or \'OFF\' or \'YES\' or \'NO\'";
                    BS_ASSERT(false) (out_s.str());
                    throw bs_exception("fi_params class",out_s.str().c_str());
                    BOSERR (section::read_data, level::error) << out_s << bs_end;
                  }
                flag = 0;
              }
          }

      }
  }


  BLUE_SKY_TYPE_STD_CREATE(fi_params)
  BLUE_SKY_TYPE_STD_COPY(fi_params)
  BLUE_SKY_TYPE_IMPL_SHORT(fi_params, objbase, "fi_params class")
}
