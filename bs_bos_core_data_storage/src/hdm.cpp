#include "bs_bos_core_data_storage_stdafx.h"

#include "hdm.h"
//#include "main_def.h"
#include "arrays.h"
//#include "localization.h"
#include "arrays_tables.h"
#include "data_class.h"

#include "write_time_to_log.h"
#include <stdio.h>


    #include <iostream>

    using namespace boost;
    using namespace std;
 
/////////////////////////////////////////////////////////////////////////////
// You can add new variable just find next Words
// for symple variable                  -- VAR_V
// for array                                    -- ARRAY_V
// for strings                                  -- STRING_V

//! Default minimal DX value allowed for active cells
#define DEFAULT_MINIMAL_DX              1.e-6
//! Default minimal DY value allowed for active cells
#define DEFAULT_MINIMAL_DY              1.e-6
//! Default minimal DZ value allowed for active cells
#define DEFAULT_MINIMAL_DZ              1.e-10
//! Default minimal PERM value allowed for active cells
#define DEFAULT_MINIMAL_PERM            1e-16
#define DEFAULT_MINIMAL_NTG             1.e-8
#define DEFAULT_MINIMAL_TOPS            1.e-8
//! default blank value
#define DEFAULT_BLANK_VALUE             1.e43
//! check blank value
#define CHECK_BLANK_VALUE(X)  (fabs ((X) - DEFAULT_BLANK_VALUE) < 10)

//! Default minimal PORO value allowed for active cells
#define DEFAULT_MINIMAL_PORO                        1.e-10
#define DEFAULT_MINIMAL_PORO_VOLUME_MULT            1.e-12
/////////////////////////////////////////////////////////////////////////////


// Please don't delete commented old style outputs .. like rep->blah blah blah
namespace blue_sky
  {
  
#define DM_ASSERT_EXCEPTION \
    BS_ASSERT(false) (out_s.str());\
    throw bs_exception("Hydrodynamic model class",out_s.str().c_str());

  
  hdm::~hdm ()
  {

  }

  
  hdm::hdm(bs_type_ctor_param /*param*/): lkeeper ("C", LC_ALL)
  {
    this->reader = BS_KERNEL.create_object("bos_reader");
    this->data = BS_KERNEL.create_object(idata::bs_type());
    init_proc_params();
    this->km = BS_KERNEL.create_object(keyword_manager::bs_type());
    
    this->pvt_3p_ = BS_KERNEL.create_object ("pvt_3p");
    this->scal_3p_ = BS_KERNEL.create_object ("scal_3p");
    this->event_manager_ = BS_KERNEL.create_object ("event_manager");
    this->well_pool_ = BS_KERNEL.create_object ("sql_well");
  }

  hdm::hdm(const hdm& src):bs_refcounter (src), lkeeper ("C", LC_ALL)
  {
    *this = src;
  }
 
  void
  hdm::init_fluids(t_int n_pvt_regions, t_int n_scal_regions)
  {
    int n_phases; 
    
    pvt_3p_->init_pvt_arrays (n_pvt_regions, data->props->get_b("oil_phase"),
											 data->props->get_b("gas_phase"),
                                             data->props->get_b("water_phase")
                                             );
    
    n_phases = data->props->get_b("oil_phase");
    n_phases += data->props->get_b("water_phase");
    n_phases += data->props->get_b("gas_phase");
    if (n_phases > 1)
      {
        scal_3p_->init_scal_input_table_arrays (n_scal_regions, data->props->get_b("oil_phase"),
                                                                data->props->get_b("gas_phase"),
                                                                data->props->get_b("water_phase"));
      }
  }
  
  void 
  hdm::init_proc_params()
  {
    data->proc_params->add_property_b (1, "PRINT_PVTO_TABLE",               "if this item is set print PVTO to LOG file");
    data->proc_params->add_property_b (1, "PRINT_PVTW_TABLE",               "print pvtw interpolation table");
    data->proc_params->add_property_b (1, "PRINT_PVDG_TABLE",               "print pvtg interpolation table");
    data->proc_params->add_property_b (0, "STORE_PANE_FLOW_RATES",          "");
    data->proc_params->add_property_b (1, "NEW_TS_SELECTION",               "use new algorithm for ts selection");
    data->proc_params->add_property_b (1, "SAVE_STEP_MAPS",                 "save maps on time step in the rst file");
    data->proc_params->add_property_b (1, "SAVE_INITIAL_DATA",              "save initial data in the rst file");
    data->proc_params->add_property_b (1, "SAVE_RESTART_DATA",              "save restart information");
    data->proc_params->add_property_b (1, "NEWTRANS",                       "correct transmissibility calculation between cells");
    data->proc_params->add_property_b (1, "CHECK_FOR_MULTY_CONN",           "check for multi connections in one cell");
    data->proc_params->add_property_b (1, "SIMPLE_GET_CELL_SOLUTION",       "save time steps into *tsl file");
    data->proc_params->add_property_b (0, "USE_TIMESTEP_FILES",             "enable simple solution restore algorithm");
    data->proc_params->add_property_b (1, "USE_LOW_SKIN_TRANS_MULT",        "if this item is set use calculation of trans mults for low skin");
    data->proc_params->add_property_b (1, "SAVE_MESH_DATA",                 "save initial data in the rst file");
    data->proc_params->add_property_b (0, "SAVE_CALC_DATA",                 "save solution vector at each large time step into .dmp file");
    data->proc_params->add_property_b (0, "LOAD_CALC_DATA",                 "load initial guesses for Newton iterations from .dmp file");
    data->proc_params->add_property_b (0, "SAVE_NORM_DATA",                 "save cells norms to .cel file");
    data->proc_params->add_property_b (1, "DENS_JAC_SCALE",                 "enable densities based multiplier's for gas and oil equation's in Jacobian");
    data->proc_params->add_property_b (0, "NEWTON_TUNING",                  "");
    data->proc_params->add_property_b (0, "SAVE_CROSSFLOW",                 "");
    data->proc_params->add_property_b (0, "G_FRACTURE",                     "");
    data->proc_params->add_property_b (0, "CLAMP_PRESSURE",                 "clamp pressure values at each Newton iteration between possible min and max");
    data->proc_params->add_property_b (0, "G_FRACTURE_FINAL",               "");
    data->proc_params->add_property_b (1, "SAVE_CONN_DATA",                 "save connection rates to rst file");
    data->proc_params->add_property_b (0, "DISABLE_FRACTURE_CHECK",         "");
    data->proc_params->add_property_b (0, "COMPRESS_IRR",                   "compress irregular matrix");
    data->proc_params->add_property_b (0, "REMOVE_ISOLATED_REGIONS",        "make isolated regions with out wells inactive");
    data->proc_params->add_property_b (1, "SAVE_WRATES_TO_ASCII_FILE",      "save well rates to ascii file");
    data->proc_params->add_property_b (1, "CREATE_HDF5_FILE",               "create hdf5 file");
    data->proc_params->add_property_b (1, "WRITE_PRESSURE_TO_HDF5",         "write pressure to hdf5");
    data->proc_params->add_property_b (1, "WRITE_SATURATION_TO_HDF5",       "write saturation to hdf5");
    data->proc_params->add_property_b (1, "FIX_SOIL_BUG",                   "");
    data->proc_params->add_property_b (0, "DISABLE_CROSSFLOW",              "");
    data->proc_params->add_property_b (0, "DEBUG_EQUIL",                    "write equil debug info");
    data->proc_params->add_property_b (1, "WRITE_GAS_OIL_RATIO_TO_HDF5",    "write gas_oil_ratio to hdf5");
    data->proc_params->add_property_b (0, "FIX_GRID",                       "calc dx dy dz as in PPP");
    data->proc_params->add_property_b (1, "WRITE_WELL_RESULTS_TO_HDF5",     "save well results to hdf5 file");
    data->proc_params->add_property_b (1, "WRITE_CONN_RESULTS_TO_HDF5",     "save connection results to hdf5 file");
    data->proc_params->add_property_b (1, "WRITE_FIP_RESULTS_TO_HDF5",      "save fip results to hdf5 file");
    data->proc_params->add_property_b (0, "WRITE_ECL_CONN_DATA",            "write connection data in ECL COMPDAT format");
    data->proc_params->add_property_b (0, "WRITE_INITIAL_DATA_TO_HDF5",     "save initial data to hdf5 file");
    data->proc_params->add_property_b (0, "READ_INITIAL_DATA_FROM_HDF5",    "read initial data from hdf5 file");
    data->proc_params->add_property_b (1, "WRITE_MESH_DATA_TO_HDF5",        "save mesh data to hdf5 file");
    data->proc_params->add_property_b (1, "WRITE_PLANE_FLOW_RATES_TO_HDF5", "write plane flow rates to hdf5 file");
    data->proc_params->add_property_b (0, "DISABLE_GRAVITY",                "if ON, disable gravitation");
    data->proc_params->add_property_b (0, "PARALLEL_ILU",                   "if ON, enable parallel ilu preconditioner");
    data->proc_params->add_property_b (0, "USE_IRR_MATRIX_IN_ILU",          "if OFF, only diagonal from irregular matrix will be used in ILU decomposition ");
    data->proc_params->add_property_b (0, "FIX_WELL_DENSITY_BUG",           "");
    data->proc_params->add_property_b (1, "WBP",                            "calc average well pressure only in well blocks");
    data->proc_params->add_property_b (0, "WBP4",                           "calc average well pressure only in 4 well neighbours");
    data->proc_params->add_property_b (0, "WBP5",                           "calc average well pressure in well blocks and 4 well neighbours");
    data->proc_params->add_property_b (1, "WBP9",                           "calc average well pressure in well blocks and 8 well neighbours");
    data->proc_params->add_property_b (0, "PURE_NEWTON",                    "disable all hacks ");
    data->proc_params->add_property_b (0, "WRITE_CNORM_TO_ASCII_FILE",      "in true C norm of all components will be save to ascii file in format i j k Cw Cg Co");
    data->proc_params->add_property_b (0, "P_INIT_APPROX",                  "if true also initialize pressure");
    data->proc_params->add_property_b (0, "PV_WEIGHTED_PRESSURE",           "pore volume weighted pressure calculation");
    data->proc_params->add_property_b (0, "CRS_SCAL_SCALE",                 "enable 3-point scal scale");
    data->proc_params->add_property_b (1, "SAVE_PCW",                       "save calculated PCW to hdf5");
    data->proc_params->add_property_b (0, "SAVE_PCG",                       "save calculated PCG to hdf5");
    data->proc_params->add_property_b (0, "NEGATIVE_BHP",                   "WELL BHP could be negative if ON");
    data->proc_params->add_property_b (0, "WRITE_PLANE_FLOW_VOLS_TO_HDF5",  "write total plane mass flow volumes to hdf5 file");
    data->proc_params->add_property_b (0, "SET_ACTIVE_FRACTURE_CELLS",      "set actnum = 1 for cells, where fractures defined");
    data->proc_params->add_property_b (0, "FRACTURE_HORIZ_WELL",            "enable fracture for horiz wells, which grows up and down by z-layer ");
    data->proc_params->add_property_b (1, "WRITE_BHP_0_FOR_CLOSED_WELLS",   "if true, write well_bhp=0 to h5 file for shutted wells");
    data->proc_params->add_property_b (1, "MULTI_CONNECTION",               "if true, allow adding multiple connections to one cell");
    data->proc_params->add_property_b (0, "FIX_PCW_BUG",                    "if true, set default psw and pcg values to 0.0");
    data->proc_params->add_property_b (0, "FIX_WPIMULT_BUG",                "if true, do not apply WPIMULT to fracture's connections");
    data->proc_params->add_property_b (1, "FIX_ARITH_BUG",                  "if true, add brackets around function and around it's arguments in arithmetic");
    data->proc_params->add_property_b (0, "USE_NITERS_VOLUMETRIC_NORM",     "use volumetric norm calculation as in Eclipse");
    data->proc_params->add_property_b (0, "FRACTURE_USE_FABS",               "use absolute coordinates of points in calculating fracture's connection factors");



    data->proc_params->add_property_i (20,    "PVT_INTERP_POINTS",        "number of interpolation points");
    data->proc_params->add_property_i (12,    "NITERS_NUM",               "maximum allowed newton iterations number");
    data->proc_params->add_property_i (30,    "LITERS_NUM",               "maximum allowed linear iterations number");
    data->proc_params->add_property_i (6,     "NITERS_INC",               "number of newton iterations to increment step length");
    data->proc_params->add_property_i (20,    "NITERS_AMG",               "number of newton iterations to build full amg setup");
    data->proc_params->add_property_i (0,     "APPROX_STEPS",             "number of approximation steps");
    data->proc_params->add_property_i (2,     "SELECT_SOL_STEPS",         "number of steps to select newton correction force term");
    data->proc_params->add_property_i (1,     "LSOLV_TYPE",               "type of linear solver");
    data->proc_params->add_property_i (10,    "GMRES_STACK_LEN",          "number of vectors in GMRES to ortonorm");
    data->proc_params->add_property_i (10,    "AMG_LITERS_NUM",           "maximum allowed AMG solver iterations");
    data->proc_params->add_property_i (20,    "WELL_NITERS_NUM",          "maximum number of well internal newton iterations");
    data->proc_params->add_property_i (0,     "PREC_TYPE",                "type of the preconditioner for linear solver");
    data->proc_params->add_property_i (0,     "MIN_CELLS_IN_REGION",      "minimum allowed cells in region for ACTNUM = 1");
    data->proc_params->add_property_i (1,     "SAVE_STEP_DATA_PERIOD",    "save step data every n step");
    data->proc_params->add_property_i (10,    "SAVE_WELL_RESULTS_PERIOD", "save well results every n step");
    data->proc_params->add_property_i (10,    "SAVE_FIP_RESULTS_PERIOD",  "save fip results every n step");
    data->proc_params->add_property_i (20000, "FRACTURE_SERIES_NUMBER",   "number of elemements in calculating of series for fracture calculation");
    data->proc_params->add_property_i (3,     "NEWTON_ITERS_GCONTROL",    "number of newton iterations to consider group control (equal to NUPCOL, used with GCONINJE and GPMAINT)");
    data->proc_params->add_property_i (0,     "FRACTURE_ANGLE_AXIS",      "fracture angle counting: 0-from mesh-based X-axis (by first cell), 1-from X axis, 2-from Y axis (azimut)");
    data->proc_params->add_property_i (2,     "1PHASE_LSOLV_TYPE",        "type of linear solver for 1 phase systems");
    data->proc_params->add_property_i (4,     "1PHASE_PREC_TYPE",         "type of the preconditioner for linear solver for 1 phase systems");
    data->proc_params->add_property_i (0,     "WRITE_NORM_TO_HDF5",       "write norm to hdf5 file flag. 1-write each large step, 2-write at each newton iteration");
    data->proc_params->add_property_i (0,     "TIMESTEP_ON_NORM_SELECT",  "timestep depends from cnorm calculated on first iteration");
                                               
                                                      
    data->proc_params->add_property_f (1000,    "PVT_PRESSURE_MAX",               "PVT tables maximum pressure");
    data->proc_params->add_property_f (1,       "PVT_PRESSURE_MIN",               "PVT tables minimal pressure");
    data->proc_params->add_property_f (400,       "PVTO_RS_MAX",                    "PVTO maximum RS for slop table");
    data->proc_params->add_property_f (1,     "PVTO_RS_MIN",                    "PVTO minimum RS for slop table");
    data->proc_params->add_property_f (100,     "TS_MAX",                         "maximum allowed time step length");
    data->proc_params->add_property_f (5.0e-6,  "TS_MIN",                         "minimum allowed time step length");
    data->proc_params->add_property_f (1,       "TS_FIRST",                       "first simulation time step length");
    data->proc_params->add_property_f (1.0e-4,  "LITERS_RESID",                   "tolerance for linear solver");
    data->proc_params->add_property_f (1.0e-2,  "NITERS_RESID",                   "tolerance for newton process");
    data->proc_params->add_property_f (2,       "TS_INC_MULT",                    "multiplier for incrementing time step length");
    data->proc_params->add_property_f (0.5,     "TS_DEC_MULT",                    "multiplier for decrementing time step length");
    data->proc_params->add_property_f (1.4,     "TS_OVERDRAFT",                   "overdraft factor (any time step could be multiplied by this factor to achieve end of report time step)");
    data->proc_params->add_property_f (300,     "P_CORR_MAX",                     "maximum allowed newton correction for pressure");
    data->proc_params->add_property_f (1,       "S_CORR_MAX",                     "maximum allowed newton correction for saturation");
    data->proc_params->add_property_f (1000,    "RS_CORR_MAX",                    "maximum allowed newton correction for gas oil ratio");
    data->proc_params->add_property_f (3000,    "WELL_P_CORR_MAX",                "maximum allowed newton correction for well pressure");
    data->proc_params->add_property_f (1.0,     "WAT_ROW_MULT",                   "multiplier for water equation in Jacobian");
    data->proc_params->add_property_f (1.0,     "GAS_ROW_MULT",                   "multiplier for gas equation in Jacobian");
    data->proc_params->add_property_f (1.0,     "OIL_ROW_MULT",                   "multiplier for oil equation in Jacobian");
    data->proc_params->add_property_f (1.0,     "P_COL_MULT",                     "multiplier for pressure derivates column in Jacobian");
    data->proc_params->add_property_f (0.5,     "TS_OMEGA",                       "omega coef in time step controling");
    data->proc_params->add_property_f (60,      "TS_DP",                          "pressure change for time step controling");
    data->proc_params->add_property_f (0.05,     "TS_DS",                          "saturation change for time step controling");
    data->proc_params->add_property_f (1e-5,    "DP_MIN_CHOP",                    "minimum pressure chop");
    data->proc_params->add_property_f (1e-7,    "DS_MIN_CHOP",                    "minimum saturation chop");
    data->proc_params->add_property_f (0.05,    "LITERS_MAX_RESID",               "maximum allowed residual");
    data->proc_params->add_property_f (1.0e-1,  "AMG_RESID",                      "tolerance for AMG preconditioner");
    data->proc_params->add_property_f (10,      "TS_DRS",                         "Rs change for time step controling");
    data->proc_params->add_property_f (10,      "MAX_NORM_ON_TS",                 "if norm on time step is greater than this value restart occur");
    data->proc_params->add_property_f (-1,      "DRSDT",                          "Maximum rate of increase of solution GOR");
    data->proc_params->add_property_f (1,       "GAS_NORM_MULT",                  "multiplier for gas norm");
    data->proc_params->add_property_f (0.99995, "FRACTURE_LAMBDA",                "the minimum value for lambda in fracture calculation, than less lambda than less accuracy, must be more than 0.99");
    data->proc_params->add_property_f (0.01,    "COMP_MASS_BALANCE",              "Scaled mass balance equation residuals. Equations are scaled by the accumulation terms (component mass dividing timestep size).");
    data->proc_params->add_property_f (0.02,    "COMP_PHASE_EQUIL",               "Scaled phase equilibrium relation residuals. These are scaled by the component fugacities in the gas phase");
    data->proc_params->add_property_f (0.0001,  "COMP_MAX_DP",                    "The maximum scaled pressure change (absolute pressure change divided by the average reservoir pressure)");
    data->proc_params->add_property_f (0.005,   "COMP_MAX_DS",                    "The maximum absolute saturation change");
    data->proc_params->add_property_f (0.001,   "COMP_MAX_DXCP",                  "The maximum absolute component mole fraction change");
    data->proc_params->add_property_f (20,      "COMP_TS_DP",                     "time step controling pressure change");
    data->proc_params->add_property_f (0.2,     "COMP_TS_DS",                     "time step controling saturation change");
    data->proc_params->add_property_f (0.02,    "COMP_TS_DXCP",                   "time step controling mole fraction change");
    data->proc_params->add_property_f (200,     "COMP_MAX_P_CORRECTION",          "maximum allowed newton correction for pressure");
    data->proc_params->add_property_f (0.5,     "COMP_MAX_S_CORRECTION",          "maximum allowed newton correction for saturation");
    data->proc_params->add_property_f (0.2,     "COMP_MAX_XCP_CORRECTION",        "maximum allowed newton correction for mole fraction");
    data->proc_params->add_property_f (1.0e-5,  "MASS_BALANS_ERROR",              "maximum allowed mass balans error");
    data->proc_params->add_property_f (1,       "MAX_NORM_ON_FIRST_N",            "maximum allowed norm on first newton iteration");
    data->proc_params->add_property_f (0.001,   "P_DIMENSION_LESS_SCALE_FACTOR",  "scale factor for dimension less pressure");
    data->proc_params->add_property_f (0.3,     "APPL_CHOP_PERC",                 "tipicaly (0.25-0.4)");
    data->proc_params->add_property_f (0.2,     "APPL_CHOP_BND",                  "tipicaly (0.1-0.25)");
    data->proc_params->add_property_f (0.0,     "FRACTURE_MIN_DELTA",             "if dist(A-B) < min_delta => connection will not added to this block, A,B- fracture and block intersection points");
    data->proc_params->add_property_f (100,     "TIMESTEP_ON_NORM_RESIDUAL_MULT", "Multiplier for NEWTON_RESIDUAL");
    data->proc_params->add_property_f (0.1,     "TIMESTEP_ON_NORM_MIN_MULT",      "minimum multiplier for current timestep");
    data->proc_params->add_property_f (0.5,     "TIMESTEP_ON_NORM_MAX_MULT",      "maximum multiplier for current timestep");
    data->proc_params->add_property_f (0.01,    "TIMESTEP_ON_NORM_RESTART_MULT",  "restart multiplier for current timestep");
    data->proc_params->add_property_f (100,     "MAX_PCW_WARNING_VALUE",          "maximum PCW value to print warning messages in log-file");
  }

  void 
  hdm::init_equil (t_int n_equil_regions)
  {
    int n_phases;
    n_phases = data->props->get_b("oil_phase");
    n_phases += data->props->get_b("water_phase");
    n_phases += data->props->get_b("gas_phase");
    
    this->equil_model_ = BS_KERNEL.create_object ("equil_model_depth"); 
    equil_model_->init_equil_model (n_equil_regions, n_phases);
  }
  
   void
  hdm::init(const std::string &model_name)
  {
    smart_ptr <hdm_iface, true> hdm = this;
    keyword_params kp;
    
    kp.hdm = this;
    if (model_name == std::string(":memory:")) {
        this->well_pool_->open_db (model_name); //"mem_well_pool.db");
        //this->well_pool_->open_db ("mem2_well_pool.db");
      }
    else
      {
        data->h5_pool->open_file (model_name + "_rex.h5");
        this->well_pool_->open_db (model_name + "_well_pool.db");
      }
    data->props->add_property_s(model_name.substr(0,model_name.find_last_of("/\\") + 1), "model_path", "model_path");
    km->init(hdm);
    
    
    switch (data->props->get_i("mesh"))
    {
      case 0:
        km->handle_keyword_reactor ("MESH_IJK", kp);
        break;
      case 1:
        km->handle_keyword_reactor ("MESH_GRDECL", kp);
        break;   
          
     default:
        bs_throw_exception ("init: wrong mesh choice");  
    }
    
    switch (data->props->get_i("init"))
    {
      case 0:
        km->handle_keyword_reactor ("EXPLICIT_MODEL", kp);
        break;
      case 1:
        km->handle_keyword_reactor ("EQUIL_MODEL", kp);
        break;   
          
     default:
        bs_throw_exception ("init: wrong mesh choice");  
    }
    
  }
  
  
  void 
  hdm::read_keyword_file(const std::string filename)
  {
    char buf[CHAR_BUF_LEN];
    char key[CHAR_BUF_LEN];
    keyword_params kp;
    int flag;
    int len;
    BS_SP(hdm_iface) sp_hdm;
    sp_hdm = this;
    std::string model_name;
    
    write_time_to_log init_time ("Read model", "");
  
    reader->open (filename.c_str (), filename.c_str ());
    model_name = filename.substr(0,filename.find_last_of("."));
    data->props->add_property_s(filename.substr(0,filename.find_last_of("/\\") + 1), "model_path", "model_path");
    
    data->h5_pool->open_file (model_name + ".h5");
    well_pool_->open_db (":memory:");
    
    km->init(sp_hdm);
    kp.hdm = this;

    // start of loop for data file reading
    flag = 1;

    for (; flag;)
      {
        // reading keyword
        len = reader->read_line (buf, CHAR_BUF_LEN);
        if (len < 0)              // break if EOF
          {
            switch (len)
              {
              case YS_NOTHING_TO_READ:
                break;
              default:
                return;
              }
            //rep->print (LOG_READ_SECTION, LOG_DEBUG, "Finish reading\n");
            BOSOUT (section::read_data, level::low) << "Finish reading" << bs_end;
            break;
          }
        if (sscanf (buf, "%s", key) != 1) // look from buffer keyword
          continue;

        std::string keywrd(key);

        if (key[0] == '/')
          {
            continue;
          }

        if (keywrd == "END")
          {
            BOSOUT (section::read_data, level::low) << "Finish reading with END keyword" << bs_end;
            break;
          }
          
        // do nat handle numbers
        if (keywrd.find_first_of("1234567890") == std::string::npos)
          {
            km->handle_keyword (keywrd, kp);
          }
      }
  }


  //! macros for checking
#define CHECK_FOR_TWO_PARAMS(VAL,ACTNUM,MIN_DEF,COUNTER)                                \
  if ((VAL) < (MIN_DEF) || CHECK_BLANK_VALUE(VAL))                                      \
    {                                                                                   \
      if ((ACTNUM))                                                                     \
        {                                                                               \
          t_long ix, iy, iz;                                                               \
          iz = i / (nx * ny);                                 \
          iy = (i - iz * (nx * ny)) / nx;        \
          ix = i - iz * (nx * ny) - iy * nx;     \
          BOSWARN (section::check_data, level::low)                                     \
            << "blocks[" << (ix + 1) << ", " << (iy + 1) << ", " << (iz + 1)            \
            << "] will be set inactive because of " << #VAL << " = " << VAL << bs_end;  \
          (ACTNUM) = 0;                                                                 \
          ++(COUNTER);                                                                  \
          continue;                                                                     \
        }                                                                               \
    }

#define CHECK_FOR_THREE_PARAMS(VAL,ACTNUM,MIN_DEF,MAX_DEF,COUNTER)                      \
  if ((VAL) < (MIN_DEF) || CHECK_BLANK_VALUE(VAL) || (VAL) > (MAX_DEF))                 \
    {                                                                                   \
      if ((ACTNUM))                                                                     \
        {                                                                               \
          int ix, iy, iz;                                                               \
          iz = i / (nx * ny);                                                           \
          iy = (i - iz * (nx * ny)) / nx;                                               \
          ix = i - iz * (nx * ny) - iy * nx;                                            \
          BOSWARN (section::check_data, level::low)                                     \
            << "blocks[" << (ix + 1) << ", " << (iy + 1) << ", " << (iz + 1)            \
            << "] will be set inactive because of " << #VAL << bs_end;                  \
          (ACTNUM) = 0;                                                                 \
          ++(COUNTER);                                                                  \
          continue;                                                                     \
        }                                                                               \
    }

#define CHECK_DIRECTION_FOR_XYZ(case,num)   \
  if ((case) && actnum [(num)])             \
    {                                       \
      n++;                                  \
      valx += dx[(num)];                    \
      valy += dy[(num)];                    \
      valz += dz[(num)];                    \
    }

#define CHECK_DIRECTION_FOR_TOPS(case,num)  \
  if ((case) && actnum[(num)])              \
    {                                       \
      n++;                                  \
      valtops += tops[(num)];               \
    }

  /*
  
  void hdm::update_geometry() const
    {
      int i, ix, iy, iz, n;
      std::ostringstream out_s;
      int nb = data->nx * data->ny * data->nz;

      int n_tops, n_xyz;
      double dx_var, dy_var, dz_var, tops_var;
      double valx, valy, valz, valtops;

      array_float16_t dx = tools::get_non_empty (((*data->d_map)[DX].array));
      array_float16_t dy = tools::get_non_empty (((*data->d_map)[DY].array));
      array_float16_t dz = tools::get_non_empty (((*data->d_map)[DZ].array));
      array_float16_t tops = tools::get_non_empty (((*data->d_map)[TOPS].array));
      array_uint8_t   actnum = tools::get_non_empty (((*data->i_map)[ACTNUM].array));

      // In case of restart we got all geometry from mesh restored from
      // the restart file. So we skip this test.
      if (data->restart)
        return;

      if (!nb)
        {
          bs_throw_exception ("update_geometry: all dimensions must be greater than 0");
        }

      dx_var = dy_var = dz_var = tops_var = 0.;
      valx = valy = valz = valtops = 0.;
      n_tops = n_xyz = 0;
      if (data->geom_def_flag == GEOM_FLAG_DX_DY_DZ)
        {
          for (iz = 0, i = 0; iz < data->nz; ++iz)
            {
              for (iy = 0; iy < data->ny; ++iy)
                {
                  for (ix = 0; ix < data->nx; ++ix, ++i)
                    {
                      if (actnum[i])
                        {
                          dx_var += dx[i];
                          dy_var += dy[i];
                          dz_var += dz[i];
                          n_xyz++;
                        }
                    }
                }
            }

          if (!n_xyz)
            {
              bs_throw_exception ("update_geometry: No active blocks");
            }

          for (iy = 0, i = 0; iy < data->ny; ++iy)
            {
              for (ix = 0; ix < data->nx; ++ix, ++i)
                {
                  if (actnum[i])
                    {
                      tops_var += tops[i];
                      ++n_tops;
                    }
                }
            }
          if (!n_tops)
            {
              bs_throw_exception ("update_geometry: No active blocks");
            }

          dx_var /= n_xyz;
          dy_var /= n_xyz;
          dz_var /= n_xyz;
          tops_var /= n_tops;
          for (iy = 0, i = 0; iy < data->ny; ++iy)
            {
              for (ix = 0; ix < data->nx; ++ix, ++i)
                {
                  if (!actnum[i])
                    {
                      valtops = 0.;
                      n = 0;
                      CHECK_DIRECTION_FOR_TOPS (ix > 0, i - 1);
                      CHECK_DIRECTION_FOR_TOPS (iy > 0, i - data->nx);
                      CHECK_DIRECTION_FOR_TOPS (ix < data->nx - 1 , i + 1);
                      CHECK_DIRECTION_FOR_TOPS (iy < data->ny - 1 , i + data->nx);
                      if (n)
                        {
                          tops[i] = valtops / n;
                        }
                      else
                        {
                          tops[i] = tops_var;

                          BOSWARN (section::check_data, level::low) 
                            << "new tops in node (" << (ix + 1)
                            << ", " << (iy + 1) << ", " << (iz + 1) << ") = " << tops[i]
                            << " (average = " << tops_var << ")" << bs_end;
                        }
                    }
                }
            }

          for (iz = 0, i = 0; iz < data->nz; ++iz)
            {
              for (iy = 0; iy < data->ny; ++iy)
                {
                  for (ix = 0; ix < data->nx; ++ix, ++i)
                    {
                      if (!actnum[i])
                        {
                          valx = valy = valz = 0.;
                          n = 0;
                          CHECK_DIRECTION_FOR_XYZ (ix > 0, i - 1);
                          CHECK_DIRECTION_FOR_XYZ (iy > 0, i - data->nx);
                          CHECK_DIRECTION_FOR_XYZ (iz > 0, i - data->nx * data->ny );
                          CHECK_DIRECTION_FOR_XYZ (ix < data->nx - 1 , i + 1);
                          CHECK_DIRECTION_FOR_XYZ (iy < data->ny - 1 , i + data->nx);
                          CHECK_DIRECTION_FOR_XYZ (iz < data->nz - 1 , i + data->nx * data->ny);
                          if (n)
                            {
                              dx[i] = valx / n;
                              dy[i] = valy / n;
                              dz[i] = valz / n;
                            }
                          else
                            {
                              dx[i] = dx_var;
                              dy[i] = dy_var;
                              dz[i] = dz_var;
                            }

                          BOSOUT (section::check_data, level::low)
                            << "Info: new dx in node (" << (ix + 1) << ", "
                              << (iy + 1) << ", " << (iz + 1) << ") = " << dx[i] << " (aver dx = "
                              << dx_var << ")" 
                              << "\n"
                            << "Info: new dy in node (" << (ix + 1) << ", "
                              << (iy + 1) << ", " << (iz + 1) << ") = " << dy[i] << " (aver dy = "
                              << dy_var << ")" 
                              << "\n"
                            << "Info: new dz in node (" << (ix + 1) << ", "
                              << (iy + 1) << ", " << (iz + 1) << ") = " << dz[i] << " (aver dz = "
                              << dz_var << ")" 
                              << bs_end;
                        }
                    }
                }
            }
        }
      else if (data->geom_def_flag == GEOM_FLAG_ZCORN_COORD)
        {
#if 0
          for (iz = 0, i = 0; iz <= nz; ++iz)
            {
              for (iy = 0; iy <= ny; ++iy)
                {
                  for (ix = 0; ix <= nx; ++ix, ++i)
                    {
                      if ((*depth)[i] < 0)
                        {
                          BOSERR (section::check_data, level::error) << "Depth of node (" << (ix + 1) << ", "
                          << (iy + 1) << ", " << (iz + 1) << ") = " << depth[i]
                          << " is out of range" << bs_end;
                          return YS_BAD_VALUE;
                        }
                    }
                }
            }
#endif
        }
      else
        {
          bs_throw_exception ("Unsupported geometry specification");
        }
    }
*/
  
  void hdm::check_arrays_for_inactive_blocks () const
    {
      spv_int actnum_  = data->get_i_array ("ACTNUM");
      spv_float permx_ = data->get_fp_array ("PERMX");
      spv_float permy_ = data->get_fp_array ("PERMY");
      spv_float permz_ = data->get_fp_array ("PERMZ");
      spv_float poro_  = data->get_fp_array ("PORO");
      spv_float ntg_   = data->get_fp_array ("NTG");
      
      t_int* actnum = actnum_->data ();
      const t_float *permx  = permx_->data ();
      const t_float *permy  = permy_->data ();
      const t_float *permz  = permz_->data ();
      const t_float *poro   = poro_->data ();
      const t_float *ntg    = ntg_->data ();

      t_long nx = data->props->get_i("nx");
      t_long ny = data->props->get_i("ny");
      t_long nz = data->props->get_i("nz");
      t_long nb = nx * ny * nz;

      t_long permx_counter = 0;
      t_long permy_counter = 0;
      t_long permz_counter = 0;
      t_long poro_counter = 0;
      t_long ntg_counter = 0;
      for (t_long i = 0; i < nb; ++i)
        {
          if (!actnum[i])
            {
              continue;
            }

          if (actnum[i] && permx[i] < DEFAULT_MINIMAL_PERM && permy[i] < DEFAULT_MINIMAL_PERM && permz[i] < DEFAULT_MINIMAL_PERM)
            {
              t_long iz = i / (nx * ny);
              t_long iy = (i - iz * (nx * ny)) / nx;
              t_long ix = i - iz * (nx * ny) - iy * nx;

              BOSWARN (section::check_data, level::low)
                << "blocks " << i << " [" << (ix + 1) << ", " << (iy + 1) << ", "
                << (iz + 1) << "] will be set inactive because of PERMEABILITY" << bs_end;

              actnum[i] = 0;
              ++permx_counter;
              continue;
            }

          CHECK_FOR_TWO_PARAMS (poro[i], actnum[i], DEFAULT_MINIMAL_PORO, poro_counter);
          if (ntg_->size () != 0)
            {
              CHECK_FOR_TWO_PARAMS (ntg[i], actnum[i], DEFAULT_MINIMAL_NTG, ntg_counter);
            }
        }

      if (permx_counter)
        BOSWARN (section::check_data, level::warning) << permx_counter << " blocks will be set inactive because of PERMX" << bs_end;
      if (permy_counter)
        BOSWARN (section::check_data, level::warning) << permy_counter << " blocks will be set inactive because of PERMY" << bs_end;
      if (permz_counter)
        BOSWARN (section::check_data, level::warning) << permz_counter << " blocks will be set inactive because of PERMZ" << bs_end;
      if (poro_counter)
        BOSWARN (section::check_data, level::warning) << poro_counter << " blocks will be set inactive because of PORO" << bs_end;
      if (ntg_counter)
        BOSWARN (section::check_data, level::warning) << ntg_counter << " blocks will be set inactive because of NTG" << bs_end;
    }
    
 
   
  //bs stuff
  BLUE_SKY_TYPE_STD_CREATE(hdm)
  BLUE_SKY_TYPE_STD_COPY(hdm)

  BLUE_SKY_TYPE_IMPL(hdm, objbase, "hdm", "BOS_Core hdm class", "BOS_Core hdm class")
}//ns bs
