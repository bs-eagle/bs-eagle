#include "bs_bos_core_data_storage_stdafx.h"

#include "bs_misc.h"
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
    this->rock_grid_ = BS_KERNEL.create_object ("rock_grid");
  }

  hdm::hdm(const hdm& src):bs_refcounter (src), lkeeper ("C", LC_ALL)
  {
    *this = src;
  }
 
  void
  hdm::init_fluids(t_int n_pvt_regions, t_int n_scal_regions)
  {
    int n_phases; 
    
    pvt_3p_->init_pvt_arrays (n_pvt_regions, data->props->get_b(L"oil_phase"),
											 data->props->get_b(L"gas_phase"),
                                             data->props->get_b(L"water_phase")
                                             );
    
    n_phases = data->props->get_b(L"oil_phase");
    n_phases += data->props->get_b(L"water_phase");
    n_phases += data->props->get_b(L"gas_phase");
    if (n_phases > 1)
      {
        scal_3p_->init_scal_input_table_arrays (n_scal_regions, data->props->get_b(L"oil_phase"),
                                                                data->props->get_b(L"gas_phase"),
                                                                data->props->get_b(L"water_phase"));
      }
  }
  
  void 
  hdm::init_proc_params()
  {
    data->proc_params->add_property_b (1, L"PRINT_PVTO_TABLE",               L"if this item is set print PVTO to LOG file");
    data->proc_params->add_property_b (1, L"PRINT_PVTW_TABLE",               L"print pvtw interpolation table");
    data->proc_params->add_property_b (1, L"PRINT_PVDG_TABLE",               L"print pvtg interpolation table");
    data->proc_params->add_property_b (0, L"STORE_PANE_FLOW_RATES",          L"");
    data->proc_params->add_property_b (1, L"NEW_TS_SELECTION",               L"use new algorithm for ts selection");
    data->proc_params->add_property_b (1, L"SAVE_STEP_MAPS",                 L"save maps on time step in the rst file");
    data->proc_params->add_property_b (1, L"SAVE_INITIAL_DATA",              L"save initial data in the rst file");
    data->proc_params->add_property_b (1, L"SAVE_RESTART_DATA",              L"save restart information");
    data->proc_params->add_property_b (1, L"NEWTRANS",                       L"correct transmissibility calculation between cells");
    data->proc_params->add_property_b (1, L"CHECK_FOR_MULTY_CONN",           L"check for multi connections in one cell");
    data->proc_params->add_property_b (1, L"SIMPLE_GET_CELL_SOLUTION",       L"save time steps into *tsl file");
    data->proc_params->add_property_b (0, L"USE_TIMESTEP_FILES",             L"enable simple solution restore algorithm");
    data->proc_params->add_property_b (1, L"USE_LOW_SKIN_TRANS_MULT",        L"if this item is set use calculation of trans mults for low skin");
    data->proc_params->add_property_b (1, L"SAVE_MESH_DATA",                 L"save initial data in the rst file");
    data->proc_params->add_property_b (0, L"SAVE_CALC_DATA",                 L"save solution vector at each large time step into .dmp file");
    data->proc_params->add_property_b (0, L"LOAD_CALC_DATA",                 L"load initial guesses for Newton iterations from .dmp file");
    data->proc_params->add_property_b (0, L"SAVE_NORM_DATA",                 L"save cells norms to .cel file");
    data->proc_params->add_property_b (1, L"DENS_JAC_SCALE",                 L"enable densities based multiplier's for gas and oil equation's in Jacobian");
    data->proc_params->add_property_b (0, L"NEWTON_TUNING",                  L"");
    data->proc_params->add_property_b (0, L"SAVE_CROSSFLOW",                 L"");
    data->proc_params->add_property_b (0, L"G_FRACTURE",                     L"");
    data->proc_params->add_property_b (0, L"CLAMP_PRESSURE",                 L"clamp pressure values at each Newton iteration between possible min and max");
    data->proc_params->add_property_b (0, L"G_FRACTURE_FINAL",               L"");
    data->proc_params->add_property_b (1, L"SAVE_CONN_DATA",                 L"save connection rates to rst file");
    data->proc_params->add_property_b (0, L"DISABLE_FRACTURE_CHECK",         L"");
    data->proc_params->add_property_b (0, L"COMPRESS_IRR",                   L"compress irregular matrix");
    data->proc_params->add_property_b (0, L"REMOVE_ISOLATED_REGIONS",        L"make isolated regions with out wells inactive");
    data->proc_params->add_property_b (1, L"SAVE_WRATES_TO_ASCII_FILE",      L"save well rates to ascii file");
    data->proc_params->add_property_b (1, L"CREATE_HDF5_FILE",               L"create hdf5 file");
    data->proc_params->add_property_b (1, L"WRITE_PRESSURE_TO_HDF5",         L"write pressure to hdf5");
    data->proc_params->add_property_b (1, L"WRITE_SATURATION_TO_HDF5",       L"write saturation to hdf5");
    data->proc_params->add_property_b (1, L"FIX_SOIL_BUG",                   L"");
    data->proc_params->add_property_b (0, L"DISABLE_CROSSFLOW",              L"");
    data->proc_params->add_property_b (0, L"DEBUG_EQUIL",                    L"write equil debug info");
    data->proc_params->add_property_b (1, L"WRITE_GAS_OIL_RATIO_TO_HDF5",    L"write gas_oil_ratio to hdf5");
    data->proc_params->add_property_b (0, L"FIX_GRID",                       L"calc dx dy dz as in PPP");
    data->proc_params->add_property_b (1, L"WRITE_WELL_RESULTS_TO_HDF5",     L"save well results to hdf5 file");
    data->proc_params->add_property_b (1, L"WRITE_CONN_RESULTS_TO_HDF5",     L"save connection results to hdf5 file");
    data->proc_params->add_property_b (1, L"WRITE_FIP_RESULTS_TO_HDF5",      L"save fip results to hdf5 file");
    data->proc_params->add_property_b (0, L"WRITE_ECL_CONN_DATA",            L"write connection data in ECL COMPDAT format");
    data->proc_params->add_property_b (0, L"WRITE_INITIAL_DATA_TO_HDF5",     L"save initial data to hdf5 file");
    data->proc_params->add_property_b (0, L"READ_INITIAL_DATA_FROM_HDF5",    L"read initial data from hdf5 file");
    data->proc_params->add_property_b (1, L"WRITE_MESH_DATA_TO_HDF5",        L"save mesh data to hdf5 file");
    data->proc_params->add_property_b (1, L"WRITE_PLANE_FLOW_RATES_TO_HDF5", L"write plane flow rates to hdf5 file");
    data->proc_params->add_property_b (0, L"DISABLE_GRAVITY",                L"if ON, disable gravitation");
    data->proc_params->add_property_b (0, L"PARALLEL_ILU",                   L"if ON, enable parallel ilu preconditioner");
    data->proc_params->add_property_b (0, L"USE_IRR_MATRIX_IN_ILU",          L"if OFF, only diagonal from irregular matrix will be used in ILU decomposition ");
    data->proc_params->add_property_b (0, L"FIX_WELL_DENSITY_BUG",           L"");
    data->proc_params->add_property_b (1, L"WBP",                            L"calc average well pressure only in well blocks");
    data->proc_params->add_property_b (0, L"WBP4",                           L"calc average well pressure only in 4 well neighbours");
    data->proc_params->add_property_b (0, L"WBP5",                           L"calc average well pressure in well blocks and 4 well neighbours");
    data->proc_params->add_property_b (1, L"WBP9",                           L"calc average well pressure in well blocks and 8 well neighbours");
    data->proc_params->add_property_b (0, L"PURE_NEWTON",                    L"disable all hacks ");
    data->proc_params->add_property_b (0, L"WRITE_CNORM_TO_ASCII_FILE",      L"in true C norm of all components will be save to ascii file in format i j k Cw Cg Co");
    data->proc_params->add_property_b (0, L"P_INIT_APPROX",                  L"if true also initialize pressure");
    data->proc_params->add_property_b (0, L"PV_WEIGHTED_PRESSURE",           L"pore volume weighted pressure calculation");
    data->proc_params->add_property_b (0, L"CRS_SCAL_SCALE",                 L"enable 3-point scal scale");
    data->proc_params->add_property_b (1, L"SAVE_PCW",                       L"save calculated PCW to hdf5");
    data->proc_params->add_property_b (0, L"SAVE_PCG",                       L"save calculated PCG to hdf5");
    data->proc_params->add_property_b (0, L"NEGATIVE_BHP",                   L"WELL BHP could be negative if ON");
    data->proc_params->add_property_b (0, L"WRITE_PLANE_FLOW_VOLS_TO_HDF5",  L"write total plane mass flow volumes to hdf5 file");
    data->proc_params->add_property_b (0, L"SET_ACTIVE_FRACTURE_CELLS",      L"set actnum = 1 for cells, where fractures defined");
    data->proc_params->add_property_b (0, L"FRACTURE_HORIZ_WELL",            L"enable fracture for horiz wells, which grows up and down by z-layer ");
    data->proc_params->add_property_b (1, L"WRITE_BHP_0_FOR_CLOSED_WELLS",   L"if true, write well_bhp=0 to h5 file for shutted wells");
    data->proc_params->add_property_b (1, L"MULTI_CONNECTION",               L"if true, allow adding multiple connections to one cell");
    data->proc_params->add_property_b (0, L"FIX_PCW_BUG",                    L"if true, set default psw and pcg values to 0.0");
    data->proc_params->add_property_b (0, L"FIX_WPIMULT_BUG",                L"if true, do not apply WPIMULT to fracture's connections");
    data->proc_params->add_property_b (1, L"FIX_ARITH_BUG",                  L"if true, add brackets around function and around it's arguments in arithmetic");
    data->proc_params->add_property_b (0, L"USE_NITERS_VOLUMETRIC_NORM",     L"use volumetric norm calculation as in Eclipse");
    data->proc_params->add_property_b (0, L"FRACTURE_USE_FABS",              L"use absolute coordinates of points in calculating fracture's connection factors");



    data->proc_params->add_property_i (20,    L"PVT_INTERP_POINTS",        L"number of interpolation points");
    data->proc_params->add_property_i (12,    L"NITERS_NUM",               L"maximum allowed newton iterations number");
    data->proc_params->add_property_i (30,    L"LITERS_NUM",               L"maximum allowed linear iterations number");
    data->proc_params->add_property_i (6,     L"NITERS_INC",               L"number of newton iterations to increment step length");
    data->proc_params->add_property_i (20,    L"NITERS_AMG",               L"number of newton iterations to build full amg setup");
    data->proc_params->add_property_i (0,     L"APPROX_STEPS",             L"number of approximation steps");
    data->proc_params->add_property_i (2,     L"SELECT_SOL_STEPS",         L"number of steps to select newton correction force term");
    data->proc_params->add_property_i (1,     L"LSOLV_TYPE",               L"type of linear solver");
    data->proc_params->add_property_i (10,    L"GMRES_STACK_LEN",          L"number of vectors in GMRES to ortonorm");
    data->proc_params->add_property_i (10,    L"AMG_LITERS_NUM",           L"maximum allowed AMG solver iterations");
    data->proc_params->add_property_i (20,    L"WELL_NITERS_NUM",          L"maximum number of well internal newton iterations");
    data->proc_params->add_property_i (0,     L"PREC_TYPE",                L"type of the preconditioner for linear solver");
    data->proc_params->add_property_i (0,     L"MIN_CELLS_IN_REGION",      L"minimum allowed cells in region for ACTNUM = 1");
    data->proc_params->add_property_i (1,     L"SAVE_STEP_DATA_PERIOD",    L"save step data every n step");
    data->proc_params->add_property_i (10,    L"SAVE_WELL_RESULTS_PERIOD", L"save well results every n step");
    data->proc_params->add_property_i (10,    L"SAVE_FIP_RESULTS_PERIOD",  L"save fip results every n step");
    data->proc_params->add_property_i (20000, L"FRACTURE_SERIES_NUMBER",   L"number of elemements in calculating of series for fracture calculation");
    data->proc_params->add_property_i (3,     L"NEWTON_ITERS_GCONTROL",    L"number of newton iterations to consider group control (equal to NUPCOL, used with GCONINJE and GPMAINT)");
    data->proc_params->add_property_i (0,     L"FRACTURE_ANGLE_AXIS",      L"fracture angle counting: 0-from mesh-based X-axis (by first cell), 1-from X axis, 2-from Y axis (azimut)");
    data->proc_params->add_property_i (2,     L"1PHASE_LSOLV_TYPE",        L"type of linear solver for 1 phase systems");
    data->proc_params->add_property_i (4,     L"1PHASE_PREC_TYPE",         L"type of the preconditioner for linear solver for 1 phase systems");
    data->proc_params->add_property_i (0,     L"WRITE_NORM_TO_HDF5",       L"write norm to hdf5 file flag. 1-write each large step, 2-write at each newton iteration");
    data->proc_params->add_property_i (0,     L"TIMESTEP_ON_NORM_SELECT",  L"timestep depends from cnorm calculated on first iteration");
                                               
                                                      
    data->proc_params->add_property_f (1000,    L"PVT_PRESSURE_MAX",               L"PVT tables maximum pressure");
    data->proc_params->add_property_f (1,       L"PVT_PRESSURE_MIN",               L"PVT tables minimal pressure");
    data->proc_params->add_property_f (400,       L"PVTO_RS_MAX",                    L"PVTO maximum RS for slop table");
    data->proc_params->add_property_f (1,     L"PVTO_RS_MIN",                    L"PVTO minimum RS for slop table");
    data->proc_params->add_property_f (100,     L"TS_MAX",                         L"maximum allowed time step length");
    data->proc_params->add_property_f (5.0e-6,  L"TS_MIN",                         L"minimum allowed time step length");
    data->proc_params->add_property_f (1,       L"TS_FIRST",                       L"first simulation time step length");
    data->proc_params->add_property_f (1.0e-4,  L"LITERS_RESID",                   L"tolerance for linear solver");
    data->proc_params->add_property_f (1.0e-2,  L"NITERS_RESID",                   L"tolerance for newton process");
    data->proc_params->add_property_f (2,       L"TS_INC_MULT",                    L"multiplier for incrementing time step length");
    data->proc_params->add_property_f (0.5,     L"TS_DEC_MULT",                    L"multiplier for decrementing time step length");
    data->proc_params->add_property_f (1.4,     L"TS_OVERDRAFT",                   L"overdraft factor (any time step could be multiplied by this factor to achieve end of report time step)");
    data->proc_params->add_property_f (300,     L"P_CORR_MAX",                     L"maximum allowed newton correction for pressure");
    data->proc_params->add_property_f (1,       L"S_CORR_MAX",                     L"maximum allowed newton correction for saturation");
    data->proc_params->add_property_f (1000,    L"RS_CORR_MAX",                    L"maximum allowed newton correction for gas oil ratio");
    data->proc_params->add_property_f (3000,    L"WELL_P_CORR_MAX",                L"maximum allowed newton correction for well pressure");
    data->proc_params->add_property_f (1.0,     L"WAT_ROW_MULT",                   L"multiplier for water equation in Jacobian");
    data->proc_params->add_property_f (1.0,     L"GAS_ROW_MULT",                   L"multiplier for gas equation in Jacobian");
    data->proc_params->add_property_f (1.0,     L"OIL_ROW_MULT",                   L"multiplier for oil equation in Jacobian");
    data->proc_params->add_property_f (1.0,     L"P_COL_MULT",                     L"multiplier for pressure derivates column in Jacobian");
    data->proc_params->add_property_f (0.5,     L"TS_OMEGA",                       L"omega coef in time step controling");
    data->proc_params->add_property_f (60,      L"TS_DP",                          L"pressure change for time step controling");
    data->proc_params->add_property_f (0.05,    L"TS_DS",                         L"saturation change for time step controling");
    data->proc_params->add_property_f (1e-5,    L"DP_MIN_CHOP",                    L"minimum pressure chop");
    data->proc_params->add_property_f (1e-7,    L"DS_MIN_CHOP",                    L"minimum saturation chop");
    data->proc_params->add_property_f (0.05,    L"LITERS_MAX_RESID",               L"maximum allowed residual");
    data->proc_params->add_property_f (1.0e-1,  L"AMG_RESID",                      L"tolerance for AMG preconditioner");
    data->proc_params->add_property_f (10,      L"TS_DRS",                         L"Rs change for time step controling");
    data->proc_params->add_property_f (10,      L"MAX_NORM_ON_TS",                 L"if norm on time step is greater than this value restart occur");
    data->proc_params->add_property_f (-1,      L"DRSDT",                          L"Maximum rate of increase of solution GOR");
    data->proc_params->add_property_f (1,       L"GAS_NORM_MULT",                  L"multiplier for gas norm");
    data->proc_params->add_property_f (0.99995, L"FRACTURE_LAMBDA",                L"the minimum value for lambda in fracture calculation, than less lambda than less accuracy, must be more than 0.99");
    data->proc_params->add_property_f (0.01,    L"COMP_MASS_BALANCE",              L"Scaled mass balance equation residuals. Equations are scaled by the accumulation terms (component mass dividing timestep size).");
    data->proc_params->add_property_f (0.02,    L"COMP_PHASE_EQUIL",               L"Scaled phase equilibrium relation residuals. These are scaled by the component fugacities in the gas phase");
    data->proc_params->add_property_f (0.0001,  L"COMP_MAX_DP",                    L"The maximum scaled pressure change (absolute pressure change divided by the average reservoir pressure)");
    data->proc_params->add_property_f (0.005,   L"COMP_MAX_DS",                    L"The maximum absolute saturation change");
    data->proc_params->add_property_f (0.001,   L"COMP_MAX_DXCP",                  L"The maximum absolute component mole fraction change");
    data->proc_params->add_property_f (20,      L"COMP_TS_DP",                     L"time step controling pressure change");
    data->proc_params->add_property_f (0.2,     L"COMP_TS_DS",                     L"time step controling saturation change");
    data->proc_params->add_property_f (0.02,    L"COMP_TS_DXCP",                   L"time step controling mole fraction change");
    data->proc_params->add_property_f (200,     L"COMP_MAX_P_CORRECTION",          L"maximum allowed newton correction for pressure");
    data->proc_params->add_property_f (0.5,     L"COMP_MAX_S_CORRECTION",          L"maximum allowed newton correction for saturation");
    data->proc_params->add_property_f (0.2,     L"COMP_MAX_XCP_CORRECTION",        L"maximum allowed newton correction for mole fraction");
    data->proc_params->add_property_f (1.0e-5,  L"MASS_BALANS_ERROR",              L"maximum allowed mass balans error");
    data->proc_params->add_property_f (1,       L"MAX_NORM_ON_FIRST_N",            L"maximum allowed norm on first newton iteration");
    data->proc_params->add_property_f (0.001,   L"P_DIMENSION_LESS_SCALE_FACTOR",  L"scale factor for dimension less pressure");
    data->proc_params->add_property_f (0.3,     L"APPL_CHOP_PERC",                 L"tipicaly (0.25-0.4)");
    data->proc_params->add_property_f (0.2,     L"APPL_CHOP_BND",                  L"tipicaly (0.1-0.25)");
    data->proc_params->add_property_f (0.0,     L"FRACTURE_MIN_DELTA",             L"if dist(A-B) < min_delta => connection will not added to this block, A,B- fracture and block intersection points");
    data->proc_params->add_property_f (100,     L"TIMESTEP_ON_NORM_RESIDUAL_MULT", L"Multiplier for NEWTON_RESIDUAL");
    data->proc_params->add_property_f (0.1,     L"TIMESTEP_ON_NORM_MIN_MULT",      L"minimum multiplier for current timestep");
    data->proc_params->add_property_f (0.5,     L"TIMESTEP_ON_NORM_MAX_MULT",      L"maximum multiplier for current timestep");
    data->proc_params->add_property_f (0.01,    L"TIMESTEP_ON_NORM_RESTART_MULT",  L"restart multiplier for current timestep");
    data->proc_params->add_property_f (100,     L"MAX_PCW_WARNING_VALUE",          L"maximum PCW value to print warning messages in log-file");
  }

  void 
  hdm::init_equil (t_int n_equil_regions)
  {
    int n_phases;
    n_phases = data->props->get_b(L"oil_phase");
    n_phases += data->props->get_b(L"water_phase");
    n_phases += data->props->get_b(L"gas_phase");
    
    this->equil_model_ = BS_KERNEL.create_object ("equil_model_depth"); 
    equil_model_->init_equil_model (n_equil_regions, n_phases);
  }
  
   void
  hdm::init(const std::wstring &model_name_)
  {
    smart_ptr <hdm_iface, true> hdm = this;
    keyword_params kp;
    std::string model_name = wstr2str (model_name_, "ru_RU.CP1251");
    kp.hdm = this;
    if (model_name == std::string(":memory:")) {
        this->well_pool_->open_db (model_name_); //"mem_well_pool.db");
        //this->well_pool_->open_db ("mem2_well_pool.db");
      }
    else
      {
        const std::wstring well_db = model_name_ + L"_well_pool.db";
        data->h5_pool->open_file (model_name + "_rex.h5");
        this->well_pool_->open_db (well_db);
      }
    data->props->add_property_s(str2wstr(model_name.substr(0, model_name.find_last_of("/\\") + 1)), L"model_path", L"model_path");
    km->init(hdm);
    
    
    switch (data->props->get_i(L"mesh"))
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
    
    switch (data->props->get_i(L"init"))
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
    data->props->add_property_s(str2wstr (filename.substr (0, filename.find_last_of("/\\") + 1)), L"model_path", L"model_path");
    
    data->h5_pool->open_file (model_name + ".h5");
    well_pool_->open_db (L":memory:");
    
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

      t_long nx = data->props->get_i(L"nx");
      t_long ny = data->props->get_i(L"ny");
      t_long nz = data->props->get_i(L"nz");
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
