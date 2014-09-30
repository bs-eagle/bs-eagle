/**
 *       \file  calc_model.cpp
 *      \brief  calc_model implementation
 *     \author  Max Nikonov
 *       \date  07.07.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "calc_model.h"
#include "fi_params.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "force_inline.h"
#include "write_time_to_log.h"
#include "vector_assign.h"
#include "memory_macroses.h"
#include "hdm_iface.h"
#include "scal_3p_iface.hpp"
//#include "scale_array_holder.h"
//#include "scal_region_info.h"
//#include "scal_region.h"
//#include "scal_2p_data_holder.h"
#include "rock_grid.h"
#include "norm_calc.h"
#include "data_class.h"
#include "pvt_base.h"
#include "pvt_oil.h"
#include "pvt_dead_oil.h"
#include "pvt_water.h"
#include "pvt_gas.h"
#include "jfunction.h"
#include "rs_mesh_iface.h"
#include "string_formater.h"
#include "init_model_iface.hpp"
#include BS_STOP_PLUGIN_IMPORT ()

//#include "well_results_storage.h"
//#include "fip_results_storage.h"

/////////////////////////////////////////////////////////////////////////////
//! if ORI < NEW than ORI = NEW
#define IF_LE_REPLACE(ORI,NEW,C,V) if ((ORI) > (NEW)) {(ORI) = (NEW);(C) = (V);}

/**
 * \todo describe
 * */
template <typename item_t>
BS_FORCE_INLINE void
CHECK_VALUE (item_t X, item_t DX, item_t MIN, item_t MAX, item_t DMAX, item_t &M)
{
  if (fabs (DX) > DMAX)
    {
      M = (std::min<item_t>) (M, DMAX / fabs (DX));
    }
  if ((X + DX) < MIN)
    {
      M = (std::min<item_t>) (M, (MIN - X) / DX);
    }
  if ((X + DX) > MAX)
    {
      M = (std::min<item_t>) (M, (MAX - X) / DX);
    }
}

namespace blue_sky
{

  /**
   * \brief  'default' ctor for calc_model
   * \param  param additional params for calc_model
   * */
  calc_model::calc_model (bs_type_ctor_param /*param*/)
      : bs_refcounter(), bs_node(bs_node::create_node())
      , pvt_regions (BS_KERNEL.create_object (v_long::bs_type ()))
      , sat_regions (BS_KERNEL.create_object (v_long::bs_type ()))
      , fip_regions (BS_KERNEL.create_object (v_long::bs_type ()))
      , rock_regions (BS_KERNEL.create_object (v_long::bs_type ()))
      , pressure (BS_KERNEL.create_object (v_double::bs_type ()))
      , saturation_3p (BS_KERNEL.create_object (v_double::bs_type ()))
      , gas_oil_ratio (BS_KERNEL.create_object (v_double::bs_type ()))
      //, well_res (BS_KERNEL.create_object ("well_results_storage"))
      //, fip_res (BS_KERNEL.create_object ("fip_results_storage"))
  {
    init();
  }

  /**
   * \brief  copy-ctor for calc_model
   * \param  src calc_model instance
   * */
  calc_model::calc_model (const this_t & /*src*/)
      : bs_refcounter(), bs_node(bs_node::create_node())
      , pvt_regions (BS_KERNEL.create_object (v_long::bs_type ()))
      , sat_regions (BS_KERNEL.create_object (v_long::bs_type ()))
      , fip_regions (BS_KERNEL.create_object (v_long::bs_type ()))
      , rock_regions (BS_KERNEL.create_object (v_long::bs_type ()))
      , pressure (BS_KERNEL.create_object (v_double::bs_type ()))
      , saturation_3p (BS_KERNEL.create_object (v_double::bs_type ()))
      , gas_oil_ratio (BS_KERNEL.create_object (v_double::bs_type ()))
      //, well_res (BS_KERNEL.create_object ("well_results_storage"))
      //, fip_res (BS_KERNEL.create_object ("fip_results_storage"))
  {
    //*this = src;
  }

  calc_model::~calc_model()
  {
  }

  void 
  calc_model::init()
  {
    n_comps = 0;
    n_phases = 0;
    phases = FI_PHASE_NULL;
    n_HCcomps = 0;
    n_HCphases = 0;
    n_pri = 0;
    n_sec = 0;
    n_vars = 0;

    assign (phase_d, 0);
    assign (sat_d, 0);

    n_pvt_regions = 0;
    n_sat_regions = 0;
    n_fip_regions = 0;

    rpo_model = RPO_DEFAULT_MODEL;

    last_c_norm = 0;
    approx_flag = 0;

    rock_grid_prop = BS_KERNEL.create_object (rock_grid::bs_type ());
    ts_params = BS_KERNEL.create_object (fi_params::bs_type ());
  }

  int 
  calc_model::init_main_arrays (const BS_SP (hdm_iface) hdm)
  {
#ifdef _DEBUG
    BOSOUT (section::init_data, level::debug) << "FI DEBUG: begin of init_main_arrays method" << bs_end;
#endif // _DEBUG
    write_time_to_log init_time ("Main arrays initialization", "");

    // set up phases
    n_phases = 0;
    phases = 0;

    BS_SP (idata)            input_data = hdm->get_data ();
    BS_SP (pvt_3p_iface)     pvt_prop_  = hdm->get_pvt ();
    BS_SP (scal_3p_iface)    scal_prop_ = hdm->get_scal ();
    BS_SP (init_model_iface) init_model = hdm->get_init_model ();
    BS_SP (rs_mesh_iface)    mesh       = hdm->get_mesh ();
    
    if (input_data->props->get_b (L"water_phase"))
      {
        phase_d[0] = n_phases++;
        phases |= 1 << FI_PHASE_WATER;
      }
    else 
      phase_d[0] = -1;

    if (input_data->props->get_b (L"gas_phase"))
      {
        phase_d[1] = n_phases++;
        phases |= 1 << FI_PHASE_GAS;
      }
    else 
      phase_d[1] = -1;

    if (input_data->props->get_b (L"oil_phase"))
      {
        phase_d[2] = n_phases++;
        phases |= 1 << FI_PHASE_OIL;
      }
    else 
      phase_d[2] = -1;

    if (n_phases > 1)
      n_sec_vars = 1;

    for (size_t i = 0, sat_counter = 0; i < FI_PHASE_TOT; ++i)
      {
        if ((this->phases & (1 << i)) && (sat_counter < this->n_phases - 1))
          this->sat_d[i] = sat_counter++;
        else
          this->sat_d[i] = -1;
      }

    if (this->n_phases < 1)
      {
        BOSERR (section::init_data, level::warning)
        << "phases for full implicit not specified. Use keywords: OIL, WATER, GAS, ..."
        << bs_end;
        return -3;
      }

    // set number of components and variables
    this->n_comps = this->n_phases;
    this->n_pri = this->n_phases;
    this->n_sec = 0;
    this->n_vars = this->n_pri + this->n_sec;

    int is_w = FI_CHK_WATER (this->phases);
    int is_g = FI_CHK_GAS (this->phases);
    int is_o = FI_CHK_OIL (this->phases);

    // fill indexes values
    if (this->n_phases == 1)
      {
        if (is_w)
          b_w_p = 0;
        else if (is_g)
          b_g_p = 0;
        else if (is_o)
          b_o_p = 0;
      }
    else if (this->n_phases == 2)
      {
        if (is_w && is_o)
          {
            b_w_w = 0;
            b_w_p = 1;
            b_o_w = 2;
            b_o_p = 3;
          }
        else if (is_g && is_o)
          {
            b_g_g = 0;
            b_g_p = 1;
            b_o_g = 2;
            b_o_p = 3;
          }
      }
    else
      {
        b_w_w = 0;
        b_w_g = 1;
        b_w_p = 2;
        b_g_w = 3;
        b_g_g = 4;
        b_g_p = 5;
        b_o_w = 6;
        b_o_g = 7;
        b_o_p = 8;
      }

    main_variable.resize(mesh->get_n_active_elements());
    old_data_.main_var.resize(mesh->get_n_active_elements());
    prev_niter_data_.main_var.resize(mesh->get_n_active_elements());

//#ifndef _DEBUG
//    data.reserve (mesh->get_n_active_elements ());
//    for (size_t i = 0, cnt = mesh->get_n_active_elements (); i < cnt; ++i)
//      {
//        data.push_back (data_t (n_phases));
//      }
//#else
    data.resize (mesh->get_n_active_elements ());
//#endif

    fip_regions->init (mesh->get_n_active_elements());
    sat_regions->init (mesh->get_n_active_elements());
    pvt_regions->init (mesh->get_n_active_elements());

    pressure->init (mesh->get_n_active_elements(), 0.0);
    old_data_.pressure->init (mesh->get_n_active_elements(), 0.0);
    prev_niter_data_.pressure->init (mesh->get_n_active_elements(), 0.0);

    if (this->n_phases > 1)
      {
        t_long count = mesh->get_n_active_elements() * this->n_phases;
        saturation_3p->init (count, 0.0);
        old_data_.saturation_3p->init (count, 0.0);
        prev_niter_data_.saturation_3p->init (count, 0.0);

        //initialize gas-oil ratio array
        if (FI_CHK_OIL_GAS(this->phases))
          {
            t_long count = mesh->get_n_active_elements();
            gas_oil_ratio->init (count, 0.0);
            old_data_.gas_oil_ratio->init (count, 0.0);
            prev_niter_data_.gas_oil_ratio->init (count, 0.0);
          }
      }

    if (input_data->props->get_i (L"rock_region") > 0)
      {
        rock_regions->init (mesh->get_n_active_elements());

        if (input_data->contains_i_array ("ROCKNUM"))
          {
            convert_arrays(mesh->get_n_active_elements (), mesh->get_int_to_ext (), *rock_regions, input_data->get_i_array ("ROCKNUM"));
            FI_DECR_ARRAY (*rock_regions, 0, mesh->get_n_active_elements(), 1);
          }
        else
          {
            rock_regions->init (mesh->get_n_active_elements ());
          }
        this->rocktab = input_data->rocktab;
      }

    // initialize fip regions
    this->n_fip_regions = input_data->props->get_i (L"fip_region");
    if (!this->n_fip_regions)
      {
        bs_throw_exception ("Number of fip regions should be not 0");
      }
    if (input_data->contains_i_array ("FIPNUM"))
      {
        convert_arrays(mesh->get_n_active_elements (), mesh->get_int_to_ext (), *fip_regions, input_data->get_i_array ("FIPNUM"));
        FI_DECR_ARRAY (*fip_regions, 0, mesh->get_n_active_elements(), 1);
      }
    else
      {
        FI_FILL_ARRAY (*fip_regions, 0, mesh->get_n_active_elements(), 0);
      }

    // initialize pvt regions
    this->n_pvt_regions = input_data->props->get_i (L"pvt_region");
    if (!this->n_pvt_regions)
      {
        bs_throw_exception ("Number of pvt regions should be not 0");
      }
    if (input_data->contains_i_array ("PVTNUM"))
      {
        convert_arrays(mesh->get_n_active_elements (), mesh->get_int_to_ext (), *pvt_regions, input_data->get_i_array ("PVTNUM"));
        FI_DECR_ARRAY (*pvt_regions, 0, mesh->get_n_active_elements(), 1);
      }
    else
      {
        FI_FILL_ARRAY (*pvt_regions, 0, mesh->get_n_active_elements(), 0);
      }

    // initialize sat regions
    this->n_sat_regions = input_data->props->get_i (L"sat_region");
    if (!this->n_sat_regions)
      {
        bs_throw_exception ("Number of sat regions should be not 0");
      }
    if (input_data->contains_i_array ("SATNUM"))
      {
        convert_arrays(mesh->get_n_active_elements (), mesh->get_int_to_ext (), *sat_regions, input_data->get_i_array ("SATNUM"));
        FI_DECR_ARRAY (*sat_regions, 0, mesh->get_n_active_elements(), 1);
      }
    else
      {
        FI_FILL_ARRAY (*sat_regions, 0, mesh->get_n_active_elements(), 0);
      }

    // initialize rpo_model
    this->rpo_model = (RPO_MODEL_ENUM)input_data->props->get_i (L"rpo_model");

    // allocate rock grid data storage
    this->rock_grid_prop->init(input_data, mesh->get_n_active_elements(), this->n_pvt_regions);

    //initialize pvt
    this->init_pvt_arrays(pvt_prop_);

#if 1
    this->init_scal_arrays (scal_prop_, input_data);
#else
    // initialize scale arrays
    scal_prop = scal_prop_;
    const BS_SP (scale_array_holder_iface) &gas_scale_ = scal_prop->get_gas_scale ();
    const BS_SP (scale_array_holder_iface) &water_scale_ = scal_prop->get_water_scale ();

    if (input_data->is_set ("SOGCR"))   gas_scale_->set (socr, "SOGCR", input_data->get_fp_array ("SOGCR"));
    if (input_data->is_set ("SGCR"))    gas_scale_->set (scr, "SGCR", input_data->get_fp_array ("SGCR"));
    if (input_data->is_set ("SGU"))     gas_scale_->set (su, "SGU", input_data->get_fp_array ("SGU"));
    if (input_data->is_set ("SGL"))     gas_scale_->set (sl, "SGL", input_data->get_fp_array ("SGL"));

    if (input_data->is_set ("SOWCR"))   water_scale_->set (socr, "SOWCR", input_data->get_fp_array ("SOWCR"));
    if (input_data->is_set ("SWCR"))    water_scale_->set (scr, "SWCR", input_data->get_fp_array ("SWCR"));
    if (input_data->is_set ("SWU"))     water_scale_->set (su, "SWU", input_data->get_fp_array ("SWU"));
    if (input_data->is_set ("SWL"))     water_scale_->set (sl, "SWL", input_data->get_fp_array ("SWL"));
    if (input_data->is_set ("PCW"))     water_scale_->set (pcp, "PCW", input_data->get_fp_array ("PCW"));

    scal_prop->set_water_jfunction (BS_KERNEL.create_object (jfunction::bs_type ()));
    scal_prop->set_gas_jfunction (BS_KERNEL.create_object (jfunction::bs_type ()));
    scal_prop->init (is_water (), is_gas (), is_oil (), phase_d, sat_d, rpo_model);
    scal_prop->update_gas_data ();
#endif 

    // initialize rock grid data
    this->rock_grid_prop->init_data(mesh->get_n_active_elements (), mesh->get_int_to_ext (), input_data);

    init_model->init (BS_SP (calc_model) (this), input_data, mesh);

#if 0
    tools::save_seq_vector ("pressure_bs.txt").save (pressure);
    if (n_phases > 1)
      tools::save_seq_vector ("saturation_bs.txt").save (saturation_3p);
    if (FI_CHK_OIL_GAS (phases))
      tools::save_seq_vector ("rs_bs.txt").save (gas_oil_ratio);

    tools::save_seq_vector ("cells_bs.txt").save (mesh->get_int_to_ext ());
#endif

    return 0;
  }

  int 
  calc_model::init_calcul_arrays (const sp_idata_t &input_data, const sp_mesh_iface_t &mesh)
  {
    write_time_to_log init_time ("Calc arrays initialization", "");
    
    if (this->ts_params->get_bool(fi_params::STORE_PANE_FLOW_RATES))
      {
        //TODO: n_planes
        this->plane_flow_rate.resize(mesh->get_n_active_elements() * this->n_phases);
        this->full_step_plane_flow_rate.resize(mesh->get_n_active_elements() * this->n_phases);
      }

    init_boundary_connections (input_data, mesh);

    return 0;
  }

  void
  calc_model::init_boundary_connections (const sp_idata_t & /*input_data*/, const sp_mesh_iface_t & /*mesh*/)
  {
    // initialize boundary connections
    // for (int i = 0; i < get_mesh()->n_boundary_connections; ++i) {
//       int i_block;
//       i_block = get_mesh()->boundary_connection_element[i];
//       if (TEST_BITS_UINT (get_mesh()->boundary_connection_type_p, i, 1) == FI_MESH_BCONN_FIX_PRESSURE)
//         {
//           bconn_pressure[i] = pressure[i_block];
//           if (n_phases == 2)
//             {
//               bconn_saturation[i] = saturation[i_block];
//             }
//           else if (n_phases == 3)
//             {
//               bconn_saturation[i * (n_phases - 1) + 0] = saturation[i_block * (n_phases - 1) + 0];
//               bconn_saturation[i * (n_phases - 1) + 1] = saturation[i_block * (n_phases - 1) + 1];
//             }
//           if (FI_CHK_OIL_GAS (phases))
//             {
//               bconn_gor[i] = gas_oil_ratio[i_block];
//               bconn_mainvar[i] = main_variable[i_block];
//             }

//           if (scal_prop.sowcr)
//             scal_prop.bconn_sowcr[i] = scal_prop.sowcr[i_block];
//           if (scal_prop.swcr)
//             scal_prop.bconn_swcr[i] = scal_prop.swcr[i_block];
//           if (scal_prop.swl)
//             scal_prop.bconn_swl[i] = scal_prop.swl[i_block];
//           if (scal_prop.swu)
//             scal_prop.bconn_swu[i] = scal_prop.swu[i_block];
//           if (scal_prop.pcw)
//             scal_prop.bconn_pcw[i] = scal_prop.pcw[i_block];
//           if (scal_prop.sogcr)
//             scal_prop.bconn_sogcr[i] = scal_prop.sogcr[i_block];
//           if (scal_prop.sgcr)
//             scal_prop.bconn_sgcr[i] = scal_prop.sgcr[i_block];
//           if (scal_prop.sgl)
//             scal_prop.bconn_sgl[i] = scal_prop.sgl[i_block];
//           if (scal_prop.sgu)
//             scal_prop.bconn_sgu[i] = scal_prop.sgu[i_block];
//         }
//     }
  }

  const calc_model &
  calc_model::operator=(const calc_model &src)
  {
    n_comps = src.n_comps;
    n_phases = src.n_phases;
    phases = src.phases;
    n_HCcomps = src.n_HCcomps;
    n_HCphases = src.n_HCphases;
    n_pri = src.n_pri;
    n_sec = src.n_sec;
    n_vars = src.n_vars;

    std::copy (src.phase_d.begin (), src.phase_d.end (), phase_d.begin ());
    std::copy (src.sat_d.begin (), src.sat_d.end (), sat_d.begin ());

    n_pvt_regions = src.n_pvt_regions;
    n_sat_regions = src.n_sat_regions;
    n_fip_regions = n_fip_regions;

    rpo_model = src.rpo_model;

    last_c_norm = src.last_c_norm;
    approx_flag = src.approx_flag;

    pvt_regions->assign (src.pvt_regions->begin (), src.pvt_regions->end ());
    sat_regions->assign (src.sat_regions->begin (), src.sat_regions->end ());
    fip_regions->assign (src.fip_regions->begin (), src.fip_regions->end ());
    rock_regions->assign (src.rock_regions->begin (), src.rock_regions->end ());

    plane_flow_rate.assign(src.plane_flow_rate.begin(),src.plane_flow_rate.end());
    full_step_plane_flow_rate.assign(src.full_step_plane_flow_rate.begin(),src.full_step_plane_flow_rate.end());
    //truns_mult.assign(src.truns_mult.begin(),src.truns_mult.end());
    //p_deriv_truns_mult.assign(src.p_deriv_truns_mult.begin(),src.p_deriv_truns_mult.end());

    //scal_prop = give_kernel::Instance().create_object_copy(src.scal_prop);
    rock_grid_prop = give_kernel::Instance().create_object_copy(src.rock_grid_prop);

    ts_params = give_kernel::Instance().create_object_copy(src.ts_params);

    mat = give_kernel::Instance().create_object_copy(src.mat);
    return *this;
  }

  /**
   * \class pvt_helper
   * \todo  describe
   * */
  struct pvt_helper
    {
      template <typename pvt_array_t, typename pvt_vector_t>
      static void
      set_array (const pvt_array_t &pvt, const pvt_vector_t &v)
      {
        for (size_t i = 0; i < v.size(); ++i)
          {
            const typename pvt_array_t::value_type &pvt__(pvt[i]);

            const typename pvt_vector_t::value_type &p = v[i];
            pvt__->insert_vector(*p.main_data_);
            if (p.has_density_)
              {
                pvt__->set_density (p.density_, p.molar_density_);
              }
          }
      }
    };

  void 
  calc_model::init_pvt_arrays (BS_SP (pvt_3p_iface) pvt_prop_)
  {
    pvt_prop = pvt_prop_;
    
    float atm_p     = this->internal_constants.atmospheric_pressure;
    float min_p     = this->ts_params->get_float (fi_params::PVT_PRESSURE_RANGE_MIN);
    float max_p     = this->ts_params->get_float (fi_params::PVT_PRESSURE_RANGE_MAX);
    int n_intervals = this->ts_params->get_int (fi_params::PVT_INTERP_POINTS);

    pvt_prop->init_pvt_calc_data (atm_p, min_p, max_p, n_intervals);

    pvt_oil_array = pvt_prop->get_pvt_oil_array ();
    pvt_gas_array = pvt_prop->get_pvt_gas_array ();
    pvt_water_array = pvt_prop->get_pvt_water_array ();
  }

  void
  calc_model::init_scal_arrays (BS_SP (scal_3p_iface) scal_prop_, const BS_SP (idata) input_data)
  {    
    // initialize scale arrays
    scal_prop = scal_prop_;
    scal_prop->init_scal_calc_data ();

    const BS_SP (scale_array_holder_iface) &gas_scale_ = scal_prop->get_gas_scale ();
    const BS_SP (scale_array_holder_iface) &water_scale_ = scal_prop->get_water_scale ();

    if (input_data->is_set ("SOGCR"))   gas_scale_->set (socr, L"SOGCR", input_data->get_fp_array ("SOGCR"));
    if (input_data->is_set ("SGCR"))    gas_scale_->set (scr, L"SGCR", input_data->get_fp_array ("SGCR"));
    if (input_data->is_set ("SGU"))     gas_scale_->set (su, L"SGU", input_data->get_fp_array ("SGU"));
    if (input_data->is_set ("SGL"))     gas_scale_->set (sl, L"SGL", input_data->get_fp_array ("SGL"));

    if (input_data->is_set ("SOWCR"))   water_scale_->set (socr, L"SOWCR", input_data->get_fp_array ("SOWCR"));
    if (input_data->is_set ("SWCR"))    water_scale_->set (scr, L"SWCR", input_data->get_fp_array ("SWCR"));
    if (input_data->is_set ("SWU"))     water_scale_->set (su, L"SWU", input_data->get_fp_array ("SWU"));
    if (input_data->is_set ("SWL"))     water_scale_->set (sl, L"SWL", input_data->get_fp_array ("SWL"));
    if (input_data->is_set ("PCW"))     water_scale_->set (pcp, L"PCW", input_data->get_fp_array ("PCW"));

  }

  calc_model::item_t
  calc_model::get_initial_rho (item_t height) const
    {
      BS_ASSERT (!pvt_oil_array.empty ());

      const double gravity = internal_constants.gravity_constant;
      const double density = 1.5 * pvt_oil_array.front ()->get_surface_density ();

      return density * gravity * height;
    }

  void
  calc_model::update_min_pressure_range (item_t min_range)
  {
    this->ts_params->set_float (fi_params::PVT_PRESSURE_RANGE_MIN, min_range);
  }

  void
  calc_model::update_max_pressure_range (item_t max_range)
  {
    this->ts_params->set_float (fi_params::PVT_PRESSURE_RANGE_MAX, max_range);
  }

  restore_solution_return_type
  calc_model::restore_solution (const sp_mesh_iface_t &mesh, const spv_double &solution, const spv_double &sec_solution)
  {
    return apply_newton_correction (1.0, 0, mesh, solution, sec_solution);
  }

  restore_solution_return_type
  calc_model::apply_newton_correction (item_t mult, index_t istart_line_search, const sp_mesh_iface_t &mesh, const spv_double &solution, const spv_double &sec_solution)
  {
    //int ret_code = 0;
    static double mult_out = 1.0;   // TODO: WTF

    if (!istart_line_search) // calculate multiplier first time only
      {
        mult_out = new_simple_get_cell_solution_mult_2 (mesh, solution);
        //mult_out = new_simple_get_cell_solution_mult (mesh, solution, sec_solution);

        //if (fabs (mult_out - mult_out_x) > 0.0)
        //  {
        //    BOSOUT (section::iters, level::debug) << "MULT_OUT_DIFF: " << mult_out << " - " << mult_out_x << ": " << mult_out - mult_out_x << bs_end;
        //  }

        if (mult_out < 0.2)
          {
            mult_out = 0.2;
          }
      }
    else
      {
        // calculate total mult from internal and external mult
        mult_out *= mult;
      }

    BOSOUT (section::iters, level::low) << "Solution MULT: " << mult_out << bs_end;
    // restore solution
    if (new_simple_get_cell_solution (mult_out, istart_line_search, mesh, solution, sec_solution))
      {
        BOSERR (section::iters, level::error) << "apply_newton_correction: Couldn't apply newton correction for cells" << bs_end;
        return SMALL_TIME_STEP_FAIL;
      }

    return SMALL_TIME_STEP_OK;
  }

//! if ORI < NEW than ORI = NEW
#define IF_LE_REPLACE(ORI,NEW,C,V) if ((ORI) > (NEW)) {(ORI) = (NEW);(C) = (V);}

  calc_model::item_t
  calc_model::new_simple_get_cell_solution_mult_2 (const sp_mesh_iface_t &msh, const spv_double &x_sol_) const
  {
    BS_ASSERT (msh);
    BS_ASSERT (x_sol_);
    BS_ASSERT (!pressure->empty ());
    if (n_phases > 1)
      {
        BS_ASSERT (!saturation_3p->empty ());
      }

    item_t mult                   = 1.0;
    const item_t max_p            = (item_t) ts_params->get_float (fi_params::MAX_P_CORRECTION);
    const item_t max_s            = (item_t) ts_params->get_float (fi_params::MAX_S_CORRECTION);
    const item_t minimal_pressure = (item_t) ts_params->get_float (fi_params::PVT_PRESSURE_RANGE_MIN);
    const item_t maximal_pressure = (item_t) ts_params->get_float (fi_params::PVT_PRESSURE_RANGE_MAX);
    const item_t max_rs           = (item_t) ts_params->get_float (fi_params::MAX_RS_CORRECTION);
    const bool is_w               = FI_CHK_WATER (phases);
    const bool is_oil_gas         = FI_CHK_OIL_GAS (phases);
    const int d_w                 = phase_d[FI_PHASE_WATER];
    const int d_g                 = phase_d[FI_PHASE_GAS];
    int condition                 = 0;
    index_t n                     = msh->get_n_active_elements ();
    const item_t *x_sol           = x_sol_->data ();
    item_t d                      = 0.0;
    item_t t_mult                 = 0.0;

    // loop through all cells
    for (index_t i = 0; i < n; ++i)
      {
        // check water saturation
        if (is_w && n_phases > 1)
          {
            index_t jj = i * n_phases + d_w;

            if (fabs (x_sol[jj]) > max_s)
              {
                t_mult = max_s / fabs (x_sol[jj]);
                IF_LE_REPLACE (mult, t_mult, condition, 1);
              }
          }
#if 1
        if (is_oil_gas)
          {
            index_t jj = i * n_phases + d_g;
            // check gas saturation
            if (main_variable[i] == FI_SG_VAR)
              {
                if (fabs (x_sol[jj]) > max_s)
                  {
                    t_mult = max_s / fabs (x_sol[jj]);
                    IF_LE_REPLACE (mult, t_mult, condition, 2);
                  }
              }
            // check gas oil ratio
            else
              {
                if (fabs (x_sol[jj]) > max_rs)
                  {
                    t_mult = max_rs / fabs (x_sol[jj]);
                    IF_LE_REPLACE (mult, t_mult, condition, 3);
                  }
                if ((*gas_oil_ratio)[i] + x_sol[jj] < 0)
                  {
                    t_mult = (*gas_oil_ratio)[i] / (-x_sol[jj]);
                    IF_LE_REPLACE (mult, t_mult, condition, 4);
                  }
              }
          }
#endif
        // check pressure
        index_t jj = (i + 1) * n_phases - 1;
        d = (*pressure)[i] + x_sol[jj];
        if (d < minimal_pressure)
          {
            t_mult = ((*pressure)[i] - minimal_pressure) / (-x_sol[jj]);
            IF_LE_REPLACE (mult, t_mult, condition, 5);
          }
        if (maximal_pressure < d)
          {
            t_mult = ((*pressure)[i] - maximal_pressure + 5) / (-x_sol[jj]);
            IF_LE_REPLACE (mult, t_mult, condition, 7);
          }
        if (fabs (x_sol[jj]) > max_p)
          {
            t_mult = max_p / fabs (x_sol[jj]);
            IF_LE_REPLACE (mult, t_mult, condition, 6);
          }
      }
#ifdef _DEBUG
    //BOSOUT (section::iters, level::debug) << boost::format ("MULT CONDITION %d") % condition << bs_end;
#endif
    return mult;
  }

  calc_model::item_t
  calc_model::new_simple_get_cell_solution_mult (const sp_mesh_iface_t &msh,
      const spv_double &x_sol_, const spv_double &sec_sol_)
  {
    BS_ASSERT (msh);
    BS_ASSERT (x_sol_);
    BS_ASSERT (sec_sol_);
    BS_ASSERT (pressure->size ());
    if (n_phases > 1)
      {
        BS_ASSERT (saturation_3p->size ());
      }

    const item_t max_p            = (item_t) ts_params->get_float (fi_params::MAX_P_CORRECTION);
    const item_t max_s            = (item_t) ts_params->get_float (fi_params::MAX_S_CORRECTION);
    const item_t minimal_pressure = (item_t) ts_params->get_float (fi_params::PVT_PRESSURE_RANGE_MIN);
    const item_t maximal_pressure = (item_t) ts_params->get_float (fi_params::PVT_PRESSURE_RANGE_MAX);
    const item_t max_rs           = (item_t) ts_params->get_float (fi_params::MAX_RS_CORRECTION);

    const index_t is_w            = (index_t) FI_CHK_WATER (phases);
    const index_t is_o            = (index_t) FI_CHK_OIL (phases);
    const index_t is_g            = (index_t) FI_CHK_GAS (phases);

    const index_t d_w             = (index_t) phase_d[FI_PHASE_WATER];
    const index_t d_g             = (index_t) phase_d[FI_PHASE_GAS];
    const index_t d_o             = (index_t) phase_d[FI_PHASE_OIL];

    item_t mult                   = 1.0;
    index_t n                     = msh->get_n_active_elements ();
    const item_t *x_sol           = &(*x_sol_)[0];
    const item_t *sec_sol         = &(*sec_sol_)[0];

    index_t n_left = 0;
    index_t n_right = n;

    // calc only internal cells
    for (index_t i = n_left; i < n_right; ++i)
      {
        item_t dsw = 0;
        item_t dpo = 0;
        item_t dsg = 0;
        item_t dro = 0;
        item_t dso = 0;

        index_t jj = (i - n_left);
        if (n_phases == 3)
          {
            dsw = sec_sol[jj];
            dpo = x_sol[jj * n_phases + p3_po];
            dsg = 0;
            dro = 0;
            dso = 0;
            if (main_variable[i] == FI_SG_VAR)
              {
                dsg = x_sol[jj * n_phases + p3_sg];
                dso = x_sol[jj * n_phases + p3_so];
              }
            else if (main_variable[i] == FI_RO_VAR)
              {
                dro = x_sol[jj * n_phases + p3_sg];
                dso = x_sol[jj * n_phases + p3_so];
              }
            CHECK_VALUE <item_t> ((*saturation_3p)[i * n_phases + d_w], dsw, -0.1f, 1.1f, max_s, mult);
            CHECK_VALUE <item_t> ((*pressure)[i], dpo, minimal_pressure, maximal_pressure, max_p, mult);
            CHECK_VALUE <item_t> ((*saturation_3p)[i * n_phases + d_g], dsg, -0.1f, 1.1f, max_s, mult);
            CHECK_VALUE <item_t> ((*saturation_3p)[i * n_phases + d_o], dso, -0.1f, 1.1f, max_s, mult);
            CHECK_VALUE <item_t> ((*gas_oil_ratio)[i], dro, 0, 10000000, max_rs, mult);
          }
        else if (is_w && is_o)
          {
            dsw = sec_sol[jj];
            dso = x_sol[jj * n_phases + p2ow_so];
            dpo = x_sol[jj * n_phases + p2ow_po];

            CHECK_VALUE <item_t> ((*saturation_3p)[i * n_phases + d_w], dsw, -0.1f, 1.1f, max_s, mult);
            CHECK_VALUE <item_t> ((*saturation_3p)[i * n_phases + d_o], dso, -0.1f, 1.1f, max_s, mult);
            CHECK_VALUE <item_t> ((*pressure)[i], dpo, minimal_pressure, maximal_pressure, max_p, mult);
          }
        else if (is_g && is_o)
          {
            dpo = x_sol[jj * n_phases + p2og_po];
            dsg = 0;
            dro = 0;
            dso = 0;
            if (main_variable[i] == FI_SG_VAR)
              {
                dso = sec_sol[jj];
                dsg = x_sol[jj * n_phases + p2og_sg];
              }
            else if (main_variable[i] == FI_RO_VAR)
              {
                dso = sec_sol[jj];
                dro = x_sol[jj * n_phases + p2og_sg];
              }
            CHECK_VALUE <item_t> ((*pressure)[i], dpo, minimal_pressure, maximal_pressure, max_p, mult);
            CHECK_VALUE <item_t> ((*saturation_3p)[i * n_phases + d_g], dsg, -0.1f, 1.1f, max_s, mult);
            CHECK_VALUE <item_t> ((*saturation_3p)[i * n_phases + d_o], dso, -0.1f, 1.1f, max_s, mult);
            CHECK_VALUE <item_t> ((*gas_oil_ratio)[i], dro, 0, 10000000, max_rs, mult);
          }
        else
          {
            dpo = x_sol[jj];

            CHECK_VALUE <item_t> ((*pressure)[i], dpo, minimal_pressure, maximal_pressure, max_p, mult);
          }
      }
    return mult;
    return 1.0; 
  }

  int
  calc_model::new_simple_get_cell_solution (const double mult, int istart_linear_search,
      const sp_mesh_iface_t &msh,
      const spv_double &x_sol_,
      const spv_double &sec_sol_)
  {
    BS_ASSERT (msh);
    BS_ASSERT (x_sol_);
    BS_ASSERT (sec_sol_);
    BS_ASSERT (pressure->size ());
    if (n_phases > 1)
      {
        BS_ASSERT (saturation_3p->size ());
      }

    const int is_w                = FI_CHK_WATER (phases);
    const int is_o                = FI_CHK_OIL (phases);
    const int is_g                = FI_CHK_GAS (phases);
    //const int is_oil_gas          = FI_CHK_OIL_GAS (phases);

    const int d_w                 = phase_d[FI_PHASE_WATER];
    const int d_g                 = phase_d[FI_PHASE_GAS];
    const int d_o                 = phase_d[FI_PHASE_OIL];

    const double minimal_pressure = ts_params->get_float_d (fi_params::PVT_PRESSURE_RANGE_MIN);
    const double maximal_pressure = ts_params->get_float_d (fi_params::PVT_PRESSURE_RANGE_MAX);
    const double max_p            = ts_params->get_float_d (fi_params::MAX_P_CORRECTION);
    //const double max_s            = ts_params->get_float_d (fi_params::MAX_S_CORRECTION);
    //const double max_rs           = ts_params->get_float_d (fi_params::MAX_RS_CORRECTION);
    //const int clamp_pres          = ts_params->get_bool_d  (fi_params::CLAMP_PRESSURE);
    //const double appl_perc        = ts_params->get_float_d (fi_params::APPL_CHOP_PERC);
    //const double appl_bnd         = ts_params->get_float_d (fi_params::APPL_CHOP_BND);

    item_t *sol_pressure              = 0;
    item_t *sol_saturation_3p         = 0;
    item_t *sol_gas_oil_ratio         = 0;
    main_var_type *sol_main_variable  = 0;

    if (!istart_linear_search)
      {
        sol_pressure              = &(*pressure)[0];
        sol_main_variable         = &main_variable[0];

        if (n_phases>1)
          {
            sol_saturation_3p     = &(*saturation_3p)[0];
            if (is_g)
              {
                sol_gas_oil_ratio = &(*gas_oil_ratio)[0];
              }
          }
      }
    else
      {
#ifdef _DEBUG
        BOSOUT (section::iters, level::debug) << ("istart_linear_search") << bs_end;
#endif
        sol_pressure              = &(*prev_niter_data_.pressure)[0];
        sol_main_variable         = &prev_niter_data_.main_var[0];

        if (n_phases>1)
          {
            sol_saturation_3p     = &(*prev_niter_data_.saturation_3p)[0];
            if (is_g)
              {
                sol_gas_oil_ratio = &(*prev_niter_data_.gas_oil_ratio)[0];
              }
          }
      }

    index_t n                     = msh->get_n_active_elements ();
    const item_t *x_sol           = &(*x_sol_)[0];
    const item_t *sec_sol         = &(*sec_sol_)[0];

    index_t n_left = 0;
    index_t n_right = n;

    // calc only internal cells
    for (index_t i = n_left; i < n_right; ++i)
      {
        item_t dsw = 0;
        item_t dpo = 0;
        item_t dsg = 0;
        item_t dro = 0;
        item_t dso = 0;

        index_t jj = (i - n_left);
        if (n_phases == 3)
          {
            dsw = sec_sol[jj];
            dpo = x_sol[jj * n_phases + p3_po];
            dsg = 0;
            dro = 0;
            dso = 0;
            if (sol_main_variable[i] == FI_SG_VAR)
              {
                dsg = x_sol[jj * n_phases + p3_sg];
                dso = x_sol[jj * n_phases + p3_so];
              }
            else if (sol_main_variable[i] == FI_RO_VAR)
              {
                dro = x_sol[jj * n_phases + p3_sg];
                dso = x_sol[jj * n_phases + p3_so];
              }
            else if (sol_main_variable[i] == FI_MOMG_VAR)
              {
                calc_approx_so_sg_ro (x_sol[jj * n_phases + p3_so],
                                      x_sol[jj * n_phases + p3_sg],
                                      data[i].porosity, data[i].invers_fvf[d_o],
                                      data[i].invers_fvf[d_g], sol_gas_oil_ratio[i],
                                      dso, dsg,
                                      dro, main_variable[i]);
              }

            (*saturation_3p)[i * n_phases + d_w] = sol_saturation_3p[i * n_phases + d_w] + dsw * mult;
            (*saturation_3p)[i * n_phases + d_o] = sol_saturation_3p[i * n_phases + d_o] + dso * mult;
            (*saturation_3p)[i * n_phases + d_g] = sol_saturation_3p[i * n_phases + d_g] + dsg * mult;

            if ((*saturation_3p)[i * n_phases + d_w] < 0)
              (*saturation_3p)[i * n_phases + d_w] = 0;
            if ((*saturation_3p)[i * n_phases + d_w] > 1)
              (*saturation_3p)[i * n_phases + d_w] = 1;
            if ((*saturation_3p)[i * n_phases + d_o] < 0)
              (*saturation_3p)[i * n_phases + d_o] = 0;
            if ((*saturation_3p)[i * n_phases + d_o] > 1)
              (*saturation_3p)[i * n_phases + d_o] = 1;

            //int fix_soil_bug = this->ts_params->get_bool(fi_params::FIX_SOIL_BUG);
            //if (fix_soil_bug)
            //  {
            //    if (fabs (saturation_3p[i * n_phases + d_o]) < EPS_DIFF)
            //      {
            //        saturation_3p[i * n_phases + d_o] = 1.0 - 0.001;
            //      }
            //  }

            if ((*saturation_3p)[i * n_phases + d_g] > 1)
              (*saturation_3p)[i * n_phases + d_g] = 1;

            (*gas_oil_ratio)[i] = sol_gas_oil_ratio[i] + dro * mult;
            if ((*gas_oil_ratio)[i] < 0)
              (*gas_oil_ratio)[i] = 0;

#ifdef _DEBUG
            if (main_variable[i] != sol_main_variable[i])
              {
                BOSOUT (section::iters, level::debug) << boost::format ("restore main_var {%d} [%d] -> [%d] (ro == %d)") % i % main_variable[i] % sol_main_variable[i] % FI_RO_VAR << bs_end;
              }
            else
              {
                //BOSOUT (section::iters, level::debug) << boost::format ("keep main_var {%d} [%d] -> [%d] (ro == %d)") % i % main_variable[i] % sol_main_variable[i] % FI_RO_VAR << bs_end;
              }
#endif
            main_variable[i] = sol_main_variable[i];
          }
        else if (is_w && is_o)
          {
            dsw = sec_sol[jj];
            dso = x_sol[jj * n_phases + p2ow_so];
            dpo = x_sol[jj * n_phases + p2ow_po];

            (*saturation_3p)[i * n_phases + d_w] = sol_saturation_3p[i * n_phases + d_w]
                                                + dsw * mult;
            (*saturation_3p)[i * n_phases + d_o] = sol_saturation_3p[i * n_phases + d_o]
                                                + dso * mult;

            if ((*saturation_3p)[i * n_phases + d_w] < 0)
              (*saturation_3p)[i * n_phases + d_w] = 0;
            if ((*saturation_3p)[i * n_phases + d_w] > 1)
              (*saturation_3p)[i * n_phases + d_w] = 1;
            if ((*saturation_3p)[i * n_phases + d_o] < 0)
              (*saturation_3p)[i * n_phases + d_o] = 0;
            if ((*saturation_3p)[i * n_phases + d_o] > 1)
              (*saturation_3p)[i * n_phases + d_o] = 1;
          }
        else if (is_g && is_o)
          {
            dpo = x_sol[jj * n_phases + p2og_po];
            dsg = 0;
            dro = 0;
            dso = 0;
            if (sol_main_variable[i] == FI_SG_VAR)
              {
                dsg = x_sol[jj * n_phases + p2og_sg];
                dso = sec_sol[jj];
              }
            else if (sol_main_variable[i] == FI_RO_VAR)
              {
                dro = x_sol[jj * n_phases + p2og_sg];
                dso = sec_sol[jj];
              }
            else if (sol_main_variable[i] == FI_MOMG_VAR)
              {
                calc_approx_so_sg_ro (sec_sol[jj],
                                      x_sol[jj * n_phases + p2og_sg],
                                      data[i].porosity, data[i].invers_fvf[d_o],
                                      data[i].invers_fvf[d_g], sol_gas_oil_ratio [i],
                                      dso, dsg,
                                      dro, main_variable[i]);
              }
            (*saturation_3p)[i * n_phases + d_o] = sol_saturation_3p[i * n_phases + d_o]
                                                + dso * mult;
            (*saturation_3p)[i * n_phases + d_g] = sol_saturation_3p[i * n_phases + d_g]
                                                + dsg * mult;
            if ((*saturation_3p)[i * n_phases + d_o] < 0)
              (*saturation_3p)[i * n_phases + d_o] = 0;
            if ((*saturation_3p)[i * n_phases + d_o] > 1)
              (*saturation_3p)[i * n_phases + d_o] = 1;
            if ((*saturation_3p)[i * n_phases + d_g] > 1)
              (*saturation_3p)[i * n_phases + d_g] = 1;

            (*gas_oil_ratio)[i] = sol_gas_oil_ratio[i] + dro * mult;
            if ((*gas_oil_ratio)[i] < 0)
              (*gas_oil_ratio)[i] = 0;
          }
        else
          {
            dpo = x_sol[jj];
          }

#ifdef _DEBUG
        //BOSOUT (section::iters, level::debug) << "pressure_mult: " << mult << bs_end;
        //BOSOUT (section::iters, level::debug) << "pressure_dpo_mult: " << dpo * mult << bs_end;
#endif

        item_t dpo_mult = dpo * mult;
        if (dpo_mult > max_p)
          dpo_mult = max_p;
        else if (dpo_mult < -max_p)
          dpo_mult = -max_p;

        (*pressure)[i] = sol_pressure[i] + dpo_mult;
        if ((*pressure)[i] > maximal_pressure)
          (*pressure)[i] = maximal_pressure - 5;
        if ((*pressure)[i] < minimal_pressure)
          (*pressure)[i] = minimal_pressure + 5;
      }

    //static int iter_counter = 0;
    //iter_counter++;
    //tools::save_seq_vector (tools::string_formater ("sat_3p.s.bs.%d.txt", iter_counter).str).save (saturation_3p);
    return 0;
  }

  int
  calc_model::calc_approx_so_sg_ro (const item_t mo_in, const item_t mg_in, const item_t poro,
      const item_t ifvf_o, const item_t ifvf_g, const item_t max_ro,
      // results
      item_t &so, item_t &sg, item_t &ro,
      main_var_type &m_var)
  {
    double d;
    const double delta_mo = 1.0e-5;
    const double delta_mg = 1.0e-5;
    double mo, mg;

    if (poro > EPS_DIV)
      {
        mo = mo_in / poro;
        mg = mg_in / poro;
      }
    else
      {
        mo = 0;
        mg = 0;
      }

    if (mo > delta_mo && mg > delta_mg)
      {
        // step 1: calculate So
        if (ifvf_o > EPS_DIV)
          so = mo / ifvf_o;
        else
          so = (item_t) 0.1;

        // step 2: calculate Sg, Ro
        d = mg / mo;
        if (d > max_ro)
          {
            // Sg > 0
            if (ifvf_g > EPS_DIV)
              sg = (mg - mo * max_ro) / ifvf_g;
            else
              sg = (item_t) 0.1;

            m_var = FI_SG_VAR;
            ro = max_ro;
          }
        else
          {
            m_var = FI_RO_VAR;
            sg = 0;
            ro = d;
          }
      }
    else if (mo > delta_mo)
      {
        // step 1: calculate So
        if (ifvf_o > EPS_DIV)
          so = mo / ifvf_o;
        else
          so = (item_t) 0.1;

        m_var = FI_RO_VAR;
        sg = 0;
        ro = 0;
      }
    else if (mg > delta_mg)
      {
        // step 1: calculate So
        if (ifvf_g > EPS_DIV)
          sg = mg / ifvf_g;
        else
          sg = (item_t) 0.1;

        m_var = FI_SG_VAR;
        so = 0;
        ro = max_ro;
      }
    return 0;
  }

  bool
  calc_model::is_water () const
    {
      return FI_CHK_WATER (phases);
    }
  bool
  calc_model::is_gas () const
    {
      return FI_CHK_GAS (phases);
    }
  bool
  calc_model::is_oil () const
    {
      return FI_CHK_OIL (phases);
    }
  calc_model::index_t
  calc_model::water_shift () const
    {
      return phase_d[FI_PHASE_WATER];
    }
  calc_model::index_t
  calc_model::gas_shift () const
    {
      return phase_d[FI_PHASE_GAS];
    }
  calc_model::index_t
  calc_model::oil_shift () const
    {
      return phase_d[FI_PHASE_OIL];
    }

  BLUE_SKY_TYPE_STD_CREATE (calc_model);
  BLUE_SKY_TYPE_STD_COPY (calc_model);

  BLUE_SKY_TYPE_IMPL (calc_model, bs_node, "calc_model", "Calc_model (double and long)", "Calc_model (double and long)");
}
