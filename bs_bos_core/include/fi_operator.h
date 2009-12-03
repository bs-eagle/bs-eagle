/**
 *       \file  fi_operator.h
 *      \brief  fi_operator (calculates full implicit model)
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  12.01.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Should be moved to src/ directory.
 * */
#ifndef BS_FI_OPERATOR_H_
#define BS_FI_OPERATOR_H_

#include BS_FORCE_PLUGIN_IMPORT ()
#include "switch_main_vars.h"
#include "pvt_water.h"
#include "pvt_dead_oil.h"
#include "pvt_gas.h"
#include "scal_region_info.h"
#include "scal_region.h"
#include "scal_2p_data_holder.h"
#include "scal_3p.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include "calc_model.h"
#include "norm_calc.h"

namespace blue_sky {

  /**
   * \class fi_operator_impl
   * \brief Implementation of fi_operator, parametrized with 
   *        is_w, is_g, is_o (for water, gas, oil phases)
   * \todo  Describe data members
   * */
  template <typename strategy_type, bool is_w, bool is_g, bool is_o>
  struct fi_operator_impl
  {
    typedef strategy_type                                   strategy_t;
    typedef typename strategy_t::item_t                     item_t;
    typedef typename strategy_t::rhs_item_t                 rhs_item_t;
    typedef typename strategy_t::index_t                    index_t;
    typedef typename strategy_t::item_array_t               item_array_t;
    typedef typename strategy_t::rhs_item_array_t           rhs_item_array_t;
    typedef typename strategy_t::index_array_t              index_array_t;
    typedef typename strategy_t::csr_matrix_t               bcsr_matrix_t;

    typedef calc_model <strategy_t>                         calc_model_t;
    typedef typename calc_model_t::data_t                   data_t;
    typedef typename calc_model_t::data_array_t             data_array_t;
    typedef typename calc_model_t::main_var_array_t         main_var_array_t;
    typedef norms_storage <strategy_t>                      norms_storage_t;
    typedef reservoir <strategy_t>                          reservoir_t;
    typedef rs_mesh_iface <strategy_t>                      mesh_iface_t;
    typedef jacobian <strategy_t>                           jacobian_t;
    typedef jacobian_matrix <strategy_t>                    jmatrix_t;

    typedef typename calc_model_t::sp_this_t                sp_calc_model_t;
    typedef typename calc_model_t::sp_reservoir_t           sp_reservoir_t;
    typedef typename calc_model_t::sp_jacobian_t            sp_jacobian_t;
    typedef typename calc_model_t::sp_jacobian_matrix_t     sp_jmatrix_t;
    typedef typename calc_model_t::sp_mesh_iface_t          sp_mesh_iface_t;
    typedef typename calc_model_t::sp_rock_grid             sp_rock_grid_prop_t;
    typedef smart_ptr <bcsr_matrix_t, true>                 sp_bcsr_matrix_t;

    typedef typename calc_model_t::sp_pvt_dead_oil_array_t  sp_pvt_dead_oil_array_t;
    typedef typename calc_model_t::sp_pvt_gas_array_t       sp_pvt_gas_array_t;
    typedef typename calc_model_t::sp_pvt_water_array_t     sp_pvt_water_array_t;

    enum {
      n_phases = is_w + is_g + is_o,
      b_sqr = n_phases * n_phases,
      is_1p = n_phases == 1,
    };

  public:
    /**
     * \brief  fi_operator_impl ctor
     * \param  calc_model
     * \param  reservoir
     * \param  mesh
     * \param  jacobian
     * \param  jmatrix
     * */
    fi_operator_impl (sp_calc_model_t &calc_model, sp_reservoir_t &reservoir, const sp_mesh_iface_t &mesh, sp_jacobian_t &jacobian, sp_jmatrix_t &jmatrix)
    : calc_model_ (calc_model),
    rock_grid_prop_ (calc_model_->rock_grid_prop),
    reservoir_ (reservoir),
    jacobian_ (jacobian),
    jmatrix_ (jmatrix),
    mesh_ (mesh),
    trns_matrix_ (jmatrix_->trns_matrix),
    trns_values_ (trns_matrix_->get_values ()),
    trns_rows_ptr_ (trns_matrix_->get_rows_ptr ()),
    trns_cols_ptr_ (trns_matrix_->get_cols_ind ()),
    reg_matrix_ (jmatrix_->get_regular_matrix ()),
    reg_values_ (reg_matrix_->get_values ()),
    reg_rows_ptr_ (reg_matrix_->get_rows_ptr ()),
    reg_cols_ptr_ (reg_matrix_->get_cols_ind ()),
    m_array_ (jmatrix_->m_array),
    p_array_ (jmatrix_->p_array),
    rhs_ (jmatrix_->get_rhs ()),
    sol_ (jmatrix_->get_solution ()),
    flux_rhs_ (jmatrix_->get_rhs_flux ()),
    sp_diag_ (jmatrix_->get_sp_diagonal ()),
    s_rhs_ (jmatrix_->get_sec_rhs ()),
    depths_ (mesh_->get_depths ()),
    n_cells_ (mesh->get_n_active_elements()),
    n_connections_ (mesh_->get_n_connections ()),
    n_sec_vars (calc_model_->n_sec_vars),
    d_w (calc_model_->phase_d[FI_PHASE_WATER]),
    d_g (calc_model_->phase_d[FI_PHASE_GAS]),
    d_o (calc_model_->phase_d[FI_PHASE_OIL]),
    ds_w (calc_model_->sat_d[FI_PHASE_WATER]),
    ds_g (calc_model_->sat_d[FI_PHASE_GAS]),
    d_gg (d_g * n_phases + d_g),
    d_gw (d_g * n_phases + d_w),
    d_go (d_g * n_phases + d_o),
    d_wg (d_w * n_phases + d_g),
    d_ww (d_w * n_phases + d_w),
    d_wo (d_w * n_phases + d_o),
    d_og (d_o * n_phases + d_g),
    d_ow (d_o * n_phases + d_w),
    d_oo (d_o * n_phases + d_o),
    data_ (calc_model_->data),
    saturation_3p_ (calc_model_->saturation_3p),
    pressure_ (calc_model_->pressure),
    gas_oil_ratio_ (calc_model_->gas_oil_ratio),
    main_vars_ (calc_model_->main_variable),
    volume_ (rock_grid_prop_->volume),
    poro_array_ (rock_grid_prop_->porosity_p_ref),
    rock_grid_comp_const_ (rock_grid_prop_->comp_const),
    rock_grid_comp_ref_pressure_ (rock_grid_prop_->comp_ref_pressure),
    sat_regions_ (calc_model_->sat_regions),
    pvt_regions_ (calc_model_->pvt_regions),
    pvt_oil_array (calc_model_->pvt_oil_array),
    pvt_water_array (calc_model_->pvt_water_array),
    pvt_gas_array (calc_model_->pvt_gas_array),
    min_p_ ((item_t) calc_model_->ts_params->get_float (fi_params::PVT_PRESSURE_RANGE_MIN)),
    max_p_ ((item_t) calc_model_->ts_params->get_float (fi_params::PVT_PRESSURE_RANGE_MAX)),
    drsdt_ ((item_t) calc_model_->ts_params->get_float (fi_params::DRSDT)),
    rhs_residual_ ((item_t) calc_model_->ts_params->get_float (fi_params::NEWTON_RESIDUAL)),
    mb_error_ ((item_t) calc_model_->ts_params->get_float (fi_params::MASS_BALANS_ERROR)),
    s_rhs_norm_ ((item_t) calc_model_->ts_params->get_float (fi_params::S_RHS_NORM)),
    norm_ (calc_model_->norm),
    cfl_ (jmatrix_->get_cfl_vector ())
    {
    }

    /**
     * \brief       Main fi_operator function
     * \param[out]  dt
     * \param       istart
     * \param       istart_well_contr
     * \param       __formal
     * \param       update_rhs_after_gauss_elimination
     * \param       save_debug_files
     * \return      
     * \todo        Describe return type
     * */
    fi_operator_return_type
    fi_operator (double &dt, index_t istart, index_t istart_well_contr, index_t &, bool update_rhs_after_gauss_elimination, bool save_debug_files)
    {
      save_debug_files;
      item_t mult = 1.0;
      bool tuning = calc_model_->ts_params->get_bool (fi_params::NEWTON_TUNING);

      static norms_storage_t base_norm;

      if (istart)
        {
          reservoir_->init_jacobian (jmatrix_, n_cells_);
        }

      index_t n_approx = 5;
      for (index_t i = 0; i < n_approx; ++i)
        {
          mult = 1.0;

          // prepare well jacobian part
          jmatrix_->init (n_cells_, n_phases, 3, 0, n_sec_vars);

          if (is_o && is_g)
            {
              fi_operator_switch_main_vars (dt);
            }

          if (n_phases > 1)
            {
              // calculate saturation properties
              calc_model_->scal_prop->process (saturation_3p_,
                sat_regions_,
                rock_grid_prop_->permeability,
                poro_array_,
                data_);
            }

          // calculate properties for cells
          // fill accumulative part of Jacobian (dM/dt)
          fi_operator_cells ((istart == 1) ? 1 : 0, dt);
          fi_operator_fill ();

          //if (!rsv_status->fip_data_is_calculated)
          //{
          //  rsv_status->fi_calculate_fip_data (this, /*precsn*/ 1e-8);
          //  rsv_status->fip_data_is_calculated = 1;
          //  copy_fip_data_to_storage (0.0);
          //}

          // calculate Flux part of jacobian (T * M* (Pi - Pj - H))
          // calculate boundary Fluxes

          //!TODO:
          if (block_connections_mpfa (dt))// || fi_operator_boundary_block_connections (*dt))
            throw bs_exception("fi_operator_block_connections_mpfa", "return -1");

          if (update_rhs_after_gauss_elimination)
            {
              // calculate well controls and rates
              reservoir_->calc_wells (((istart_well_contr == 1) ? 1 : 0), dt, calc_model_, mesh_, jmatrix_);

              // fill rhs from wells
              reservoir_->fill_rhs_wells (dt, calc_model_, flux_rhs_, update_rhs_after_gauss_elimination);
            }

          norm_calc ();

          if (update_rhs_after_gauss_elimination)
            {
              // fill jacobian (irregular matrix) by wells
              reservoir_->end_jacobian (dt, calc_model_, jacobian_);
            }

#ifdef _DEBUG
          //if (save_debug_files)
            {
              debug_save_data (dt);
            }
#endif

          if (check_norm (istart))
            return FI_OPERATOR_RETURN_FAIL;

          if (check_solution_mult_cell (istart, i, n_approx, base_norm))
            continue;

          do_dt_reduce (dt, istart);
          do_dt_tuning (dt, tuning);

          break;
        }

      base_norm = norm_;

      jmatrix_->clear_solution ();        // clear solution vector
      jmatrix_->summ_rhs ();              // summarize rhs

      return FI_OPERATOR_RETURN_OK;
    }

    /**
     * \brief  Inits fi_operator process, called from main_loop_calc
     * \param  istart
     * \param  dt
     * \return 
     * */
    void
    fi_operator_init (index_t istart, double dt)
    {
      if (is_o && is_g)
        {
          fi_operator_switch_main_vars (dt);
        }

      if (n_phases > 1)
        {
          // calculate saturation properties
          calc_model_->scal_prop->process (saturation_3p_,
            sat_regions_,
            rock_grid_prop_->permeability,
            poro_array_,
            data_);
        }

      fi_operator_cells (istart, dt);
    }

    /**
     * \brief  Calculates physical parameters for all cells
     * \param  
     * \return 
     * */
    void
    fi_operator_cells (index_t istart, const item_t dt)
    {
      // local declaration
      index_t i;
      rhs_item_t *jac_block     = 0;
      rhs_item_t *rhs_block     = 0;
      index_t b_sqr         = n_phases * n_phases;

      rhs_item_array_t &main_diag_acc = jmatrix_->get_regular_acc_diag();//get_main_diagonal_accumulative ();
      rhs_item_array_t &ss_diag       = jmatrix_->get_ss_diagonal ();
      rhs_item_array_t &sp_diag       = jmatrix_->get_sp_diagonal ();
      rhs_item_array_t &s_rhs         = jmatrix_->get_sec_rhs ();
      rhs_item_array_t &rhs           = jmatrix_->get_rhs ();

      BS_ERROR (!main_diag_acc.empty (), "fi_operator_cells");
      BS_ERROR (!ss_diag.empty (), "fi_operator_cells");
      BS_ERROR (!sp_diag.empty (), "fi_operator_cells");
      BS_ERROR (!s_rhs.empty (), "fi_operator_cells");
      BS_ERROR (!rhs.empty (), "fi_operator_cells");

      index_t switch_to_sg_count = 0;
      index_t switch_to_ro_count = 0;
      index_t switch_to_momg_count = 0;
      item_t total_volume = 0.0;
      rhs_item_t *ss_block = 0;
      rhs_item_t *s_rhs_block = 0;

#ifdef _MPI
      int n_left = 0;///mpi_decomp->get_recv_l ();
      int n_own = 0;///mpi_decomp->get_n_local_own () + n_left;
      ///rhs -= n_left * n_phases;
      double mpi_invers_fvf_average[3];
      int n_procs;
#endif //_MPI

#ifdef FI_OPERATOR_CELLS_PARALLEL
      double invers_fvf_average_w = 0.;
      double invers_fvf_average_g = 0.;
      double invers_fvf_average_o = 0.;
#endif //FI_OPERATOR_CELLS_PARALLEL


      // set flag to 0 at newton iteration
      calc_model_->lsearch_force_newton_step = 0;

      // set equal to 0
      assign (calc_model_->invers_fvf_average, 0);

      // set flag to 0 at newton iteration
#ifdef FI_OPERATOR_CELLS_PARALLEL
#pragma omp parallel for private  (jac_block, rhs_block) \
    reduction (+:invers_fvf_average_w, invers_fvf_average_g, invers_fvf_average_o,total_volume)
#endif //FI_OPERATOR_CELLS_PARALLEL
      // loop through all cells
      for (i = 0; i < n_cells_; ++i)
        {
          jac_block   = &main_diag_acc[b_sqr * i];
          ss_block    = &ss_diag[n_sec_vars * n_sec_vars * i];
          s_rhs_block = &s_rhs[n_sec_vars * i];
          rhs_block   = &rhs[i * n_phases];
#ifdef _MPI
          BS_ASSERT (false && "MPI: NOT IMPL YET");
          if ((i < n_left) || (i >= n_own))
            {
              ///jac_block = just_double;
              ///rhs_block = just_double;
            }
          else
            {
              jac_block -= n_left * b_sqr;
            }
#endif //_MPI

          // calculate properties for cell i
#ifdef FI_OPERATOR_CELLS_PARALLEL
          fi_operator_cell (istart, dt, i, jac_block, rhs_block,
                            ss_block, s_rhs_block,
                            switch_to_sg_count, switch_to_ro_count, switch_to_momg_count,
                            invers_fvf_average_w, invers_fvf_average_g, invers_fvf_average_o, total_volume);
#else //FI_OPERATOR_CELLS_PARALLEL
          fi_operator_cell (istart, dt, i, jac_block, rhs_block,
                            ss_block, s_rhs_block,
                            switch_to_sg_count, switch_to_ro_count, switch_to_momg_count,
                            total_volume);
#endif //FI_OPERATOR_CELLS_PARALLEL
        }

#ifdef _MPI
      BS_ASSERT (false && "MPI: NOT IMPL YET");
      double mpi_total_volume;
      MPI_Allreduce (&total_volume, &mpi_total_volume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      total_volume = mpi_total_volume;

      MPI_Allreduce (&invers_fvf_average, &mpi_invers_fvf_average, FI_PHASE_TOT, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      if (is_w)
        invers_fvf_average[d_w] = mpi_invers_fvf_average[d_w];
      if (is_g)
        invers_fvf_average[d_g] = mpi_invers_fvf_average[d_g];
      if (is_o)
        invers_fvf_average[d_o] = mpi_invers_fvf_average[d_o];

      ///n_elements = mpi_decomp->get_n_elements ();
#endif //_MPI

#ifdef FI_OPERATOR_CELLS_PARALLEL
      if (is_w)
        invers_fvf_average[d_w] = invers_fvf_average_w / (double) n_elements;
      if (is_g)
        invers_fvf_average[d_g] = invers_fvf_average_g / (double) n_elements;
      if (is_o)
        invers_fvf_average[d_o] = invers_fvf_average_o / (double) n_elements;
#else //FI_OPERATOR_CELLS_PARALLEL

      if (is_w)
        calc_model_->invers_fvf_average[d_w] /= (item_t) n_cells_;
      if (is_g)
        calc_model_->invers_fvf_average[d_g] /= (item_t) n_cells_;
      if (is_o)
        calc_model_->invers_fvf_average[d_o] /= (item_t) n_cells_;
#endif //FI_OPERATOR_CELLS_PARALLEL

      calc_model_->ave_volume = total_volume / (item_t) n_cells_;
    }

#define POROSITY                data_i.porosity
#define P_DERIV_POROSITY        data_i.p_deriv_porosity
#define INVERS_FVF_O            data_i.invers_fvf[d_o]
#define INVERS_FVF_G            data_i.invers_fvf[d_g]
#define INVERS_FVF_W            data_i.invers_fvf[d_w]
#define P_DERIV_INVERS_FVF_O    data_i.p_deriv_invers_fvf[d_o]
#define P_DERIV_INVERS_FVF_G    data_i.p_deriv_invers_fvf[d_g]
#define P_DERIV_INVERS_FVF_W    data_i.p_deriv_invers_fvf[d_w]
#define GAS_OIL_RATIO           gas_oil_ratio_[i]
#define P_DERIV_GAS_OIL_RATIO   data_i.p_deriv_gas_oil_ratio
#define VOLUME                  volume_[i]
#define S_DERIV_CAP_PRESSURE_G  data_i.s_deriv_cap_pressure[ds_g]
#define S_DERIV_CAP_PRESSURE_W  data_i.s_deriv_cap_pressure[ds_w]
#define GOR_DERIV_INVERS_FVF    data_i.gor_deriv_invers_fvf
#define PREV_FLUID_VOLUME_O     data_i.prev_fluid_volume[d_o]
#define PREV_FLUID_VOLUME_G     data_i.prev_fluid_volume[d_g]
#define PREV_FLUID_VOLUME_W     data_i.prev_fluid_volume[d_w]

    /**
     * \class local_data
     * \brief Local data that used in calculation of each cell
     * */
    struct local_data
    {
      item_t                    sat_w;      //!< Saturation value for water phase
      item_t                    sat_g;      //!< Saturation value for gas phase
      item_t                    sat_o;      //!< Saturation value for oil phase
      boost::array <item_t, 3>  ps_block;   //!< PS block
    };

    /**
     * \brief  Fills Jacobain part for cell i
     * \param  i Index of mesh cell
     * \param  data_i 
     * \param  jac_block Jacobian value for cell
     * \param  local_data_
     * */
    BS_FORCE_INLINE void
    fill_jacobian (index_t i, const data_t &data_i, rhs_item_t *jac_block, local_data &local_data_)
    {
      // water equation
      if (is_w)
        {
          item_t p_der  = VOLUME * (P_DERIV_POROSITY * INVERS_FVF_W + POROSITY * P_DERIV_INVERS_FVF_W);
          item_t sw_der = VOLUME * POROSITY * (INVERS_FVF_W + local_data_.sat_w * P_DERIV_INVERS_FVF_W * S_DERIV_CAP_PRESSURE_W);

          if (n_phases == 3)
            {
              jac_block[p3_wat_po] = p_der * local_data_.sat_w;
              jac_block[p3_wat_sg] = 0;
              jac_block[p3_wat_so] = 0;
              local_data_.ps_block[p3_wat] = sw_der;
            }
          else if (n_phases == 2)
            {
              jac_block[p2ow_wat_po] = p_der * local_data_.sat_w;
              jac_block[p2ow_wat_so] = 0;
              local_data_.ps_block[p2ow_wat] = sw_der;
            }
          else
            {
              jac_block[0] = p_der;
            }
        }
      // gas equation
      if (is_g)
        {
          if (n_phases > 1)
            {
              item_t drg_p = 0, drg_sg = 0, drg_so = 0;
              if (FI_CHK_SG (main_vars_, i))
                {
                  //drg_p  = VOLUME * (sat_g * (P_DERIV_POROSITY * INVERS_FVF_G + POROSITY * P_DERIV_INVERS_FVF_G) +
                  //                   sat_o * (GAS_OIL_RATIO * (P_DERIV_POROSITY * INVERS_FVF_O + POROSITY * P_DERIV_INVERS_FVF_O)
                  //                            + POROSITY * INVERS_FVF_O * P_DERIV_GAS_OIL_RATIO));
                  drg_p  = VOLUME * (local_data_.sat_g * (P_DERIV_POROSITY * INVERS_FVF_G + POROSITY * P_DERIV_INVERS_FVF_G) +
                                     local_data_.sat_o * (P_DERIV_POROSITY * INVERS_FVF_O * GAS_OIL_RATIO + POROSITY * (P_DERIV_GAS_OIL_RATIO * INVERS_FVF_O + GAS_OIL_RATIO * P_DERIV_INVERS_FVF_O)));
                  drg_sg = VOLUME * POROSITY * (INVERS_FVF_G + local_data_.sat_g * P_DERIV_INVERS_FVF_G * S_DERIV_CAP_PRESSURE_G - GAS_OIL_RATIO * INVERS_FVF_O);
                }
              else if (FI_CHK_RO (main_vars_, i))
                {
                  drg_p  = VOLUME * local_data_.sat_o * GAS_OIL_RATIO * (P_DERIV_POROSITY * INVERS_FVF_O + P_DERIV_INVERS_FVF_O * POROSITY);
                  drg_sg = VOLUME * POROSITY * local_data_.sat_o * (INVERS_FVF_O + GAS_OIL_RATIO * GOR_DERIV_INVERS_FVF);
                }
              drg_so = VOLUME * POROSITY * GAS_OIL_RATIO * INVERS_FVF_O;
              if (FI_CHK_MOMG (main_vars_, i))
                {
                  drg_p   = 0;
                  drg_sg  = VOLUME;
                  drg_so  = 0;
                }
              if (n_phases == 3)
                {
                  jac_block[p3_gas_po] = drg_p;
                  jac_block[p3_gas_sg] = drg_sg;
                  jac_block[p3_gas_so] = drg_so;
                  local_data_.ps_block[p3_gas] = 0;
                }
              else
                {
                  jac_block[p2og_gas_po] = drg_p;
                  jac_block[p2og_gas_sg] = drg_sg;
                  local_data_.ps_block[p2og_gas] = drg_so;
                }
            }
          else
            {
              jac_block[0] = VOLUME * (P_DERIV_POROSITY * INVERS_FVF_G + P_DERIV_INVERS_FVF_G * POROSITY);
            }

        }
      // oil equation
      if (is_o)
        {
          if (n_phases == 3)
            {
              jac_block[p3_oil_po] = VOLUME * local_data_.sat_o * (P_DERIV_POROSITY * INVERS_FVF_O + P_DERIV_INVERS_FVF_O * POROSITY);
              jac_block[p3_oil_so] = /*-*/VOLUME * POROSITY * INVERS_FVF_O;

              if (FI_CHK_SG (main_vars_, i))
                {
                  //jac_block[p3_oil_sg] = 0.0;
                  jac_block[p3_oil_sg] = -VOLUME * POROSITY * INVERS_FVF_O;
                }
              else if (FI_CHK_RO (main_vars_, i))
                {
                  jac_block[p3_oil_sg] = VOLUME * POROSITY * local_data_.sat_o * GOR_DERIV_INVERS_FVF;
                }
              local_data_.ps_block[p3_oil] = 0.0;

              if (FI_CHK_MOMG (main_vars_, i))
                {
                  jac_block[p3_oil_po] = 0;
                  jac_block[p3_oil_so] = VOLUME;
                  jac_block[p3_oil_sg] = 0.0;
                }
            }
          else if (n_phases == 2 && is_g)
            {
              jac_block[p2og_oil_po] = VOLUME * local_data_.sat_o * (P_DERIV_POROSITY * INVERS_FVF_O + P_DERIV_INVERS_FVF_O * POROSITY);

              if (FI_CHK_SG (main_vars_, i))
                jac_block[p2og_oil_sg] = 0.0;
              else if (FI_CHK_RO (main_vars_, i))
                {
                  jac_block[p2og_oil_sg] = VOLUME * POROSITY * local_data_.sat_o * GOR_DERIV_INVERS_FVF;
                }
              local_data_.ps_block[p2og_oil] = /*-*/VOLUME * POROSITY * INVERS_FVF_O;
              if (FI_CHK_MOMG (main_vars_, i))
                {
                  jac_block[p2og_oil_po] = 0;
                  local_data_.ps_block[p2og_oil] = VOLUME;
                  jac_block[p2og_oil_sg] = 0.0;
                }
            }
          else if (n_phases == 2 && is_w)
            {
              jac_block[p2ow_oil_po] = VOLUME * local_data_.sat_o * (P_DERIV_POROSITY * INVERS_FVF_O + P_DERIV_INVERS_FVF_O * POROSITY);
              jac_block[p2ow_oil_so] = /*-*/VOLUME * POROSITY * INVERS_FVF_O;
              local_data_.ps_block[p2ow_oil] = 0.0;
            }
          else
            {
              jac_block[0] = VOLUME * (P_DERIV_POROSITY * INVERS_FVF_O + P_DERIV_INVERS_FVF_O * POROSITY);
            }
        }
    }

    /**
     * \brief Switches main vars for gas phase (Undersaturated conditon or free gas exist)
     * \param dt
     * */
    void
    fi_operator_switch_main_vars (double dt)
    {
      typedef typename calc_model_t::sp_pvt_oil sp_pvt_oil_t;

      sp_pvt_oil_t pvt_oil;
      index_t prev_pvt_reg = -1;

      index_t switch_to_sg_count = 0;
      index_t switch_to_ro_count = 0;
      index_t switch_to_momg_count = 0;

      for (index_t i = 0, cnt = n_cells_; i < cnt; ++i)
        {
          index_t pvt_reg = pvt_regions_[i];
          if (pvt_reg != prev_pvt_reg)
            {
              prev_pvt_reg = pvt_reg;
              pvt_oil = pvt_oil_array [pvt_reg];
            }

          switch_main_vars <strategy_t>::do_switch (
            is_w, is_g, is_o,
            d_o, d_g, d_w,
            pvt_oil_array [pvt_reg],
            pressure_[i],
            main_vars_[i],
            &saturation_3p_[i * n_phases],
            gas_oil_ratio_[i],
            calc_model_->lsearch_force_newton_step,
            drsdt_,
            dt,
            calc_model_->old_data_.gas_oil_ratio[i],
            switch_to_sg_count,
            switch_to_ro_count,
            switch_to_momg_count,
            i);
        }

      BOSOUT (section::iters, level::medium)
        << "Switch to Sg: " << switch_to_sg_count
        << "\tSwitch to Ro: "     << switch_to_ro_count
        << "\tSwitch to MoMg: "   << switch_to_momg_count
        << bs_end;
    }

    /**
     * \brief  Fills accumulative rhs part for cell i
     * \param  i Index of cell
     * \param  data_i
     * \param  rhs_block
     * \param  local_data_
     * */
    BS_FORCE_INLINE void
    fill_acc_rhs (index_t i, const data_t &data_i, rhs_item_t *rhs_block, local_data &local_data_)
    {
      // water
      if (is_w)
        {
          item_t r_wat = 0;
          item_t sat_w__ = n_phases > 1 ? local_data_.sat_w : 1.0;

          r_wat = -VOLUME * (POROSITY * INVERS_FVF_W * sat_w__ - PREV_FLUID_VOLUME_W);
          //if (istart)
          //  r_wat = 0.0;

          if (n_phases == 3)
            rhs_block[p3_wat] = r_wat;
          else if (n_phases == 2)
            rhs_block[p2ow_wat] = r_wat;
          else
            rhs_block[0] = r_wat;
        }
      // gas
      if (is_g)
        {
          item_t r_gas = 0;

          if (n_phases > 1)
            {
              if (FI_CHK_SG (main_vars_, i))
                {
                  r_gas = -VOLUME * (POROSITY * INVERS_FVF_G * local_data_.sat_g + GAS_OIL_RATIO * POROSITY * local_data_.sat_o * INVERS_FVF_O - PREV_FLUID_VOLUME_G);
                }
              else if (FI_CHK_RO (main_vars_, i))
                {
                  r_gas = -VOLUME * (GAS_OIL_RATIO * POROSITY * local_data_.sat_o * INVERS_FVF_O - PREV_FLUID_VOLUME_G);
                }
              else if (FI_CHK_MOMG (main_vars_, i))
                {
                  r_gas = VOLUME * PREV_FLUID_VOLUME_G;
                }
            }
          else
            {
              r_gas = -VOLUME * (POROSITY * INVERS_FVF_G - PREV_FLUID_VOLUME_G);
            }
          //if (istart)
          //  r_gas = 0.0;

          if (n_phases == 3)
            rhs_block[p3_gas] = r_gas;
          else if (n_phases == 2)
            rhs_block[p2og_gas] = r_gas;
          else
            rhs_block[0] = r_gas;
        }
      // oil
      if (is_o)
        {
          item_t r_oil = 0;
          item_t sat_o__ = n_phases > 1 ? local_data_.sat_o : 1.0;

          r_oil = -VOLUME * (POROSITY * INVERS_FVF_O * sat_o__ - PREV_FLUID_VOLUME_O);
          //if (istart)
          //  r_oil = 0.0;

          if (n_phases == 3)
            rhs_block[p3_oil] = r_oil;
          else if (n_phases == 2 && is_w)
            rhs_block[p2ow_oil] = r_oil;
          else if (n_phases == 2 && is_g)
            rhs_block[p2og_oil] = r_oil;
          else
            rhs_block[0] = r_oil;
        }
    }

    /**
     * \brief  Eliminates cell
     * \param  local_data_
     * \param  jac_block
     * \param  rhs_block
     * \param  sp_block
     * \param  s_rhs_block
     * */
    BS_FORCE_INLINE void
    eliminate_cell (const local_data &local_data_, rhs_item_t *jac_block, rhs_item_t *rhs_block, rhs_item_t *sp_block, rhs_item_t *s_rhs_block)
    {
      // Schur complemet 1) App-ApsDssDsp and 2) Bp-ApsDssBs (Dss=Ass^(-1); Dsp=Asp)
      // 1)
      M_MINUS_VV_PROD (n_phases, local_data_.ps_block, sp_block, jac_block);
      // 2)
      V_MINUS_VS_PROD (n_phases, local_data_.ps_block, s_rhs_block, rhs_block);
    }

    /**
     * \brief  Fills Jacobian, accumulative rhs and eliminates all cells
     * */
    void
    fi_operator_fill ()
    {
      index_t b_sqr                   = n_phases * n_phases;
      rhs_item_array_t &main_diag_acc = jmatrix_->get_regular_acc_diag();//get_main_diagonal_accumulative ();
      rhs_item_array_t &ss_diag       = jmatrix_->get_ss_diagonal ();
      rhs_item_array_t &sp_diag       = jmatrix_->get_sp_diagonal ();
      rhs_item_array_t &s_rhs         = jmatrix_->get_sec_rhs ();
      rhs_item_array_t &rhs           = jmatrix_->get_rhs ();

      BS_ASSERT (!main_diag_acc.empty ());
      BS_ASSERT (!ss_diag.empty ());
      BS_ASSERT (!sp_diag.empty ());
      BS_ASSERT (!s_rhs.empty ());
      BS_ASSERT (!rhs.empty ());

      // loop through all cells
      for (index_t i = 0; i < n_cells_; ++i)
        {
          rhs_item_t *jac_block   = &main_diag_acc[b_sqr * i];
          rhs_item_t *ss_block    = &ss_diag[n_sec_vars * n_sec_vars * i];
          rhs_item_t *sp_block    = &sp_diag[n_sec_vars * n_phases * i];
          rhs_item_t *s_rhs_block = &s_rhs[n_sec_vars * i];
          rhs_item_t *rhs_block   = &rhs[i * n_phases];

          index_t i_w = FI_PH_IND (i, d_w, n_phases);
          index_t i_g = FI_PH_IND (i, d_g, n_phases);
          index_t i_o = FI_PH_IND (i, d_o, n_phases);

          local_data local_data_;
          const data_t &data_i = data_[i];

          // fill ss and sp part of blocks
          if (n_phases > 1)
            ss_block[0] = 1;
          if (n_phases == 3)
            {
              s_rhs_block[0] = 1.0 - saturation_3p_[i_w] - saturation_3p_[i_g] - saturation_3p_[i_o];
              if (FI_CHK_SG (main_vars_, i))
                {
                  sp_block[p3_sg] = 1.0;
                  sp_block[p3_so] = 1.0;
                  sp_block[p3_po] = 0;
                }
              else if (FI_CHK_RO (main_vars_, i))
                {
                  sp_block[p3_sg] = 0;
                  sp_block[p3_so] = 1.0;
                  sp_block[p3_po] = 0;
                }
              else if (FI_CHK_MOMG (main_vars_, i))
                {
                  sp_block[p3_sg] = 0;
                  sp_block[p3_so] = 0;
                  sp_block[p3_po] = 0;
                }
            }
          else if (n_phases == 2 && is_w && is_o)
            {
              s_rhs_block[0] = 1.0 - saturation_3p_[i_w] - saturation_3p_[i_o];
              sp_block[p2ow_so] = 1;
              sp_block[p2ow_po] = 0;
            }
          else if (n_phases == 2 && is_g && is_o)
            {
              s_rhs_block[0] = 1.0 - saturation_3p_[i_g] - saturation_3p_[i_o];
              if (FI_CHK_SG (main_vars_, i))
                {
                  sp_block[p2og_sg] = 1;
                  sp_block[p2og_po] = 0;
                }
              else if (FI_CHK_RO (main_vars_, i))
                {
                  sp_block[p2og_sg] = 0;
                  sp_block[p2og_po] = 0;
                }
            }

          assign (local_data_.ps_block, 0);
          local_data_.sat_w = local_data_.sat_g = local_data_.sat_o = 0.0;
          // set up saturation
          if (is_w && n_phases > 1)
            local_data_.sat_w = saturation_3p_[i_w];
          if (is_g && n_phases > 1)
            local_data_.sat_g = saturation_3p_[i_g];
          if (is_o && n_phases > 1)
            local_data_.sat_o = saturation_3p_[i_o];

          fill_jacobian (i, data_i, jac_block, local_data_);
          fill_acc_rhs (i, data_i, rhs_block, local_data_);
          eliminate_cell (local_data_, jac_block, rhs_block, sp_block, s_rhs_block);

//#ifdef _DEBUG
//          for (index_t j = 0; j < b_sqr; ++j)
//            {
//              if (jac_block[j] == 0)
//                {
//                  index_t i_coord = 0, j_coord = 0, k_coord = 0;
//                  //mesh_->get_element_int_to_ijk (i, i_coord, j_coord, k_coord);
//
//                  size_t iiiiiiiiii = 0;
//                }
//            }
//#endif

        }
    }

    /**
     * \brief Calculates physical params for cell
     * \param istart
     * \param dt
     * \param cell_ind
     * \param jac_block
     * \param rhs_block
     * \param ss_block
     * \param s_rhs_block
     * \param switch_to_sg_count
     * \param switch_to_ro_count
     * \param switch_to_momg_count
     * \param total_volume
     * \todo  Remove obsolete params
     * */
#ifdef FI_OPERATOR_CELLS_PARALLEL
    inline void
    fi_operator_cell (index_t istart, const item_t dt,
        index_t cell_ind, item_t *jac_block,
        item_t *rhs_block,
        item_t *ss_block,
        item_t *s_rhs_block,
        int &switch_to_sg_count,
        int &switch_to_ro_count,
        int &switch_to_momg_count,
        item_t &invers_fvf_average_w,
        item_t &invers_fvf_average_g,
        item_t &invers_fvf_average_o,
        item_t &total_volume,
        const sp_mesh_iface_t &mesh)
#else //FI_OPERATOR_CELLS_PARALLEL
    inline void
    fi_operator_cell (index_t /*istart*/, const item_t dt,
      index_t cell_ind, rhs_item_t * /*jac_block*/,
      rhs_item_t * /*rhs_block*/,
        rhs_item_t * /*ss_block*/,
        rhs_item_t * /*s_rhs_block*/,
        int &/*switch_to_sg_count*/,
        int &/*switch_to_ro_count*/,
        int &/*switch_to_momg_count*/,
        item_t &total_volume)
#endif //FI_OPERATOR_CELLS_PARALLEL

    {
      local_data local_data_;
      item_t p_water;
      item_t p_gas;
      item_t p_oil;
      index_t i_temp;

#ifdef _MPI
      int n_left = 0;///mpi_decomp->get_recv_l ();
      int n_own = 0;///mpi_decomp->get_n_local_own () + n_left;
#endif //_MPI

      // base constants
      index_t i = cell_ind;

      item_t gor = 0.0;
      item_t d_gor = 0.0;
      // get saturation region number
      index_t pvt_reg = pvt_regions_[i];

      assign (local_data_.ps_block, 0);
      // check pressure for consistency
      if (pressure_[i] < min_p_ || pressure_[i] > max_p_)
        {
          BOSERR (section::iters, level::error) << boost::format ("Pressure %f in cell [%d] is out of range") %
            pressure_[i] % i << bs_end;
        }

      data_t &data_i = data_[i];

      // calculate pressures for all phases
      p_oil = p_water = p_gas = pressure_[i];
      if (n_phases > 1)
        {
          if (is_w)
            p_water += data_i.cap_pressure[ds_w];
          if (is_g)
            p_gas += data_i.cap_pressure[ds_g];
        }

      // calculate pvt properties
      i_temp = n_phases * i;
      if (is_w)
        {
          pvt_water_array[pvt_reg]->calc(
            p_water,
            &data_i.invers_fvf[d_w],
            &data_i.p_deriv_invers_fvf[d_w],
            &data_i.invers_viscosity[d_w],
            &data_i.p_deriv_invers_viscosity[d_w],
            &data_i.invers_visc_fvf[d_w],
            &data_i.p_deriv_invers_visc_fvf[d_w]);
        }

      if (is_g)
        {
          pvt_gas_array[pvt_reg]->calc(
            p_gas,
            &data_i.invers_fvf[d_g],
            &data_i.p_deriv_invers_fvf[d_g],
            &data_i.invers_viscosity[d_g],
            &data_i.p_deriv_invers_viscosity[d_g],
            &data_i.invers_visc_fvf[d_g],
            &data_i.p_deriv_invers_visc_fvf[d_g]);
        }

      if (is_o && is_g) // gas and oil
        {
          pvt_oil_array[pvt_reg]->calc(
            is_g, main_vars_[i],p_oil, gas_oil_ratio_[i],
            &data_i.invers_fvf[d_o],
            &data_i.p_deriv_invers_fvf[d_o],
            &data_i.gor_deriv_invers_fvf,
            &data_i.invers_viscosity[d_o],
            &data_i.p_deriv_invers_viscosity[d_o],
            &data_i.gor_deriv_invers_viscosity,
            &data_i.invers_visc_fvf[d_o],
            &data_i.p_deriv_invers_visc_fvf[d_o],
            &data_i.gor_deriv_invers_visc_fvf,
            &gor, &d_gor, drsdt_, dt, calc_model_->old_data_.gas_oil_ratio[i]
          );
        }
      else if (is_o) // oil only
        {
          pvt_oil_array[pvt_reg]->calc(
            is_g, -1, p_oil, 0,
            &data_i.invers_fvf[d_o],
            &data_i.p_deriv_invers_fvf[d_o],
            0,
            &data_i.invers_viscosity[d_o],
            &data_i.p_deriv_invers_viscosity[d_o],
            0,
            &data_i.invers_visc_fvf[d_o],
            &data_i.p_deriv_invers_visc_fvf[d_o],
            0,
            0, 0
          );
        }

      // calculate porosity and derivate
      calc_porosity_and_deriv (i,
                               pvt_reg,
                               &data_i.porosity,
                               &data_i.p_deriv_porosity,
                               &data_i.truns_mult,
                               &data_i.p_deriv_truns_mult);

      total_volume += data_i.porosity * volume_[i];

      item_t oil_surface_density    = pvt_oil_array[pvt_reg]->get_surface_density ();
      item_t gas_surface_density    = pvt_gas_array[pvt_reg]->get_surface_density ();
      item_t water_surface_density  = pvt_water_array[pvt_reg]->get_surface_density ();

      // calculate density
      // oil
      if (is_o)
        {
          if (is_g)
            {
              // Kamaltinova A.
              if (FI_CHK_SG (main_vars_, i))
                {
                  gas_oil_ratio_[i]             = gor;
                  data_i.p_deriv_gas_oil_ratio  = d_gor;
                  data_i.density[d_o]           = data_i.invers_fvf[d_o] * (oil_surface_density + gor * gas_surface_density);
                  //data_i.p_deriv_density[d_o]   = data_i.p_deriv_invers_fvf[d_o] * oil_surface_density + gas_surface_density;
                  data_i.p_deriv_density[d_o]   = data_i.p_deriv_invers_fvf[d_o] * oil_surface_density + gas_surface_density * (gas_oil_ratio_[i] * data_i.p_deriv_invers_fvf[d_o] + data_i.invers_fvf[d_o] * data_i.p_deriv_gas_oil_ratio);
                  data_i.gor_deriv_density      = 0.0;
                }
              else if (FI_CHK_RO (main_vars_, i))
                {
                  data_i.density[d_o]           = data_i.invers_fvf[d_o] * (oil_surface_density + gas_oil_ratio_[i] * gas_surface_density);
                  data_i.p_deriv_gas_oil_ratio  = 0.0;
                  data_i.p_deriv_density[d_o]   = data_i.p_deriv_invers_fvf[d_o] * (oil_surface_density + gas_oil_ratio_[i] * gas_surface_density);
                  data_i.gor_deriv_density      = data_i.gor_deriv_invers_fvf * (oil_surface_density + gas_oil_ratio_[i] * gas_surface_density) + data_i.invers_fvf[d_o] * gas_surface_density;
                }
              else if (FI_CHK_MOMG (main_vars_, i))
                {
                  gas_oil_ratio_[i]             = gor;
                  data_i.p_deriv_gas_oil_ratio  = d_gor;
                  data_i.density[d_o]           = data_i.invers_fvf[d_o] * oil_surface_density;
                }
            }
          else
            {
              data_i.density[d_o]               = data_i.invers_fvf[d_o] * oil_surface_density;
              data_i.p_deriv_density[d_o]       = data_i.p_deriv_invers_fvf[d_o] * oil_surface_density;
            }
        }

      // calculate water density
      if (is_w)
        {
          data_i.density[d_w]           = data_i.invers_fvf[d_w] * water_surface_density;
          data_i.p_deriv_density[d_w]   = data_i.p_deriv_invers_fvf[d_w] * water_surface_density;
        }
      // calculate gas density
      if (is_g)
        {
          data_i.density[d_g]           = data_i.invers_fvf[d_g] * gas_surface_density;
          data_i.p_deriv_density[d_g]   = data_i.p_deriv_invers_fvf[d_g] * gas_surface_density;
        }

      //Kamaltinova A.
      // calculate mobility
      //water phase
      if (is_w)
        {
          // mobility = 1/(B*mu)
          data_i.mobility[d_w]                = data_i.invers_visc_fvf[d_w];
          data_i.p_deriv_mobility[d_w]        = data_i.p_deriv_invers_visc_fvf[d_w];
          if (n_phases > 1)
            {
              // mobility = mobility*k'
              data_i.mobility[d_w]            *= data_i.relative_perm[d_w];
              data_i.p_deriv_mobility[d_w]    *= data_i.relative_perm[d_w];
              data_i.s_deriv_mobility[d_ww]   = data_i.s_deriv_relative_perm[d_ww] * data_i.invers_visc_fvf[d_w]
                                              + data_i.relative_perm[d_w] * data_i.p_deriv_invers_visc_fvf[d_w] * data_i.s_deriv_cap_pressure[ds_w];
              if (is_o)
                data_i.s_deriv_mobility[d_wo] = 0.0;
              if (is_g)
                data_i.s_deriv_mobility[d_wg] = 0.0;
            }
        }
      //gas phase
      if (is_g)
        {
          data_i.mobility[d_g]                = data_i.invers_visc_fvf[d_g];
          data_i.p_deriv_mobility[d_g]        = data_i.p_deriv_invers_visc_fvf[d_g];
          if (n_phases > 1)
            {
              data_i.mobility[d_g]            *= data_i.relative_perm[d_g];
              data_i.p_deriv_mobility[d_g]    *= data_i.relative_perm[d_g];

              if (FI_CHK_SG (main_vars_, i))
                {
                  data_i.s_deriv_mobility[d_gg] = data_i.s_deriv_relative_perm[d_gg] * data_i.invers_visc_fvf[d_g]
                                                + data_i.relative_perm[d_g] * data_i.p_deriv_invers_visc_fvf[d_g] * data_i.s_deriv_cap_pressure[ds_g];
                }
              else
                {
                  data_i.s_deriv_mobility[d_gg] = 0;
                  data_i.mobility[d_g] = 0;
                  data_i.p_deriv_mobility[d_g] = 0;
                }

              data_i.s_deriv_mobility[d_go]   = 0;

              if (is_w)
                data_i.s_deriv_mobility[d_gw] = 0;

              if (FI_CHK_MOMG (main_vars_, i))
                {
                  data_i.mobility[d_g]        = 0;
                  data_i.p_deriv_mobility[d_g] = 0;
                }
            }
        }
      //oil phase
      if (is_o)
        {
          data_i.mobility[d_o]                = data_i.invers_visc_fvf[d_o];
          data_i.p_deriv_mobility[d_o]        = data_i.p_deriv_invers_visc_fvf[d_o];
          if (n_phases > 1)
            {
              data_i.mobility[d_o]            *= data_i.relative_perm[d_o];
              data_i.p_deriv_mobility[d_o]    *= data_i.relative_perm[d_o];

              data_i.s_deriv_mobility[d_oo]   = data_i.s_deriv_relative_perm[d_oo] * data_i.invers_visc_fvf[d_o];

              if (is_w)
                data_i.s_deriv_mobility[d_ow] = data_i.s_deriv_relative_perm[d_ow] * data_i.invers_visc_fvf[d_o];

              if (is_g)
                {
                  if (FI_CHK_SG (main_vars_, i))
                    data_i.s_deriv_mobility[d_og] = data_i.s_deriv_relative_perm[d_og] * data_i.invers_visc_fvf[d_o];
                  else if (FI_CHK_RO (main_vars_, i))
                    data_i.s_deriv_mobility[d_og] = data_i.relative_perm[d_o] * data_i.gor_deriv_invers_visc_fvf;
                  else if (FI_CHK_MOMG (main_vars_, i))
                    {
                      data_i.mobility[d_o] = 0;
                      data_i.p_deriv_mobility[d_o] = 0;
                      data_i.s_deriv_mobility[d_oo] = 0;
                      data_i.s_deriv_mobility[d_og] = 0;
                      if (is_w)
                        data_i.s_deriv_mobility[d_ow] = 0;
                    }
                }
            }
        }

#ifdef FI_OPERATOR_CELLS_PARALLEL
      if (is_w)
        invers_fvf_average_w += invers_fvf[i_w];
      if (is_g)
        invers_fvf_average_g += invers_fvf[i_g];
      if (is_o)
        invers_fvf_average_o += invers_fvf[i_o];
#else //FI_OPERATOR_CELLS_PARALLEL

#ifdef _MPI
      if (i >= n_left && i < n_own)
        {
#endif //_MPI
          if (is_w)
            calc_model_->invers_fvf_average[d_w] += data_i.invers_fvf[d_w];
          if (is_g)
            calc_model_->invers_fvf_average[d_g] += data_i.invers_fvf[d_g];
          if (is_o)
            calc_model_->invers_fvf_average[d_o] += data_i.invers_fvf[d_o];
#ifdef _MPI
        }
#endif //_MPI

#endif //FI_OPERATOR_CELLS_PARALLEL

    }

    /**
     * \brief  Checks norm
     * \param  istart
     * \return True if norm smaller than threshold
     * */
    bool
    check_norm (index_t istart)
    {
      if (!istart)
        {
#ifdef _DEBUG
          BOSOUT (section::iters, level::debug) << "rhs_residual: " << fabs (norm_.val[norms::C_CPV]) << " / " << rhs_residual_ << bs_end;
          BOSOUT (section::iters, level::debug) << "mb_error: " << fabs (norm_.val[norms::MB_ERR]) << " / " << mb_error_ << bs_end;
          BOSOUT (section::iters, level::debug) << "s_rhs: " << fabs (norm_.val[norms::S_RHS]) << " / " << s_rhs_norm_ << bs_end;
#endif

          return fabs (norm_.val[norms::C_CPV]) < rhs_residual_
            && fabs (norm_.val[norms::MB_ERR]) < mb_error_
            //&& fabs (norm_.val[norms::S_RHS]) < s_rhs_norm_
            ;
        }

      return false;
    }
    /**
     * \brief  Checks multiplier
     * \todo   Obsolete
     * */
    bool
    check_solution_mult_cell (index_t istart, index_t iteration, index_t n_approx, const norms_storage_t &base_norm)
    {
      // if we return false here - we got smaller number of newtonian iterations
      // if we compute solution_mult_cell - we got bigger number of newtonian iterations
      return false;

      if (!istart && !iteration)
        {
          BOSOUT (section::iters, level::debug) << "check_solution_mult_cell" << bs_end;

          item_t mult = calc_solution_mult_cell (base_norm);
          if (mult < 0.9 && iteration < n_approx - 1)
            {
              if (mult < item_t (0.1))
                mult = item_t (0.1);

              BOSOUT (section::iters, level::debug) << "check_solution_mult_cell 2, mult = " << mult << bs_end;

              restore_prev_niter_vars ();
              if (calc_model_->apply_newton_correction (mult, 2, mesh_, jmatrix_))
                bs_throw_exception ("apply_newton_correction failed");

              return true;
            }
#ifdef _DEBUG
          else
            {
              BOSOUT (section::iters, level::debug) << "check_solution_mult_cell 3, mult = " << mult << bs_end;
            }
#endif
        }

      return false;
    }

    /**
     * \brief       Reduces dt and multiplies Jacobian with new dt value
     * \param[out]  dt
     * \param[in]   istart
     * */
    void
    do_dt_reduce (double &dt, index_t istart)
    {
      if (/*false && */istart)
        {
          item_t dt_mult = calc_step_dt_mult (0, (item_t) calc_model_->ts_params->get_float (fi_params::MAX_NORM_ON_FIRST_N));
          if (dt_mult < 0.95)
            {
              jmatrix_->mult_flux_part (dt_mult);
              dt *= dt_mult;

              BOSOUT (section::iters, level::medium) << "dT reduced by factor " << dt_mult << ": " << dt << bs_end;
            }
          dt_mult = 0;
        }
    }
    /**
     * \brief      Multiplies Jacobian with new calculated dt value
     * \param[in]  dt
     * \param[in]  istart
     * */
    void
    do_dt_tuning (const double &dt, bool tuning)
    {
      if (tuning)
        {
          item_t dt_mult = calc_step_dt_mult (0, (item_t) calc_model_->ts_params->get_float (fi_params::MAX_NORM_ON_TS));
          if (dt_mult < 0.95)
            {
              jmatrix_->mult_flux_part (dt_mult);

              BOSOUT (section::iters, level::medium) << "Use dT " << dt_mult * (dt) << " instead of original dT " << dt << bs_end;
            }
          else
            dt_mult = 1.0;
        }
    }

    /**
     * \brief  Saves data before start of new newton iteration
     * */
    void
    save_prev_niter_vars ()
    {
      int is_line_search = calc_model_->ts_params->get_int (fi_params::SELECT_SOL_STEPS);
      if (is_line_search != 0)
        {
          calc_model_->prev_niter_data_.save (calc_model_);
          reservoir_->pre_newton_step ();
        }
    }

    /**
     * \brief  Restores data if newton iteration failed
     * */
    void
    restore_prev_niter_vars ()
    {
      reservoir_->restart_newton_step ();
    }

    /**
     * \brief Calculates norms
     * */
    void
    norm_calc ();

    /**
     * \brief  Calculates norm in one cell and update norm storage
     * \param  i Index of cell
     * \param  ns Norms storage
     * */
    void
    update_norm_by_cell (index_t i, norms_storage_t &ns);

    /**
     * \brief  Calculates dt at the first newton iteration on time step
     * \param  prev_mult Previous dt multiplier
     * \param  max_res
     * \return New dt multiplier
     * */
    item_t
    calc_step_dt_mult (item_t prev_mult, item_t max_res);

    /**
     * \brief      Calculates porosity and derivativies, also calculates 
     *             trunsmissibility multipliers
     *
     * \param[in]  i Cell index
     * \param[in]  pvt_reg PVT region index
     * \param[out] poro Calculated porosity
     * \param[out] dp_poro Calculated porosity derivative
     * \param[out] t_mult  Calculated trunsmissibility multiplier
     * \param[out] dp_t_mult Calculated trunsmissibility multiplier derivative
     * */
    void
    calc_porosity_and_deriv (index_t i,
                             index_t pvt_reg,
                             item_t *poro,
                             item_t *dp_poro,
                             item_t *t_mult,
                             item_t *dp_t_mult);

    /**
     * \brief  Calculates fluid volume on previous step
     * */
    void
    calc_prev_fluid_volume ();

    /**
     * \brief  Saves debug data
     * \param  dt
     * */
    void
    debug_save_data (item_t dt);

    /**
     * \brief  Calculates multipler for solution
     * \param  base_norm Norm storage
     * \return Calculated multiplier
     * */
    item_t
    calc_solution_mult_cell (const norms_storage_t &base_norm);

    /**
     * \brief  Calculates flux part of Jacobian
     * \param  dt
     * \return True
     * */
    bool
    block_connections_mpfa (const item_t &dt);

  public:
    sp_calc_model_t               &calc_model_;
    const sp_rock_grid_prop_t     &rock_grid_prop_;
    sp_reservoir_t                &reservoir_;
    sp_jacobian_t                 &jacobian_;
    sp_jmatrix_t                  &jmatrix_;
    const sp_mesh_iface_t         &mesh_;

    sp_bcsr_matrix_t              trns_matrix_;
    const rhs_item_array_t        &trns_values_;
    const index_array_t           &trns_rows_ptr_;
    const index_array_t           &trns_cols_ptr_;

    sp_bcsr_matrix_t              reg_matrix_;
    rhs_item_array_t              &reg_values_;
    const index_array_t           &reg_rows_ptr_;
    const index_array_t           &reg_cols_ptr_;

    const index_array_t           &m_array_;
    const index_array_t           &p_array_;

    rhs_item_array_t              &rhs_;
    item_array_t                  &sol_;
    rhs_item_array_t              &flux_rhs_;
    rhs_item_array_t              &sp_diag_;
    rhs_item_array_t              &s_rhs_;

    const item_array_t            &depths_;

    index_t                       n_cells_;
    index_t                       n_connections_;

    index_t                       n_sec_vars;

    index_t                       d_w;
    index_t                       d_g;
    index_t                       d_o;
    index_t                       ds_w;
    index_t                       ds_g;

    index_t                       d_gg, d_gw, d_go;
    index_t                       d_wg, d_ww, d_wo;
    index_t                       d_og, d_ow, d_oo;

    data_array_t                  &data_;
    item_array_t                  &saturation_3p_;
    const item_array_t            &pressure_;
    item_array_t                  &gas_oil_ratio_;
    main_var_array_t              &main_vars_;
    const item_array_t            &volume_;
    const item_array_t            &poro_array_;
    const item_array_t            &rock_grid_comp_const_;
    const item_array_t            &rock_grid_comp_ref_pressure_;
    const index_array_t           &sat_regions_;
    const index_array_t           &pvt_regions_;

    const sp_pvt_dead_oil_array_t &pvt_oil_array;                  //!< (n_pvt_regions)
    const sp_pvt_water_array_t    &pvt_water_array;
    const sp_pvt_gas_array_t      &pvt_gas_array;

    item_t                        min_p_;
    item_t                        max_p_;
    item_t                        drsdt_;
    item_t                        rhs_residual_;
    item_t                        mb_error_;
    item_t                        s_rhs_norm_;

    norms_storage_t               &norm_;
    rhs_item_array_t              &cfl_;
  };

} // namespace blue_sky

#include "fi_operator_norm_calc.h"
#include "fi_operator_calc_step_dt_mult.h"
#include "fi_operator_calc_porosity_and_deriv.h"
#include "fi_operator_calc_prev_fluid_volume.h"
#include "fi_operator_calc_solution_mult_cell.h"

#ifdef BS_USE_TPFA_
#include "fi_operator_block_connections_mpfa.h"
#else
#include "fi_operator_block_connections_mpfa_2.h"
#endif

#endif  // #ifndef BS_FI_OPERATOR_H_

