/**
 *       \file  main_loop_calc.h
 *      \brief  Main calculation loop implementation
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  04.09.2008
 *  \copyright  This source code is released under the terms of
 *              the BSD License. See LICENSE for more details.
 *       \todo  Split and move implementation to src/
 *       \todo  Describe data members
 * */
#ifndef BS_MAIN_LOOP_CALC_H_
#define BS_MAIN_LOOP_CALC_H_

#include "data_storage_interface.h"
#include "event_manager.h"

#include "string_formater.h"
#include "well_type_helper.h"
#include "well_connection.h"
#include "calc_well.h"
#include "rr_rw_wr_saver.h"

#include "fi_operator.h"
//#include "well_results_storage.h"
//#include "fip_results_storage.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "scale_array_holder.h"
#include "print_process_memory_info.h"
#include "path_tools.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include "well_reporter.h"

namespace blue_sky
  {

    /**
     * \class main_loop_calc_base
     * \brief Base interface for main_loop calculation
     *        implementation
     * \todo  Should be renamed to main_loop_calc_iface
     * */
  struct main_loop_calc_base
  {
    typedef event_manager event_manager_t;
    typedef event_manager_t::event_list_t  event_list_t;

  public:
    /**
     * \brief  dtor
     * */
    virtual ~main_loop_calc_base () {}

    //! Prepares calculation
    virtual void ready () = 0;
    //! Starts calculation loop
    virtual void go ()    = 0;
    //! Ends calculation
    virtual void end ()   = 0;

    //! Applies model events on each step
    /**
     * \brief  Applies model events on each step
     * \param  event_list
     * */
    virtual void
    apply_events (const event_list_t &) = 0;
  };

  /**
   * \class main_loop_calc
   * \brief Main calculation loop implementation
   * */
  template <bool is_w, bool is_g, bool is_o>
  struct main_loop_calc : public main_loop_calc_base 
    {
      typedef main_loop_calc_base               base_t;
      typedef main_loop_calc <is_w, is_g, is_o> this_t;

      typedef t_float                           item_t;
      typedef t_long                            index_t;
      typedef spv_float                         item_array_t;

      typedef calc_model                        calc_model_t;
      typedef event_manager                     event_manager_t;
      typedef rock_grid                         rock_grid_t;
      typedef reservoir                         reservoir_t;
      typedef rs_mesh_iface                     mesh_iface_t;
      typedef jacobian                          jacobian_t;
      typedef reservoir_simulator               reservoir_simulator_t;
      typedef idata                             idata_t;

      typedef event_manager_t::event_list_t     event_list_t;
      typedef trans_multipliers_calc            trans_multipliers_calc_t;

      typedef smart_ptr <calc_model_t, true>                sp_calc_model_t;
      typedef smart_ptr <event_manager_t, true>             sp_event_manager_t;
      typedef smart_ptr <rock_grid_t, true>                 sp_rock_grid_t;
      typedef smart_ptr <reservoir_t, true>                 sp_reservoir_t;
      typedef smart_ptr <mesh_iface_t, true>                sp_mesh_iface_t;
      typedef smart_ptr <jacobian_t, true>                  sp_jacobian_t;
      typedef smart_ptr <fi_params, true>                   sp_fi_params_t;
      typedef smart_ptr <data_storage_interface, true>      sp_storage_t;
      typedef smart_ptr <reservoir_simulator_t, true>       sp_rs_t;
      typedef smart_ptr <idata_t, true>                     sp_idata_t;

      typedef boost::posix_time::ptime                      ptime;

    public:
      /**
       * \brief  ctor
       * \param  rs Instance of reservoir_simulator to which
       *            calculation will be performed
       * */
      main_loop_calc (const sp_rs_t &rs)
          : rs_ (rs)
          , calc_model_ (rs->get_calc_model ())
          , event_manager_ (rs->get_event_manager ())
          , rock_grid_prop_ (calc_model_->rock_grid_prop)
          , facility_storage_ (rs->get_facility_storage ())
          , reservoir_ (rs->get_reservoir ())
          , mesh_ (rs->get_hdm ()->get_mesh ())
          , jacobian_ (rs->get_jacobian ())
          , params_ (calc_model_->ts_params)
          , data_map_ (rs->get_hdm ()->get_data ())
          , height_ (0)
          , rho_ (0)
      {
        number_of_newtonian_iterations  = 0;
        number_of_linear_iterations     = 0;
        number_of_restarts              = 0;
        number_of_max_iters_restarts    = 0;
        number_of_fi_operator_restarts  = 0;
        number_of_close_wells_restarts  = 0;

        large_time_step_length_         = 0;
        large_time_step_start_          = 0;
        large_time_step_end_            = 0;

        num_last_newton_iters_          = 0;
        num_last_lin_iters_             = 0;

        number_of_small_time_steps      = 0;
        number_of_large_time_steps      = 1;

        height_                         = 0;
        rho_                            = 0;
        min_z_                          = 0;
        max_z_                          = 0;
        dt_                             = 0;
        current_time_                   = 0;

        do_calc_prev_fluid_volume_      = true;

        if (calc_model_->n_phases > 1)
          {
            if (calc_model_->phase_d[FI_PHASE_WATER] != calc_model_->sat_d[FI_PHASE_WATER])
              bs_throw_exception ((boost::format ("Phase index (%d) for water should be equal to saturation index (%d)") % calc_model_->phase_d[FI_PHASE_WATER] % calc_model_->sat_d[FI_PHASE_WATER]).str ());

            if (calc_model_->phase_d[FI_PHASE_GAS] != calc_model_->sat_d[FI_PHASE_GAS])
              bs_throw_exception ((boost::format ("Phase index (%d) for gas should be equal to saturation index (%d)") % calc_model_->phase_d[FI_PHASE_GAS] % calc_model_->sat_d[FI_PHASE_GAS]).str ());
          }
      }

      /**
       * \brief  Inits fi_operator
       * \param  i
       * */
      inline void
      fi_operator_cells (index_t i)
      {
        fi_operator_impl <is_w, is_g, is_o> fi_operator_impl_ (calc_model_, 
          reservoir_, 
          mesh_, 
          jacobian_); 
        fi_operator_impl_.fi_operator_init (i, dt_);
      }

      /**
       * \brief  Applies model events
       * \param  event_list
       * */
      inline void
      apply_events (const event_list_t &event_list)
      {
        typename event_list_t::const_iterator it = event_list.begin ();
        typename event_list_t::const_iterator e  = event_list.end ();

        for (; it != e; ++it)
          {
            (*it)->apply (reservoir_, mesh_, calc_model_, data_map_);
          }
      }

      /**
       * \brief  Returns new time-step value
       * \return New time-step value
       * */
      inline item_t
      get_dt () const
        {
          //if (ts_it != ts_it_end && (int)ts_it->size () > small_time_step_num && use_timestep_file)
          //  {
          //    dt = (*ts_it)[small_time_step_num].step;
          //    if (t + dt > large_time_step_end)
          //      dt = large_time_step_end - t;
          //  }
          //else if (ts->ts_params.check_d_params (FI_PARAMS_D_FIRST_TS)
          //          || (small_time_step_num == 0 && !l_time_step_num))

          //  {
          //    dt = model_base->ts_params.get_d_param (FI_PARAMS_D_FIRST_TS);
          //    printf ("hello\n");
          //  }
          //else
          //  dt = model_base->increase_ts (dt, large_time_step_end - t, nniters);

          // miryanov
          if (/*get_first_ts () || */(number_of_small_time_steps == 0 && !number_of_large_time_steps))
            {
              BOSOUT (section::app_info, level::debug) << "number_of_small_time_steps == 0 && !number_of_large_time_steps" << bs_end;
              return (std::min <item_t>) (1.0/*get_first_ts ()*/, (item_t)large_time_step_end_);
            }
          else
            {
              BOSOUT (section::main_loop, level::low) << "increase_ts" << bs_end;
              return increase_ts (dt_, large_time_step_end_ - current_time_, num_last_newton_iters_);
            }
        }

      /**
       * \brief  Small-step loop
       * */
      inline void
      process_small_steps ()
      {
        setup_jacobian_solver_params ();

        fi_operator_impl <is_w, is_g, is_o> fi_operator (calc_model_, reservoir_, mesh_, jacobian_);

        for (number_of_small_time_steps = 0; current_time_ < large_time_step_end_ - EPS_DIFF; ++number_of_small_time_steps)
          {
            do_calc_prev_fluid_volume_ = true;
            dt_ = get_dt ();
            save_base_vars ();

            for (;;)
              {
                // no neighbours
                calc_model_->approx_flag = false;
                if (get_clamp_pressure ())
                  compute_solution_range ();

                rs_->on_small_step_start ();
                index_t ret_code = compute_small_time_step (fi_operator, num_last_newton_iters_, num_last_lin_iters_);
                if (ret_code >= 0 && ret_code <= get_newton_iters_num ())
                  {
                    if (!process_well_model ())
                      break;
                    else
                      {
                        do_calc_prev_fluid_volume_ = false;
                        continue;
                      }
                  }
                else if (ret_code >= 0)
                  {
                    newton_iter_fail(ret_code);
                  }
              }

            newton_iter_success();
          }
      }

      /**
       * \brief  Called if newton iteration failed
       * \param  ret_code Return code from newton iteration
       * */
      BS_FORCE_INLINE
      void
      newton_iter_fail (size_t ret_code)
      {
        restore_base_vars ();
        item_t old_dt = dt_;
        dt_ = decrease_ts (dt_, large_time_step_end_ - current_time_);

        BOSOUT (section::main_loop, level::low) << "have to decrease dt " << old_dt << " -> " << (double)dt_ << bs_end;
        if (fabs (old_dt - dt_) < 1.0e-12)
          {
            throw bs_exception ("process_small_step", "Newton process failed");
          }

        if (ret_code > get_newton_iters_num ())
          {
            update_number_of_fi_operator_restarts ();
            update_number_of_restarts ();
          }
        else if (ret_code == get_newton_iters_num ())
          {
            update_number_of_max_iters_restarts ();
            update_number_of_restarts ();
          }

        do_calc_prev_fluid_volume_ = false;
        rs_->on_newton_iter_fail ();
       }
      /**
       * \brief  Called if newton iteration successed
       * */
      BS_FORCE_INLINE
      void
      newton_iter_success ()
      {
        update_number_of_newtonian_iterations (num_last_newton_iters_);
        update_number_of_linear_iterations (num_last_lin_iters_);
        update_current_time (dt_);

        fi_borates_total ();

        //save_newton_iter_info ();

        //print_saturation ();
        //print_pressure ();
        //print_volume ();


        BOSOUT (section::main_loop, level::high)
          << (boost::format ("Info: Newtonian iteration [%d], linear interations [%d]") % (size_t)num_last_newton_iters_ % (size_t)num_last_lin_iters_).str () << "\n"
          << (boost::format ("      dt: [%10.20lf] on time [%10.20lf]") % (double)dt_ % (double)current_time_).str () << bs_end;

        rs_->on_newton_iter_success ();
      }

      /**
       * \brief  Saves info about newton iteration
       * \todo   Obsolete, should be removed
       * */
      inline void
      save_newton_iter_info ()
      {
        static FILE *file = fopen ("newton_iter_info.bs.txt", "wt");
        if (!file)
          {
            throw bs_exception ("can't open file", "");
          }

        fprintf (file, "newtonian iters: %d, linear iters: %d, dt = %.20e on time = %.20e\n", (int)num_last_newton_iters_, (int)num_last_lin_iters_, (double)dt_, (double)current_time_);
        fflush (file);
      }

      /**
       * \brief  Saves saturation values to file
       * \todo   Obsolete, should be removed
       * */
      void
      print_saturation ()
      {
        static size_t iter_counter = 0;
        ++iter_counter;
        //tools::save_seq_vector (tools::string_formater ("nw_saturation.bs.%d.txt", iter_counter).str).save (calc_model_->saturation_3p);
      }
      /**
       * \brief  Saves pressure values to file
       * \todo   Obsolete, should be removed
       * */
      void
      print_pressure ()
      {
        static size_t iter_counter = 0;
        ++iter_counter;
        //tools::save_seq_vector (tools::string_formater ("nw_pressure.bs.%d.txt", iter_counter).str).save (calc_model_->pressure);
      }

      /**
       * \brief  Saves volume values to file
       * \todo   Obsolete, should be removed
       * */
      void
      print_volume ()
      {
        static size_t iter_counter = 0;
        ++iter_counter;
        //tools::save_seq_vector (tools::string_formater ("nw_volumes.bs.%d.txt", iter_counter).str).save (calc_model_->rock_grid_prop->volume);
      }

      inline bool
      process_well_model ()
      {
        if (calc_model_->well_model_var_ == WELL_MODEL_1VAR)
          {
            return process_well_model_1var ();
          }
        if (calc_model_->well_model_var_ == WELL_MODEL_3VAR)
          {
            return process_well_model_3var ();
          }

        return false;
      }

      inline bool
      process_well_model_1var ()
      {
        int ret_code = fi_borates_check_well_consistensy ();
        if (ret_code < 0)
          {
            restore_base_vars ();

            //TODO: LOG
            update_number_of_restarts ();
          }
        else if (ret_code > 0)
          {
            restore_base_vars ();

            item_t old_dt = dt_;
            dt_ = decrease_ts (dt_, large_time_step_end_ - current_time_);

            if (fabs (old_dt - dt_) < 1.0e-12)
              {
                throw bs_exception ("process_well_model_1var", "Newton process failed");
              }

            update_number_of_restarts ();
          }

        return ret_code != 0;
      }

      inline bool
      process_well_model_3var ()
      {
        if (false)
          {
            restore_base_vars ();

            // TODO: LOG

            update_number_of_close_wells_restarts ();
            update_number_of_restarts ();
            return true;
          }

        return false;
      }

      inline void
      fi_borates_total ()
      {
        // TODO: MPI staff
        compute_acc_rates ();
      }

      inline index_t
      fi_borates_check_well_consistensy ()
      {
        BS_ASSERT (false && "NOT IMPL YET");
        return 0;
      }

      inline void
      set_approx_solution ()
      {
        BS_ASSERT (false && "NOT IMPL YET. IN OLD CODE TOO.");
      }

      /**
       * \brief  Returns number of maximum newton iterations
       *         (fi_params::NEWTON_ITERS_NUM)
       * \return Number of newton iterations
       * */
      inline item_t
      get_newton_iters_num () const
        {
          return params_->get_int (fi_params::NEWTON_ITERS_NUM);
        }

      /**
       * \brief  Returns maximum (?) pressure value
       * \return Maximum (?) pressure value
       * */
      inline item_t
      get_clamp_pressure () const
        {
          return params_->get_bool (fi_params::CLAMP_PRESSURE);
        }

      /**
       * \brief  Returns maximum linear solver tolerance
       *         (fi_params::MAX_ALLOWED_LIN_SOLV_RESID)
       * \return Maximum linear solver tolerance
       * */
      inline item_t
      get_max_tolerance () const
        {
          return params_->get_float (fi_params::MAX_ALLOWED_LIN_SOLV_RESID);
        }

      /**
       * \brief  Returns minimum value of time-step
       *         (fi_params::MIN_TS)
       * \return Minimum value of time-step
       * */
      inline item_t
      get_small_ts () const
        {
          return 10 * params_->get_float (fi_params::MIN_TS);
        }
      /**
       * \brief  Return value of first time-step
       *         (fi_params::FIRST_TS)
       * \return Value of first time-step
       * */
      inline item_t
      get_first_ts () const
        {
          BOSOUT (section::main_loop, level::low) << "get_first_ts: " << params_->get_float (fi_params::FIRST_TS) << bs_end;
          return params_->get_float (fi_params::FIRST_TS);
          //return params_->check_value (fi_params::FIRST_TS);
        }

      /**
       * \brief  Returns number of maximum newton iterations
       *         (fi_params::APPROX_STEPS or fi_params::NEWTON_ITERS_NUM)
       * \return Number of maximum newton iterations
       * */
      inline index_t
      get_n_max_iters () const
      {
        if (calc_model_->approx_flag)
          return params_->get_int (fi_params::APPROX_STEPS);
        else
          return params_->get_int (fi_params::NEWTON_ITERS_NUM);
      }
      inline index_t
      get_istart () const
      {
        return 1;
      }

      /**
       * \brief Setups and solves jacobians, also raises related events. Increment nliters if tolerance ok.
       * \param nliters
       * \return true if tolerance <= [solver->tolerance, this->get_max_tolerance]
       * */
      bool
      solve_jacobian (t_long &nliters)
      {
        spv_int filter = BS_KERNEL.create_object (v_int::bs_type ());
        BS_SP (bcsr_matrix_iface) m = jacobian_->get_matrix ()->merge (filter);
        BS_SP (lsolver_iface) solver = jacobian_->get_solver ();
        rs_->on_before_jacobian_setup ();
        if (solver->setup (m))
          {
            bs_throw_exception ("Can't setup solver");
          }

        rs_->on_before_jacobian_solve ();
        t_long solved = solver->solve (m, jacobian_->get_rhs_dbl (), jacobian_->get_solution ());
        t_double tolerance = 0;
        t_long n_current_liters = 0;

        // FIXME: should we restore secondary solution?
        jacobian_->restore_sec_solution ();

        bool valid_tolerance = false;
        if (solved < 0)
          {
            BOSERR (section::solvers, level::error) << "Linear solver failed with retcode = " << solved << bs_end;
            // FIXME: we just return false
            //tolerance = 10.0f;
          }
        else
          {
            n_current_liters = solver->get_prop ()->get_i (iters_idx);
            tolerance = solver->get_prop ()->get_f (final_res_idx);
            if (tolerance > solver->get_prop ()->get_f (tol_idx))
              {
                BOSERR (section::solvers, level::error)
                  << "Linear solver failed with tolerance "
                  << tolerance
                  << bs_end;
              }
            else
              {
                BOSOUT (section::solvers, level::medium)
                  << "Linear solver iterations " 
                  << n_current_liters 
                  << ", tol = " << tolerance
                  << ", ret_code = " << solved 
                  << bs_end;

                valid_tolerance = tolerance <= get_max_tolerance ();
                if (valid_tolerance)
                  {
                    nliters += n_current_liters;
                  }
              }
          }

        return valid_tolerance;
      }

      /**
       * \brief  Calculates one small step
       * \param  fi_operator Instance of fi_operator_impl
       * \param  nniters
       * \param  nliters
       * \return Number of newton iteration was performed
       * */
      inline index_t
      compute_small_time_step (fi_operator_impl <is_w, is_g, is_o> &fi_operator,
        index_t &nniters, index_t &nliters)
      {
        index_t max_n_iters = get_n_max_iters ();
        index_t istart = get_istart ();
        nniters = 0;
        nliters = 0;

        for (index_t i = 0; i < max_n_iters + 1; ++i, ++nniters)
          {
            int init = (i == 0) ? istart : 0;

            if (do_calc_prev_fluid_volume_ && init)
              fi_operator.calc_prev_fluid_volume ();

            if (istart && !number_of_small_time_steps)
              {
                set_approx_solution ();
              }

            if (false)
              generate_numeric_jacobian (init);

            rs_->on_before_fi_operator ();
            fi_operator_return_type finish_flag = fi_operator.fi_operator (dt_, init, init, number_of_small_time_steps, true, false);
            if (finish_flag == FI_OPERATOR_RETURN_RESTART || i == max_n_iters)
              {
                return max_n_iters + 1;
              }
            else if (finish_flag == FI_OPERATOR_RETURN_FAIL)
              {
                return i;
              }

            if (!solve_jacobian (nliters))
              {
                return max_n_iters + 1;
              }

            rs_->on_before_restore_solution ();
            fi_operator.save_prev_niter_vars ();
            restore_solution_return_type ret_code = calc_model_->restore_solution (fi_operator.mesh_, 
              jacobian_->get_solution (), 
              jacobian_->get_sec_solution ());

            if (ret_code == SMALL_TIME_STEP_CHOP)
              {
                check_time_step ();
              }
            else if (ret_code == SMALL_TIME_STEP_FAIL)
              {
                return max_n_iters + 1;
              }

            reservoir_->restore_wells_solution (dt_, fi_operator.sol_, jacobian_->get_sec_solution (), calc_model_->n_phases);
          }

        return max_n_iters + 1;
      }

      void print_memory_at_the_end ()
      {
#ifdef BS_BOS_CORE_DEBUG_MEMORY
        BOSOUT (section::app_info, level::debug) << "at the end" << bs_end;
        BS_KERNEL.get_memory_manager ().print_info (false);
#endif
      }
      void print_memory_before_jacobian_solve ()
      {
#ifdef BS_BOS_CORE_DEBUG_MEMORY
        BOSOUT (section::app_info, level::debug) << "before solve_jacobian" << bs_end;
        BS_KERNEL.get_memory_manager ().print_info (false);
#endif
      }
      void print_memory_before_jacobian_setup ()
      {
#ifdef BS_BOS_CORE_DEBUG_MEMORY
        BOSOUT (section::app_info, level::debug) << "before setup_jacobian" << bs_end;
        BS_KERNEL.get_memory_manager ().print_info (false);
#endif
      }
      void print_memory_before_fi_operator ()
      {
#ifdef BS_BOS_CORE_DEBUG_MEMORY
        BOSOUT (section::app_info, level::debug) << "before fi_operator" << bs_end;
        BS_KERNEL.get_memory_manager ().print_info (false);
#endif
      }

      /**
       * \brief  Fills jacobian with derivs that numerically calculated
       * \param  init
       * \todo   Obsolete
       * */
      inline void
      generate_numeric_jacobian (int init);

      /**
       * \brief  Calculates total rates for well and
       *         rate and total rate for reservoir
       * */
      inline void
      compute_acc_rates ();

      /**
       * \brief  Checks time-step
       * \return Throws exception is time-step smaller
       *         than threshold
       * */
      inline void
      check_time_step ()
      {
        if (dt_ < get_small_ts ())
          {
            throw bs_exception ("compute_small_step", "Non-Linear Equation Convergence Failure");
          }
      }

      /**
       * \todo   Obsolete
       * */
      inline void
      compute_solution_range ()
      {
        get_min_max_z ();                                   // TODO: should be called only once

        height_ = max_z_ - min_z_;                          // TODO: should be called only once
        rho_    = calc_model_->get_initial_rho (height_);   // TODO: should be called only once

        BS_ASSERT (min_z_ <= max_z_) (min_z_) (max_z_);
        calc_model_->update_min_pressure_range (min_z_ - 10);
        calc_model_->update_max_pressure_range (max_z_ + 10);
      }

      /**
       * \brief  Setups Jacobian solvers params
       * */
      inline void
      setup_jacobian_solver_params ()
      {
        jacobian_->setup_solver_params (calc_model_->well_model_type_, calc_model_->n_phases, params_);
      }
      ///**
      // * \brief  Setups Jacobian
      // * */
      //inline void
      //setup_jacobian ()
      //{
      //  if (jacobian_->setup ())
      //    {
      //      throw bs_exception ("compute_small_time_step", "return -1");
      //    }
      //}
      ///**
      // * \brief  Solves Jacobian
      // * */
      //inline item_t
      //solve_jacobian (index_t &n)
      //{
      //  return jacobian_->solve (n);
      //}

      inline void
      get_min_max_z ()
      {
        item_t foo;
        mesh_->get_dimensions_range (foo, foo, foo, foo, max_z_, min_z_);
      }

      /**
       * \brief  Saves facilities data to facility storage
       * */
      inline void
      save_data ()
      {
        reservoir_->save_data (facility_storage_);
      }

      /**
       * \brief Checks limits
       * \return True for stop simulation and false in other cases
       * */
      inline bool
      check_limits ()
      {
        return reservoir_->check_limits (params_);
      }

      /**
       * \brief  Decreases time-step if newton iteration failed
       * \param  old_ts
       * \param  max_ts
       * \return New time-step value
       * */
      inline item_t
      decrease_ts (double old_ts, double max_ts) const
        {
          double new_ts = old_ts * params_->get_float (fi_params::TS_DEC_MULT);
          if (new_ts < params_->get_float (fi_params::MIN_TS))
            new_ts = params_->get_float (fi_params::MIN_TS);
          if (new_ts > params_->get_float (fi_params::MAX_TS))
            new_ts = params_->get_float (fi_params::MAX_TS);

          if (new_ts * params_->get_float (fi_params::OVERDRAFT) > max_ts)
            new_ts = max_ts;

          return new_ts;
        }

      /**
       * \brief  Increases time-step if previous newton iteration successed
       * \param  old_ts
       * \param  max_ts
       * \param  n_iters
       * \return New time-step value
       * */
      inline item_t
      increase_ts (item_t old_ts, item_t max_ts, index_t n_iters) const
        {
          item_t new_ts;
          index_t i, n;
          item_t dp = 0;
          item_t ds = 0;
          item_t drs = 0;
          item_t dp_max = 0;
          item_t ds_max = 0;
          item_t drs_max = 0;
          item_t user_dp;
          item_t user_ds;
          item_t user_drs;
          item_t user_omega;
          item_t alpha;
          item_t alpha2;
          //index_t ds_w, ds_g;
          //index_t flag = 0;
          index_t dp_i = 0, ds_i = 0, drs_i = 0;
          index_t n_left = 0, n_right = 0;

#ifdef _MPI
          item_t mpi_dp_max, mpi_ds_max, mpi_drs_max;
#endif //_MPI

          if (!params_->get_bool (fi_params::NEW_TS_SELECTION_ALGORITHM))
            {
              if (n_iters <= params_->get_int (fi_params::NEWTON_ITERS_NUM)
                  - params_->get_int (fi_params::NEWTON_ITERS_INC))
                new_ts = old_ts * params_->get_float (fi_params::TS_INC_MULT);
              else
                new_ts = old_ts;
            }
          else
            {
              BS_ASSERT (mesh_);

              user_dp = params_->get_float (fi_params::TS_DP);
              user_ds = params_->get_float (fi_params::TS_DS);
              user_drs = params_->get_float (fi_params::TS_DRS);
              user_omega = params_->get_float (fi_params::TS_OMEGA);
              //ds_w = calc_model_->sat_d[FI_PHASE_WATER];
              //ds_g = calc_model_->sat_d[FI_PHASE_GAS];

              n = (int)mesh_->get_n_active_elements();

#ifdef _MPI
              n_left = 0;///mpi_decomp->get_recv_l ();
              n_right = n;///mpi_decomp->get_n_local_own () + n_left;
#else //_MPI
              n_left = 0;
              n_right = n;
#endif //_MPI

              // calculate pressure and saturation changing
              for (i = n_left; i < n_right; ++i)
                {
                  dp = ((*calc_model_->pressure)[i] - (*calc_model_->old_data_.pressure)[i]);
                  if (fabs (dp_max) < fabs(dp))
                    {
                      dp_max = dp;
                      dp_i = i;
                    }
                }

              n_left *= (calc_model_->n_phases);
              n_right *= (calc_model_->n_phases);

              if (calc_model_->n_phases > 1)
                for (i = n_left; i < n_right; ++i)
                  {
                    ds = ((*calc_model_->saturation_3p)[i] - (*calc_model_->old_data_.saturation_3p)[i]);
                    if (fabs (ds_max) < fabs(ds))
                      {
                        ds_max = ds;
                        ds_i = (int)(i / (calc_model_->n_phases));
                      }
                  }

              if (FI_CHK_OIL_GAS (calc_model_->phases))
                {
#ifdef _MPI
                  n_left = 0;///mpi_decomp->get_recv_l ();
                  n_right = n;///mpi_decomp->get_n_local_own () + n_left;
#else //_MPI
                  n_left = 0;
                  n_right = n;
#endif //_MPI
                  for (i = n_left; i < n_right; ++i)
                    {
                      drs = ((*calc_model_->old_data_.gas_oil_ratio)[i] - (*calc_model_->gas_oil_ratio)[i]);
                      if (fabs(drs_max) < fabs(drs))
                        {
                          drs_max = drs;
                          drs_i = i;
                        }

                    }
                }

#ifdef _MPI
              BS_ASSERT(false&&"MPI: NOT IMPL YET");
              //dp_max = fabs (dp_max);
              //ds_max = fabs (ds_max);
              //drs_max = fabs (drs_max);
              //MPI_Allreduce (&dp_max, &mpi_dp_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
              //MPI_Allreduce (&ds_max, &mpi_ds_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
              //MPI_Allreduce (&drs_max, &mpi_drs_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
              //dp_max = mpi_dp_max;
              //ds_max = mpi_ds_max;
              //drs_max = mpi_drs_max;
#endif //_MPI

              dp = (1.0 + user_omega) * user_dp / (dp_max + user_omega * user_dp);
              if (calc_model_->n_phases > 1)
                ds = (1.0 + user_omega) * user_ds / (ds_max + user_omega * user_ds);
              drs = (1.0 + user_omega) * user_drs / (drs_max + user_omega * user_drs);

              BOSOUT (section::main_loop, level::low)
                << "dp_max[ " << dp_i << "] = " << dp_max
                << ", ds_max[" << ds_i << "] = " << ds_max
                << ", drs_max[" << drs_i << "] = " << drs_max
                << bs_end;

              alpha = dp;

              if (calc_model_->n_phases > 1)
                if (ds < alpha)
                  {
                    alpha = ds;
                    //flag = 1;
                  }
              if (drs < alpha)
                {
                  alpha = drs;
                  //flag = 2;
                }
              //printf ("FLAG = %d\n", flag);
              if (alpha < 1)
                alpha = 1;
              alpha2 = 1.0;
              if (n_iters <= params_->get_int (fi_params::NEWTON_ITERS_NUM) - params_->get_int (fi_params::NEWTON_ITERS_INC))
                alpha2 = params_->get_float (fi_params::TS_INC_MULT);
              if (alpha > alpha2)
                alpha = alpha2;
              new_ts = old_ts * alpha;
            }

          if (new_ts < params_->get_float (fi_params::MIN_TS))
            new_ts = params_->get_float (fi_params::MIN_TS);
          if (new_ts > params_->get_float (fi_params::MAX_TS))
            new_ts = params_->get_float (fi_params::MAX_TS);

          // select time step length
          if (new_ts * params_->get_float (fi_params::OVERDRAFT) > max_ts)
            new_ts = max_ts;
          else if (new_ts > max_ts)
            new_ts = max_ts;
          else if (max_ts - new_ts < new_ts)
            new_ts = max_ts * 0.5;

          BOSOUT (section::main_loop, level::low) << "increase ts: " << old_ts << " -> " << new_ts << bs_end;
          return new_ts;
        }

      /**
       * \brief  Saves base vars
       * \todo   Bad design
       * */
      inline void
      save_base_vars ()
      {
        calc_model_->old_data_.save (calc_model_);
        reservoir_->pre_small_step ();
      }
      /**
       * \brief  Restores base vars
       * \todo   Bad design
       * */
      inline void
      restore_base_vars ()
      {
#ifdef _DEBUG
        BOSOUT (section::main_loop, level::debug) << "restore_base_vars" << bs_end;
#endif
        calc_model_->old_data_.restore (calc_model_);
        reservoir_->restart_small_step ();
      }

      /**
       * \brief  Updates total number of newton iterations
       * \param  nniters
       * */
      inline void
      update_number_of_newtonian_iterations (index_t nniters)
      {
        number_of_newtonian_iterations += nniters;
      }
      /**
       * \brief  Updates total number of linear iterations
       * \param  nliters
       * */
      inline void
      update_number_of_linear_iterations (index_t nliters)
      {
        number_of_linear_iterations += nliters;
      }
      /**
       * \brief  Updates current time-step
       * \param  step_
       * */
      inline void
      update_current_time (item_t step_)
      {
        current_time_ += step_;
      }
      /**
       * \brief  Updates total number of fi_operator restarts
       * */
      inline void
      update_number_of_fi_operator_restarts ()
      {
        ++number_of_fi_operator_restarts;
      }
      /**
       * \brief  Updates total number of restarts
       * */
      inline void
      update_number_of_restarts ()
      {
        ++number_of_restarts;
      }
      inline void
      update_number_of_max_iters_restarts ()
      {
        ++number_of_max_iters_restarts;
      }
      inline void
      update_number_of_close_wells_restarts ()
      {
        ++number_of_close_wells_restarts;
      }


      /**
       * \brief  Updates length of large time-step (from model)
       * \param  current_time
       * \param  next_time
       * \return True if length of time-step positive
       * */
      BS_FORCE_INLINE bool
      update_large_time_step_length (double current_time, double next_time)
      {
        //using namespace boost::posix_time;

        //time_period tp (current_time, next_time);
        //time_duration td (tp.length ());

        //large_time_step_length_ = item_t (td.total_milliseconds ()) / item_t (24.0 * 60 * 60 * 1000);
        large_time_step_length_ = next_time - current_time;

        if (large_time_step_length_ > 0)
          {
            large_time_step_start_ = large_time_step_end_;
            large_time_step_end_ = large_time_step_end_ + large_time_step_length_;
            return true;
          }

        return false;
      }

      /**
       * \brief  Updates number of processed large time-steps
       * */
      inline void
      update_large_time_step_num ()
      {
        ++number_of_large_time_steps;
      }

      /**
       * \brief  Prints mesh info
       * */
      void
      print_mesh_info ()
      {
        BOSOUT (section::mesh, level::medium) << "mesh cells: " << mesh_->get_n_active_elements () << bs_end;
        BOSOUT (section::mesh, level::medium) << "mesh connections: " << mesh_->get_n_connections () << bs_end;
      }

      /**
       * \brief  Prints PVT tables
       * */
      void
      print_pvt_info ()
      {
        for (size_t i = 0, cnt = calc_model_->pvt_water_array.size (); i < cnt; ++i)
          {
            calc_model_->pvt_water_array[i]->print ();
          }
        for (size_t i = 0, cnt = calc_model_->pvt_oil_array.size (); i < cnt; ++i)
          {
            calc_model_->pvt_oil_array[i]->print ();
          }
        for (size_t i = 0, cnt = calc_model_->pvt_gas_array.size (); i < cnt; ++i)
          {
            calc_model_->pvt_gas_array[i]->print ();
          }
      }

      /**
       * \brief  Large time-step iteration
       * \param  current_time
       * \param  next_time
       * \param  event_list
       * */
      inline void
      iteration (double current_time, double next_time, const event_list_t &event_list)
      {
        using namespace boost::posix_time;
        BOSOUT (section::main_loop, level::high) << "\nTIMESTEP " /*<< boost::posix_time::to_simple_string (current_time) */<< "\n" << bs_end;

        if (update_large_time_step_length (current_time, next_time))
          {
            rs_->pre_large_step (event_list);
            //trans_multipliers_.apply ();
            rs_->on_large_step_start ();

            process_small_steps ();
            well_data_printer::print_prod (data_map_, reservoir_);
            well_data_printer::print_inj  (data_map_, reservoir_);
            well_data_printer::print_total_prod  (data_map_, reservoir_);

            time_step_end ((int)event_list.size ());
          }
      }

      /**
       * \brief  Performs actions on time-step end
       * \param  total_number_of_time_steps
       * */
      inline void
      time_step_end (int total_number_of_time_steps)
      {
        save_data ();
        check_limits ();
        update_large_time_step_num ();

        reservoir_->write_step_to_storage (calc_model_, mesh_, jacobian_, number_of_large_time_steps, total_number_of_time_steps, current_time_);

        BOSOUT (section::main_loop, level::high) << "number_of_large_time_steps: "      << number_of_large_time_steps << bs_end;
        BOSOUT (section::main_loop, level::high) << "number_of_small_time_steps: "      << number_of_small_time_steps << bs_end;
        BOSOUT (section::main_loop, level::high) << "number_of_newtonian_iterations: "  << number_of_newtonian_iterations << bs_end;
        BOSOUT (section::main_loop, level::high) << "number_of_linear_iterations: "     << number_of_linear_iterations << bs_end;
        BOSOUT (section::main_loop, level::high) << "number_of_restarts: "              << number_of_restarts << bs_end;
        BOSOUT (section::main_loop, level::high) << "number_of_max_iters_restarts: "    << number_of_max_iters_restarts << bs_end;
        BOSOUT (section::main_loop, level::high) << "number_of_fi_operator_restarts: "  << number_of_fi_operator_restarts << bs_end;
        BOSOUT (section::main_loop, level::high) << "number_of_close_wells_restarts: "  << number_of_close_wells_restarts << bs_end;
      }

      /**
       * \brief  Prepares calculation
       * */
      inline void
      ready ()
      {
        print_mesh_info ();
        print_pvt_info ();

        fi_operator_cells (0);
        trans_multipliers_.apply ();

        dt_ = get_first_ts ();
        save_base_vars ();

//#ifdef _HDF5
//        reservoir_->open_hdf5_file (path::join (path::dirname (rs_->model_filename ()), "results.h5"));
//        reservoir_->write_mesh_to_hdf5 (mesh_);
//        boost::gregorian::date start_date = rs_->get_hdm ()->get_keyword_manager ()->get_starting_date ().date ();
//        boost::gregorian::date base_date (1900, 1, 1);
//        double starting_date = (start_date - base_date).days () + 2;
//        reservoir_->get_hdf5_file ()->write_array ("/initial_data", "starting_date", &starting_date, 1);
//#endif

        reservoir_->open_storage (path::join (path::dirname (rs_->model_filename ()), "results-v2.h5"));
        reservoir_->write_mesh_to_storage (mesh_);
        reservoir_->write_starting_date_to_storage (rs_->get_hdm ()->get_keyword_manager ()->get_starting_date ());
      }

      /**
       * \brief  Starts calculation loop
       * */
      inline void
      go ()
      {
        event_manager::event_map::iterator it = event_manager_->event_list.begin ();
        event_manager::event_map::iterator e  = event_manager_->event_list.end ();

        rs_->on_simulation_start ();
        for (--e; it != e; ++it)
          {
            event_manager::event_map::iterator it2 = it;
            ++it2;

            iteration (it->first, it2->first, it->second);
          }
      }

      /**
       * \brief  Ends calculation
       * */
      inline void
      end ()
      {
#ifdef _HDF5
        reservoir_->close_hdf5_file ();
#endif

        BOSOUT (section::main_loop, level::high) << "TOTAL: " << bs_end;
        BOSOUT (section::main_loop, level::high) << "number_of_large_time_steps: "      << number_of_large_time_steps << bs_end;
        BOSOUT (section::main_loop, level::high) << "number_of_small_time_steps: "      << number_of_small_time_steps << bs_end;

        BOSOUT (section::main_loop, level::high) << "number_of_newtonian_iterations: "  << number_of_newtonian_iterations << bs_end;
        BOSOUT (section::main_loop, level::high) << "number_of_linear_iterations: "     << number_of_linear_iterations << bs_end;
        BOSOUT (section::main_loop, level::high) << "number_of_restarts: "              << number_of_restarts << bs_end;
        BOSOUT (section::main_loop, level::high) << "number_of_max_iters_restarts: "    << number_of_max_iters_restarts << bs_end;
        BOSOUT (section::main_loop, level::high) << "number_of_fi_operator_restarts: "  << number_of_fi_operator_restarts << bs_end;
        BOSOUT (section::main_loop, level::high) << "number_of_close_wells_restarts: "  << number_of_close_wells_restarts << bs_end;
      }

public:
      // don't change line order. never.
      sp_rs_t                                 rs_;
      sp_calc_model_t                         calc_model_;
      sp_event_manager_t                      event_manager_;
      sp_rock_grid_t                          rock_grid_prop_;
      sp_storage_t                            facility_storage_;
      sp_reservoir_t                          reservoir_;
      sp_mesh_iface_t                         mesh_;
      sp_jacobian_t                           jacobian_;
      sp_fi_params_t                          params_;
      sp_idata_t                              data_map_;
      trans_multipliers_calc_t                trans_multipliers_;

public:
      item_t                                  height_;
      item_t                                  rho_;
      item_t                                  min_z_;
      item_t                                  max_z_;

      double                                  dt_;
      double                                  current_time_;

      item_t                                  large_time_step_length_;
      item_t                                  large_time_step_start_;
      item_t                                  large_time_step_end_;

      index_t                                 num_last_newton_iters_;
      index_t                                 num_last_lin_iters_;

      index_t                                 number_of_small_time_steps;
      index_t                                 number_of_large_time_steps;

      index_t                                 number_of_newtonian_iterations;
      index_t                                 number_of_linear_iterations;
      index_t                                 number_of_restarts;
      index_t                                 number_of_max_iters_restarts;
      index_t                                 number_of_fi_operator_restarts;
      index_t                                 number_of_close_wells_restarts;

      bool                                    do_calc_prev_fluid_volume_;

      //save_well_data                          well_data;
    };


} // namespace blue_sky

#include "generate_numeric_jacobian.h"
#include "compute_acc_rates.h"

#endif  // #ifndef BS_MAIN_LOOP_CALC_H_

