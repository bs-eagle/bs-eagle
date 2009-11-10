/**
 * \file save_connection_data.h
 * \brief save connection data (for debug)
 * \author Sergey Miryanov
 * \date 15.12.2008
 * */
#ifndef BS_WELLS_SAVE_CONNECTION_DATA_H_
#define BS_WELLS_SAVE_CONNECTION_DATA_H_

#include "well_rate_control.h"
#include "apply_wefac.h"

namespace blue_sky
  {

  template <typename strategy_t>
  struct save_connection_data
    {
      typedef calc_model <strategy_t>                                                 calc_model_t;
      typedef typename strategy_t::index_t                                            index_t;

      typedef typename calc_model_t::data_array_t                                     data_array_t;
      typedef typename calc_model_t::item_t                                           item_t;
      typedef typename calc_model_t::item_array_t                                     item_array_t;
      typedef typename strategy_t::rhs_item_array_t                                   rhs_item_array_t;
      typedef typename calc_model_t::helper_t::item_rr_block_t                        item_rr_block_t;
      typedef typename calc_model_t::helper_t::item_rw_block_t                        item_rw_block_t;
      typedef typename calc_model_t::helper_t::item_wr_block_t                        item_wr_block_t;
      typedef typename calc_model_t::connection_t                                     connection_t;
      typedef typename calc_model_t::well_t                                           well_t;
      typedef typename calc_model_t::reservoir_t::facility_manager_t::well_const_iterator_t well_iterator_t;
      typedef typename calc_model_t::strategy_type                                    strategy_type;

      typedef smart_ptr <calc_model_t, true>                                          sp_calc_model_t;
      typedef typename calc_model_t::sp_well_t                                        sp_well_t;
      typedef typename calc_model_t::sp_connection_t                                  sp_connection_t;

      void
      save (const char *filename, const sp_calc_model_t &calc_model, item_t dt, well_iterator_t wb, const well_iterator_t &we, size_t iter_counter, item_t time, const item_array_t &sol, const item_array_t &rhs)
      {
        FILE *con_data = fopen (tools::string_formater (filename, iter_counter).str, "wt");
        index_t n_phases = calc_model->n_phases;
        bool is_o = calc_model->is_oil ();
        bool is_w = calc_model->is_water ();
        bool is_g = calc_model->is_gas ();
        for (well_iterator_t it = wb; it != we; ++it)
          {
            sp_well_t well (it->second, bs_dynamic_cast ());
            item_t h = well->get_connection_list ().front ()->connection_depth;
            item_t diff_h = 0;
            item_t g = calc_model->internal_constants.gravity_constant;
            for (size_t i = 0, cnt = well->get_connections_count (); i < cnt; ++i)
              {
                const sp_connection_t &c = (well->get_connection_list ())[i];

                index_t n_block = c->n_block ();
                fprintf (con_data, "%d\n", n_block);
                fprintf (con_data, "%10.20lf\n", apply_wefac (c->get_rate_value () [0], well->exploitation_factor_));
                fprintf (con_data, "%10.20lf\n", apply_wefac (c->get_rate_value () [1], well->exploitation_factor_));
                fprintf (con_data, "%10.20lf\n", apply_wefac (c->get_rate_value () [2], well->exploitation_factor_));
                fprintf (con_data, "%10.20lf\n", c->get_cur_bhp ());
                fprintf (con_data, "%10.20lf\n", calc_model->pressure [n_block]);
                if (calc_model->n_phases > 1)
                  {
                    fprintf (con_data, "%10.20lf\n", calc_model->saturation_3p [n_block]);
                  }
                else
                  {
                    fprintf (con_data, "0\n");
                  }
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].mobility [0]);
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].mobility [1]);
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].mobility [2]);
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].relative_perm [0]);
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].relative_perm [1]);
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].relative_perm [2]);
                fprintf (con_data, "%10.20lf\n", well->bhp ());
                fprintf (con_data, "%10.20lf\n", well->oil_rate ());
                fprintf (con_data, "%10.20lf\n", well->water_rate ());
                fprintf (con_data, "%10.20lf\n", well->gas_rate ());
                fprintf (con_data, "%10.20lf\n", well->get_well_controller ()->bhp ());
                fprintf (con_data, "%10.20lf\n", well->get_well_controller ()->oil_rate ());
                fprintf (con_data, "%10.20lf\n", well->get_well_controller ()->water_rate ());
                fprintf (con_data, "%10.20lf\n", well->get_well_controller ()->gas_rate ());
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].invers_fvf [0]);
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].invers_fvf [1]);
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].invers_fvf [2]);
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].prev_fluid_volume [0]);
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].prev_fluid_volume [1]);
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].prev_fluid_volume [2]);
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].porosity);
                fprintf (con_data, "%10.20lf\n", calc_model->rock_grid_prop->volume [n_block]);

                fprintf (con_data, "%10.20lf\n", well->acc_rate_prod.oil);
                fprintf (con_data, "%10.20lf\n", well->acc_rate_prod.water);
                fprintf (con_data, "%10.20lf\n", well->acc_rate_prod.gas);
                fprintf (con_data, "%10.20lf\n", well->acc_rate_inj.oil);
                fprintf (con_data, "%10.20lf\n", well->acc_rate_inj.water);
                fprintf (con_data, "%10.20lf\n", well->acc_rate_inj.gas);

                fprintf (con_data, "%10.20lf\n", c->get_fact ());

                using namespace wells;
                if (well->is_shut ())
                  {
                    fprintf (con_data, "%d\n", 0);
                  }
                else if (well->get_well_controller ()->is_production ())
                  {
                    if (well->is_bhp ())
                      fprintf (con_data, "%d\n", -2);
                    else if (well->get_well_controller ()->get_control_type () == oil_rate_control)
                      fprintf (con_data, "%d\n", -4);
                    else if (well->get_well_controller ()->get_control_type () == water_rate_control)
                      fprintf (con_data, "%d\n", -3);
                    else if (well->get_well_controller ()->get_control_type () == liquid_rate_control)
                      fprintf (con_data, "%d\n", -1);
                    else
                      fprintf (con_data, "%d\n", -5);
                  }
                else
                  {
                    if (well->is_bhp ())
                      fprintf (con_data, "%d\n", 2);
                    else if (well->get_well_controller ()->injection () == injection_water)
                      fprintf (con_data, "%d\n", 1);
                    else
                      fprintf (con_data, "%d\n", 5);
                  }

                fprintf (con_data, "%10.20lf\n", well->exploitation_factor_.data ());
                fprintf (con_data, "%10.20lf\n", dt);
                fprintf (con_data, "%10.20lf\n", time);
                fprintf (con_data, "%10.20lf\n", well->bw_value);
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].invers_viscosity [0]);
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].invers_viscosity [1]);
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].invers_viscosity [2]);
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].invers_visc_fvf [0]);
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].invers_visc_fvf [1]);
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].invers_visc_fvf [2]);

                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].p_deriv_invers_fvf [0]);
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].p_deriv_invers_fvf [1]);
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].p_deriv_invers_fvf [2]);
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].p_deriv_invers_viscosity [0]);
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].p_deriv_invers_viscosity [1]);
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].p_deriv_invers_viscosity [2]);
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].p_deriv_invers_visc_fvf [0]);
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].p_deriv_invers_visc_fvf [1]);
                fprintf (con_data, "%10.20lf\n", calc_model->data [n_block].p_deriv_invers_visc_fvf [2]);

                fprintf (con_data, "%10.20lf\n", well->rate_prod.oil);
                fprintf (con_data, "%10.20lf\n", well->rate_prod.water);
                fprintf (con_data, "%10.20lf\n", well->rate_prod.gas);
                fprintf (con_data, "%10.20lf\n", well->rate_inj.oil);
                fprintf (con_data, "%10.20lf\n", well->rate_inj.water);
                fprintf (con_data, "%10.20lf\n", well->rate_inj.gas);

                diff_h    = c->connection_depth - h;
                h         = c->connection_depth;
                item_t Po = c->cur_bhp + (c->density * g * diff_h) - calc_model->pressure[c->n_block ()];
                item_t Pw = Po;
                item_t Pg = Po;

                if (calc_model->is_water ())
                  Pw += calc_model->data[n_block].cap_pressure[0];

                if (calc_model->is_gas ())
                  Pg += calc_model->data[n_block].cap_pressure[1];

                fprintf (con_data, "%10.20lf\n", Po);
                fprintf (con_data, "%10.20lf\n", Pw);
                fprintf (con_data, "%10.20lf\n", Pg);

                fprintf (con_data, "%10.20lf\n", (double)c->density);
                fprintf (con_data, "%10.20lf\n", (double)c->head_term);
                fprintf (con_data, "%10.20lf\n", (double)c->connection_depth);

                fprintf (con_data, "%10.20lf\n", sol[n_block * n_phases]);
                fprintf (con_data, "%10.20lf\n", sol[n_block * n_phases + 1]);

                fprintf (con_data, "%10.20lf\n", rhs[n_block * n_phases]);
                fprintf (con_data, "%10.20lf\n", rhs[n_block * n_phases + 1]);
              }
          }
        fclose (con_data);
      }

      void
      save_acc (const char *filename, const sp_calc_model_t &/*calc_model*/, item_t /*dt*/, well_iterator_t wb, const well_iterator_t &we, size_t iter_counter, item_t /*time*/, const item_array_t &/*sol*/, const rhs_item_array_t &/*rhs*/)
      {
        FILE *con_data = fopen (tools::string_formater (filename, iter_counter).str, "wt");
        for (well_iterator_t it = wb; it != we; ++it)
          {
            const sp_well_t &well = it->second;
            fprintf (con_data, "well: %s\n", well->name ().c_str ());
            fprintf (con_data, "prod_oil:   %10.20lf\n", well->acc_rate_prod.oil);
            fprintf (con_data, "prod_water: %10.20lf\n", well->acc_rate_prod.water);
            fprintf (con_data, "prod_gas:   %10.20lf\n", well->acc_rate_prod.gas);
            fprintf (con_data, "inj_oil:    %10.20lf\n", well->acc_rate_inj.oil);
            fprintf (con_data, "inj_water:  %10.20lf\n", well->acc_rate_inj.water);
            fprintf (con_data, "inj_gas:    %10.20lf\n", well->acc_rate_inj.gas);

            //for (size_t i = 0, cnt = well->get_connections_count (); i < cnt; ++i)
            //  {
            //    const sp_connection_t &c = (well->get_connection_list ())[i];
            //    index_t n_block = c->n_block ();
            //    fprintf (con_data, "n_block: %d\n", n_block);
            //    fprintf (con_data, "cur_bhp: %10.20lf\n", c->get_cur_bhp ());
            //  }
          }

        fclose (con_data);
      }
    };


} // namespace blue_sky



#endif  // #ifndef BS_WELLS_SAVE_CONNECTION_DATA_H_

