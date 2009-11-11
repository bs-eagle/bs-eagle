/**
 * \file well_results_storage.cpp
 * \brief Storage for well and connection results
 * \author Ilshat Sayfullin
 * \date 25.09.2009
 * */
#ifndef __BS_WELL_RESULTS_STORAGE_H
#define __BS_WELL_RESULTS_STORAGE_H

#include "well_rate_control.h"
#include "apply_wefac.h"
#include "rate_control_type.h"
#include "fi_params.h"
//#include "well_connection.h"

namespace blue_sky
  {

  enum well_d_params
  {
    WELL_D_PARAM_COR = 0,         //0 - current oil rate
    WELL_D_PARAM_CWR,             //1 - current water rate
    WELL_D_PARAM_CGR,             //2 - current gas rate
    WELL_D_PARAM_CLR,             //3 - current liquid rate
    WELL_D_PARAM_COI,             //4 - current oil injection
    WELL_D_PARAM_CWI,             //5 - current water injection
    WELL_D_PARAM_CGI,             //6 - current gas injection

    WELL_D_PARAM_HCOR,            //7 - current history oil rate
    WELL_D_PARAM_HCWR,            //8 - current history water rate
    WELL_D_PARAM_HCGR,            //9 - current history gas rate
    WELL_D_PARAM_HCLR,            //10- current history liquid rate
    WELL_D_PARAM_HCOI,            //11- current history oil injection
    WELL_D_PARAM_HCWI,            //12- current history water injection
    WELL_D_PARAM_HCGI,            //13- current history gas injection

    WELL_D_PARAM_CBHP,            //14- BHP for connection
    WELL_D_PARAM_HBHP,            //15- BHP defined for well

    WELL_D_PARAM_TOR,             //16- total oil rate
    WELL_D_PARAM_TWR,             //17- total water rate
    WELL_D_PARAM_TGR,             //18- total gas rate
    WELL_D_PARAM_TLR,             //19- total liquid rate
    WELL_D_PARAM_TOI,             //20- total oil injection
    WELL_D_PARAM_TWI,             //21- total water injection
    WELL_D_PARAM_TGI,             //22- total gas injection

    WELL_D_PARAM_HTOR,            //23- total history oil rate
    WELL_D_PARAM_HTWR,            //24- total history water rate
    WELL_D_PARAM_HTGR,            //25- total history gas rate
    WELL_D_PARAM_HTLR,            //26- total history liquid rate
    WELL_D_PARAM_HTOI,            //27- total history oil injection
    WELL_D_PARAM_HTWI,            //28- total history water injection
    WELL_D_PARAM_HTGI,            //29- total history gas injection
    WELL_D_PARAM_WEFAC,           //30- wefac (effectivity of well)
    WELL_D_PARAM_WBP1,            //31- average well's pressure for 1 block
    WELL_D_PARAM_WBP4,            //32- average well's pressure for 4 block
    WELL_D_PARAM_WBP5,            //33- average well's pressure for 5 block
    WELL_D_PARAM_WBP9,            //34- average well's pressure for 9 block

    // should be the last
    WELL_D_PARAM_COUNT
  };

  enum conn_d_params
  {
    CONN_D_PARAM_COR = 0,         //0 - current oil rate
    CONN_D_PARAM_CWR,             //1 - current water rate
    CONN_D_PARAM_CGR,             //2 - current gas rate
    CONN_D_PARAM_CLR,             //3 - current liquid rate
    CONN_D_PARAM_COI,             //4 - current oil injection
    CONN_D_PARAM_CWI,             //5 - current water injection
    CONN_D_PARAM_CGI,             //6 - current gas injection
    CONN_D_PARAM_BHP,             //7 - BHP for connection
    CONN_D_PARAM_BULKP,           //8 - bulk pressure
    CONN_D_PARAM_TOR,             //9 - total oil rate
    CONN_D_PARAM_TWR,             //10- total water rate
    CONN_D_PARAM_TGR,             //11- total gas rate
    CONN_D_PARAM_TLR,             //12- total liquid rate
    CONN_D_PARAM_TOI,             //13- total oil injection
    CONN_D_PARAM_TWI,             //14- total water injection
    CONN_D_PARAM_TGI,             //15- total gas injection
    CONN_D_PARAM_FACTOR,          //16- connection factor
    CONN_D_PARAM_DEPTH,           //17- connection depth

    // should be the last
    CONN_D_PARAM_COUNT
  };

  enum conn_i_params
  {
    CONN_I_PARAM_STATUS = 0,      //0 - state of connection: 0 - close, 1 - open
    CONN_I_PARAM_GRP_STATUS,      //1 - type of connection: 0 - usual, 1 - grp_primary, 2 - grp_secondary
    CONN_I_PARAM_I,               //2 - index of connection block by X axes in mesh
    CONN_I_PARAM_J,               //3 - index of connection block by Y axes in mesh
    CONN_I_PARAM_K1,              //4 - start layer index of connection block by Z
    CONN_I_PARAM_K2,              //5 - end layer index of connection block by Z
    CONN_I_PARAM_NBLOCK,          //6 - block number in mesh

    // should be the last
    CONN_I_PARAM_COUNT
  };

  enum well_i_params
  {
    WELL_I_PARAM_HSTATUS = 0,
    WELL_I_PARAM_STATUS,

    // should be the last
    WELL_I_PARAM_COUNT
  };

  class BS_API_PLUGIN connection_results
    {
      // methods
    public:
      connection_results ();
      ~connection_results ();

      void
      clear ();

      // variables
    public:
      typedef std::vector<double> dates_type;
      typedef std::vector<float> d_params_internal_type;
      typedef std::vector<d_params_internal_type> d_params_type;
      typedef std::vector<int> i_params_internal_type;
      typedef std::vector<i_params_internal_type> i_params_type;

      dates_type dates;                           //!< dates array
      d_params_type d_params;                     //!< set of float parameter arrays
      i_params_type i_params;                     //!< set of int   parameter arrays
    };

  class BS_API_PLUGIN well_results
    {
      // methods
    public:
      well_results ();
      ~well_results ();

      void
      clear ();

      // variables
    public:
      typedef std::vector<double> dates_type;
      typedef std::vector<float> d_params_internal_type;
      typedef std::vector<d_params_internal_type> d_params_type;
      typedef std::vector<int> i_params_internal_type;
      typedef std::vector<i_params_internal_type> i_params_type;
      typedef std::map<int, connection_results> conn_type;
      typedef std::pair<int, connection_results> conn_pair_type;

      conn_type connections;
      dates_type dates;
      d_params_type d_params;
      i_params_type i_params;
      std::string group;
    };

  class BS_API_PLUGIN well_results_storage : public objbase
    {
      BLUE_SKY_TYPE_DECL (well_results_storage);
      // methods
    public:
      ~well_results_storage ();

      // variables
    public:
      typedef std::map<std::string, well_results> wells_type;
      typedef std::pair<std::string, well_results> wells_pair_type;
      wells_type wells;
    };


  template <typename strategy_t>
  struct save_well_data
    {
      typedef calc_model <strategy_t>                                                 calc_model_t;
      typedef typename strategy_t::index_t                                            index_t;

      typedef typename calc_model_t::data_array_t                                     data_array_t;
      typedef typename calc_model_t::item_t                                           item_t;
      typedef typename calc_model_t::connection_t                                     connection_t;
      typedef typename calc_model_t::well_t                                           well_t;
      typedef typename calc_model_t::reservoir_t::facility_manager_t::well_const_iterator_t well_iterator_t;
      typedef typename calc_model_t::strategy_type                                    strategy_type;

      typedef smart_ptr <calc_model_t, true>                                          sp_calc_model_t;
      typedef typename calc_model_t::sp_well_t                                        sp_well_t;
      typedef typename calc_model_t::sp_connection_t                                sp_connection_t;
      typedef smart_ptr<well_results_storage, true>     sp_well_results_storage;

      void
      copy_well_data_to_storage (sp_calc_model_t &calc_model, item_t /*dt*/, well_iterator_t wb, const well_iterator_t &we, size_t /*iter_counter*/, item_t time)
      {
        const sp_well_results_storage &w_res (calc_model->well_res);

        //TODO: groups
        for (well_iterator_t well = wb; well != we; ++well)
          {
            sp_well_t ws (well->second, bs_dynamic_cast ());
            well_results &wr = w_res->wells[ws->get_name ()];

            wr.group = std::string ("FIELD");//TODO ws->get_group_name ()
            // add time value
            wr.dates.push_back (time);

            // add current rates
            wr.d_params[WELL_D_PARAM_COR].push_back (-(float)ws->rate ().prod.oil);
            wr.d_params[WELL_D_PARAM_CWR].push_back (-(float)ws->rate ().prod.water);
            wr.d_params[WELL_D_PARAM_CGR].push_back (-(float)ws->rate ().prod.gas);
            wr.d_params[WELL_D_PARAM_CLR].push_back (-(float)ws->rate ().prod.liquid);
            wr.d_params[WELL_D_PARAM_COI].push_back ((float)ws->rate ().inj.oil);
            wr.d_params[WELL_D_PARAM_CWI].push_back ((float)ws->rate ().inj.water);
            wr.d_params[WELL_D_PARAM_CGI].push_back ((float)ws->rate ().inj.gas);


            // add initial rates
            wr.d_params[WELL_D_PARAM_HCOR].push_back (0);
            wr.d_params[WELL_D_PARAM_HCWR].push_back (0);
            wr.d_params[WELL_D_PARAM_HCGR].push_back (0);
            wr.d_params[WELL_D_PARAM_HCLR].push_back (0);
            wr.d_params[WELL_D_PARAM_HCOI].push_back (0);
            wr.d_params[WELL_D_PARAM_HCWI].push_back (0);
            wr.d_params[WELL_D_PARAM_HCGI].push_back (0);

            wr.d_params[WELL_D_PARAM_TOR].push_back (-(float)ws->rate_total ().prod.oil);
            wr.d_params[WELL_D_PARAM_TWR].push_back (-(float)ws->rate_total ().prod.water);
            wr.d_params[WELL_D_PARAM_TGR].push_back (-(float)ws->rate_total ().prod.gas);
            wr.d_params[WELL_D_PARAM_TLR].push_back (-(float)ws->rate_total ().prod.liquid);
            wr.d_params[WELL_D_PARAM_TOI].push_back ((float)ws->rate_total ().inj.oil);
            wr.d_params[WELL_D_PARAM_TWI].push_back ((float)ws->rate_total ().inj.water);
            wr.d_params[WELL_D_PARAM_TGI].push_back ((float)ws->rate_total ().inj.gas);

            /*          //TODO
                        wr.d_params[WELL_D_PARAM_HTOR].push_back (-(float)ws->total_initial_rate.oil);
                        wr.d_params[WELL_D_PARAM_HTWR].push_back (-(float)ws->total_initial_rate.water);
                        wr.d_params[WELL_D_PARAM_HTGR].push_back (-(float)ws->total_initial_rate.gas);
                        wr.d_params[WELL_D_PARAM_HTLR].push_back ((float)ws->total_initial_rate.liquid);
                        wr.d_params[WELL_D_PARAM_HTOI].push_back ((float)ws->total_initial_rate.oil_inj);
                        wr.d_params[WELL_D_PARAM_HTWI].push_back ((float)ws->total_initial_rate.water_inj);
                        wr.d_params[WELL_D_PARAM_HTGI].push_back ((float)ws->total_initial_rate.gas_inj);
                        wr.d_params[WELL_D_PARAM_WEFAC].push_back ((float)ws->wefac);
            */

            wr.d_params[WELL_D_PARAM_CBHP].push_back ((float)ws->get_well_controller ()->bhp ());
            wr.d_params[WELL_D_PARAM_HBHP].push_back ((float)ws->get_well_controller ()->bhp_history ());

            // add integer information
            using namespace wells;
            int status = -999, producer = 0;
            if (ws->get_well_controller ()->is_production ())
              producer = 1;
            rate_control_type control_type = ws->get_well_controller ()->get_control_type ();
            if (!ws->is_open ())
              {
                status = 0;
              }
            else if (control_type == bhp_control)
              {
                if (producer)
                  status = -2;
                else
                  status = 2;
              }
            else if (control_type == oil_rate_control) // (assume is producer)
              {
                status = -4;
              }
            else if (control_type == water_rate_control)
              {
                if (producer)
                  status = -3;
                else
                  status = 1;
              }
            else if (control_type == gas_rate_control) // (assume is producer)
              {
                status = -5;
              }
            else if (control_type == rate_control) // (assume is injector)
              {
                status = 3;
              }
            else if (control_type == liquid_rate_control)// (assume is producer)
              {
                status = -1;
              }

            wr.i_params[WELL_I_PARAM_HSTATUS].push_back (0);//TODO
            wr.i_params[WELL_I_PARAM_STATUS].push_back (status);

            if (calc_model->ts_params->get_bool (fi_params::WRITE_CONN_RESULTS_TO_HDF5))
              {
                for (size_t i = 0, cnt = ws->get_connections_count (); i < cnt; ++i)
                  {
                    const sp_connection_t &ci = ws->get_connection (i);
                    int cell = ci->n_block ();
                    connection_results &cr = wr.connections[cell];
                    cr.dates.push_back (time);

                    if (!ci->is_shut ())
                      {/*
                                                //TODO
                                                cr.d_params[CONN_D_PARAM_COR].push_back (-(float)ci->second.rate.water);
                                                cr.d_params[CONN_D_PARAM_CWR].push_back (-(float)ci->second.rate.water);
                                                cr.d_params[CONN_D_PARAM_CGR].push_back (-(float)ci->second.rate.gas);
                                                cr.d_params[CONN_D_PARAM_CLR].push_back (-(float)ci->second.rate.liquid);
                                                cr.d_params[CONN_D_PARAM_COI].push_back ((float)ci->second.rate.oil_inj);
                                                cr.d_params[CONN_D_PARAM_CWI].push_back ((float)ci->second.rate.water_inj);
                                                cr.d_params[CONN_D_PARAM_CGI].push_back ((float)ci->second.rate.gas_inj);
                                                */
                      }
                    else
                      {
                        cr.d_params[CONN_D_PARAM_COR].push_back (0);
                        cr.d_params[CONN_D_PARAM_CWR].push_back (0);
                        cr.d_params[CONN_D_PARAM_CGR].push_back (0);
                        cr.d_params[CONN_D_PARAM_CLR].push_back (0);
                        cr.d_params[CONN_D_PARAM_COI].push_back (0);
                        cr.d_params[CONN_D_PARAM_CWI].push_back (0);
                        cr.d_params[CONN_D_PARAM_CGI].push_back (0);
                      }

                    cr.d_params[CONN_D_PARAM_FACTOR].push_back ((float)ci->get_fact ());
                    cr.d_params[CONN_D_PARAM_BHP].push_back ((float)ci->get_cur_bhp ());
                    cr.d_params[CONN_D_PARAM_BULKP].push_back ((float)ci->get_bulkp ());
                    cr.d_params[CONN_D_PARAM_DEPTH].push_back ((float)ci->connection_depth);

                    cr.i_params[CONN_I_PARAM_STATUS].push_back (ci->get_status ());
                    cr.i_params[CONN_I_PARAM_GRP_STATUS].push_back (0);//TODO
                    cr.i_params[CONN_I_PARAM_I].push_back (ci->i_coord ());
                    cr.i_params[CONN_I_PARAM_J].push_back (ci->j_coord ());
                    cr.i_params[CONN_I_PARAM_K1].push_back (ci->k_coord ());
                    cr.i_params[CONN_I_PARAM_K2].push_back (ci->k_coord ());

                    /*                  //TODO
                                        cr.d_params[CONN_D_PARAM_TOR].push_back (-(float)ci->second.commulative_rate.oil);
                                        cr.d_params[CONN_D_PARAM_TWR].push_back (-(float)ci->second.commulative_rate.water);
                                        cr.d_params[CONN_D_PARAM_TGR].push_back (-(float)ci->second.commulative_rate.gas);
                                        cr.d_params[CONN_D_PARAM_TLR].push_back (-(float)ci->second.commulative_rate.liquid);
                                        cr.d_params[CONN_D_PARAM_TOI].push_back ((float)ci->second.commulative_rate.oil_inj);
                                        cr.d_params[CONN_D_PARAM_TWI].push_back ((float)ci->second.commulative_rate.water_inj);
                                        cr.d_params[CONN_D_PARAM_TGI].push_back ((float)ci->second.commulative_rate.gas_inj);
                    */
                  }

              }
          }
      }
    };

} //blue_sky
#endif //__BS_WELL_RESULTS_STORAGE_H
