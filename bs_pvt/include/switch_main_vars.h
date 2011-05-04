/**
 * \file switch_main_vars.h
 * \brief declaration of switch_main_vars function
 * \author Sergey Miryanov (from bos::ver_2_2)
 * \date 25.08.2008
 * */
#ifndef BS_SWITCH_MAIN_VARS_H_
#define BS_SWITCH_MAIN_VARS_H_

#include "pvt_oil.h"
#include "scal_interpolate.h"

namespace blue_sky
  {

  // TODO: pvt_oil::public

  struct switch_main_vars
    {
      //typedef boost::array <double, PHASE_TOTAL>  sat_3p_t;
      typedef t_double         item_t;
      typedef t_long        index_t;

      typedef pvt_oil 								pvt_oil_t;
      typedef pvt_dead_oil 						pvt_dead_oil_t;

      typedef smart_ptr <pvt_oil_t, true>         sp_pvt_oil_t;
      typedef smart_ptr <pvt_dead_oil_t, true>    sp_pvt_dead_oil_t;

      // TODO: in the future we can remove static modifier and made this class a stateful
      static
      void
      do_switch (bool is_w, bool is_g, bool is_o,
                 index_t i_o, index_t i_g, index_t i_w,
                 const sp_pvt_oil_t &pvt, item_t p,
                 // results
                 main_var_type &main_var,
                 item_t *sat_3p,
                 item_t &gas_oil_ratio,
                 bool &lsearch_force_newton_step,
                 item_t drsdt, item_t dt,
                 item_t old_gas_oil_ratio,
                 index_t &switch_to_sg_count,
                 index_t &switch_to_ro_count,
                 index_t &switch_to_momg_count,
                 index_t cell_index)
      {
        if (is_o & is_g)
          {
            /*item_t sat_w = is_w ? sat_3p[i_w] : 0;*/
            item_t sat_g = sat_3p[i_g];
            item_t sat_o = sat_3p[i_o];

            if (sat_g < 1e-7 && sat_o < 1e-7)
              momg_case (i_o, i_g, i_w, sat_3p, gas_oil_ratio, main_var, switch_to_momg_count, cell_index);
            else if (sat_o < EPS_DIFF || sat_g > 1.0 - EPS_DIFF)
              sg_case (gas_oil_ratio, main_var, switch_to_sg_count, cell_index);
            else if (main_var == FI_SG_VAR && sat_g >= 0)
              return;
            else
              other_case (i_o, i_w, i_g, pvt, p, main_var, sat_3p, gas_oil_ratio,
                          lsearch_force_newton_step, drsdt, dt, old_gas_oil_ratio,
                          switch_to_sg_count, switch_to_ro_count, cell_index);
          }
      }

private:

      static void
      momg_case (index_t i_o, index_t i_g, index_t i_w, item_t *sat_3p, item_t &gas_oil_ratio, main_var_type &main_var, index_t &count, index_t cell_index)
      {
        sat_3p[i_g] = 0;
        sat_3p[i_w] = 1;
        sat_3p[i_o] = 0;
        gas_oil_ratio = 0;
        if (main_var != FI_MOMG_VAR)
          {
#ifdef _DEBUG
            //BOSOUT (section::iters, level::debug) << boost::format ("switch-> momg: to sg [%d]") % cell_index << bs_end;
#endif
            main_var = FI_MOMG_VAR;
            ++count;
          }
      }

      static void
      sg_case (item_t &gas_oil_ratio, main_var_type &main_var, index_t &count, index_t /*cell_index*/)
      {
        gas_oil_ratio = 0;
        if (main_var != FI_SG_VAR)
          {
#ifdef _DEBUG
            //BOSOUT (section::iters, level::debug) << boost::format ("switch-> ro(1): to sg [%d]") % cell_index << bs_end;
#endif
            main_var = FI_SG_VAR;
            ++count;
          }
      }

      static void
      other_case (index_t i_o, index_t /*i_w*/, index_t i_g,
                  const sp_pvt_oil_t &pvt, item_t p,
                  main_var_type &main_var,
                  item_t *sat_3p,
                  item_t &gas_oil_ratio,
                  bool &lsearch_force_newton_step,
                  item_t /*drsdt*/, item_t /*dt*/,
                  item_t /*old_gas_oil_ratio*/,
                  index_t &switch_to_sg_count,
                  index_t &switch_to_ro_count, index_t /*cell_index*/)
      {
        BS_ASSERT (pvt);

        const pvt_oil_t::vector_t &pressure_ = pvt->get_pressure ();
        const pvt_oil_t::vector_t &gor_      = pvt->get_gor ();

        size_t il = 0, iu = 1;
        il = binary_search (p, pressure_, std::less <item_t> ());
        if (il != 0)
          {
            iu = il;
            --il;
          }

        item_t gor = gor_[il] + (gor_[iu] - gor_[il]) / (pressure_[iu] - pressure_[il]) * (p - pressure_[il]);
        if (main_var == FI_RO_VAR && gas_oil_ratio - gor <= 0)
          {
            sat_3p[i_g] = 0;
            return ;
          }

        if (main_var == FI_SG_VAR && sat_3p[i_g] < 0)
          {
            gas_oil_ratio = gor - 1e-4;
            if (gas_oil_ratio < EPS_DIFF)
              gas_oil_ratio = (item_t) EPS_DIFF;

#if 0
            if (drsdt > -EPS_DIFF)
              {
                if ((gas_oil_ratio - old_gas_oil_ratio) > drsdt * dt - EPS_DIFF)
                  {
                    gas_oil_ratio = old_gas_oil_ratio + drsdt * dt - EPS_DIFF;
                  }
              }
#endif //0

#ifdef _DEBUG
            //BOSOUT (section::iters, level::debug) << boost::format ("switch-> sg: to ro [%d]") % cell_index << bs_end;
#endif
            main_var = FI_RO_VAR;
            ++switch_to_ro_count;

            sat_3p[i_o] += sat_3p[i_g];
            sat_3p[i_g] = 0;
          }
        else if (main_var == FI_RO_VAR && gas_oil_ratio - gor > 0)
          {
#ifdef _DEBUG
            //BOSOUT (section::iters, level::debug) << boost::format ("switch-> ro: to sg [%d]") % cell_index << bs_end;
#endif
            main_var = FI_SG_VAR;
            sat_3p[i_g] = item_t (0.0001);
            //if (sat_3p[i_g] > (1 - sat_3p[i_w]) * 0.5)
            //  sat_3p[i_g] = (1 - sat_3p[i_w]) * 0.5;
            sat_3p[i_o] -= item_t (0.0001);

            gas_oil_ratio = gor;
            ++switch_to_sg_count;
          }

        lsearch_force_newton_step = 1;
        return ;
      }
    };

} // namespace blue_sky

#endif  // #ifndef BS_SWITCH_MAIN_VARS_H_

