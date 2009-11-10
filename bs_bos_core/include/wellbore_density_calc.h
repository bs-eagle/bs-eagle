/**
 * \file wellbore_density_calc.h
 * \brief perforation's density calculation (like as wellbore_density_calc in old code)
 * \author Sergey Miryanov
 * \date 18.11.2008
 * */
#ifndef BS_WELLBORE_DENSITY_CALC_H_
#define BS_WELLBORE_DENSITY_CALC_H_

#include "calc_perf_density_base.h"

namespace blue_sky
  {

  template <typename strategy_t>
  class BS_API_PLUGIN wellbore_density_calc : public calc_perf_density_base <strategy_t>
    {
    public:

      typedef calc_perf_density_base <strategy_t>     base_t;
      typedef wellbore_density_calc <strategy_t>      this_t;

      typedef typename base_t::item_t                 item_t;
      typedef typename base_t::sp_well_t              sp_well_t;
      typedef typename base_t::sp_calc_model_t        sp_calc_model_t;

    public:

      void
      calculate (sp_well_t &well, const sp_calc_model_t &calc_model) const;

      BLUE_SKY_TYPE_DECL_T (wellbore_density_calc);

    protected:

      bool
      density_calc (sp_well_t &well, const sp_calc_model_t &calc_model, item_t bhp) const;

    };

} // namespace blue_sky


#endif  // #ifndef BS_WELLBORE_DENSITY_CALC_H_

