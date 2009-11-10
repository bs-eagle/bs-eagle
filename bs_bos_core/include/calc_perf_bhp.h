/**
 * \file calc_perf_bhp.h
 * \brief perforation's bhp calculation
 * \author Sergey Miryanov
 * \date 18.11.2008
 * */
#ifndef BS_CALC_PERF_BHP_H_
#define BS_CALC_PERF_BHP_H_

#include "calc_perf_bhp_base.h"

namespace blue_sky
  {


  template <typename strategy_t>
  class BS_API_PLUGIN calc_perf_bhp : public calc_perf_bhp_base <strategy_t>
    {
    public:

      typedef calc_perf_bhp_base <strategy_t>     base_t;
      typedef calc_perf_bhp <strategy_t>          this_t;

      typedef typename base_t::sp_calc_model_t    sp_calc_model_t;
      typedef typename base_t::sp_well_t          sp_well_t;
      typedef typename base_t::sp_mesh_iface_t    sp_mesh_iface_t;
      typedef typename base_t::item_t             item_t;
      typedef typename base_t::index_t            index_t;

    public:

      void
      calculate (sp_well_t &well, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh) const;

      BLUE_SKY_TYPE_DECL_T (calc_perf_bhp <strategy_t>);

    };


} // namespace blue_sky


#endif  // #ifndef BS_CALC_PERF_BHP_H_

