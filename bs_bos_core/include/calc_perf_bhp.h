/**
 *       \file  calc_perf_bhp.h
 *      \brief  Calculates perforation bhp value
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  18.11.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_CALC_PERF_BHP_H_
#define BS_CALC_PERF_BHP_H_

#include "calc_perf_bhp_base.h"

namespace blue_sky
  {


    /**
     * \class calc_perf_bhp
     * \brief Calculates perforation (well connection) bhp value, 
     *        implements calc_perf_bhp_base interface
     * */
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

      /**
       * \brief  For each well perforation calculates bhp value
       * \param  well well to calculate perforation bhp values
       * \param  calc_model
       * \param  mesh
       * */
      void
      calculate (sp_well_t &well, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh) const;

      //! blue-sky type declaration
      BLUE_SKY_TYPE_DECL_T (calc_perf_bhp <strategy_t>);

    };


} // namespace blue_sky


#endif  // #ifndef BS_CALC_PERF_BHP_H_

