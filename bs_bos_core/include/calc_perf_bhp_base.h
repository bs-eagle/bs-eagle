/**
 *       \file  calc_perf_bhp_base.h
 *      \brief  Base class (interface) for well perforation bhp calculation
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  18.11.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_CALC_PERF_BHP_BASE_H_
#define BS_CALC_PERF_BHP_BASE_H_

// WTF??
#include "well_results_storage.h"
#include "fip_results_storage.h"

namespace blue_sky
  {

  template <typename strategy_t>
  class calc_model;

  template <typename strategy_t>
  class well;

  /**
   * \class calc_perf_bhp_base
   * \brief Base class (interface) for well perforation bhp calculation
   * \todo  Should be renamed to calc_perf_bhp_iface
   * */
  template <typename strategy_t>
  class calc_perf_bhp_base : public objbase
    {
    public:

      typedef typename strategy_t::index_t      index_t;
      typedef typename strategy_t::item_t       item_t;

      typedef calc_model <strategy_t>           calc_model_t;
      typedef well <strategy_t>                 well_t;
      typedef rs_mesh_iface <strategy_t>        mesh_iface_t;

      typedef smart_ptr <calc_model_t, true>    sp_calc_model_t;
      typedef smart_ptr <well_t, true>          sp_well_t;
      typedef smart_ptr <mesh_iface_t, true>    sp_mesh_iface_t;

    public:
      /**
       * \brief  calc_perf_bhp_base dtor
       * */
      virtual ~calc_perf_bhp_base ()
      {}

      /**
       * \brief  For each well perforation calculates bhp value
       * \param  well well to calculate perforation bhp values
       * \param  calc_model
       * \param  mesh
       * */
      virtual void
      calculate (sp_well_t &well, const sp_calc_model_t &calc_model, const sp_mesh_iface_t &mesh) const = 0;
    };


} // namespace blue_sky


#endif  // #ifndef BS_CALC_PERF_BHP_BASE_H_

