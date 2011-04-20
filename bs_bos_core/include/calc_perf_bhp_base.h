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

namespace blue_sky
  {

  class calc_model;
  class well;

  /**
   * \class calc_perf_bhp_base
   * \brief Base class (interface) for well perforation bhp calculation
   * \todo  Should be renamed to calc_perf_bhp_iface
   * */
  class calc_perf_bhp_base : public objbase
    {
    public:

      typedef t_long                            index_t;
      typedef t_double                          item_t;

      typedef calc_model                        calc_model_t;
      typedef well                              well_t;
      typedef rs_mesh_iface                     mesh_iface_t;

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

