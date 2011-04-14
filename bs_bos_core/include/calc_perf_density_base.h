/**
 *       \file  calc_perf_density_base.h
 *      \brief  Base class (interface) for well perforation density calculation
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  18.11.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_CALC_PERF_DENSITY_BASE_H_
#define BS_CALC_PERF_DENSITY_BASE_H_


namespace blue_sky
  {
  class calc_model;
  class well;

  /**
   * \class calc_perf_density_base
   * \brief Base class (interface) for well perforation density calculation
   * \todo  Should be renamed to calc_perf_density_iface
   * */
  class calc_perf_density_base : public objbase
    {
    public:

      typedef t_long                            index_t;
      typedef t_double                          item_t;

      typedef calc_model                        calc_model_t;
      typedef well                              well_t;

      typedef smart_ptr <calc_model_t, true>    sp_calc_model_t;
      typedef smart_ptr <well_t, true>          sp_well_t;

    public:
      /**
       * \brief  calc_perf_density_base dtor
       * */
      virtual ~calc_perf_density_base () {}

      /**
       * \brief  For each well perforation calculates density value
       * \param  well well to calculate perforation density value
       * \param  calc_model
       * */
      virtual void
      calculate (sp_well_t &well, const sp_calc_model_t &calc_model) const = 0;
    };


} // namespace blue_sky


#endif // #ifndef BS_CALC_PERF_DENSITY_BASE_H_

