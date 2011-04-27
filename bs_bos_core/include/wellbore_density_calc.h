/**
 *       \file  wellbore_density_calc.h
 *      \brief  Calculates perforation density
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  18.11.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_WELLBORE_DENSITY_CALC_H_
#define BS_WELLBORE_DENSITY_CALC_H_

#include "calc_perf_density_base.h"

namespace blue_sky
  {

  /**
   * \class wellbore_density_calc
   * \brief Calculates density for all well perforations
   * */
  class BS_API_PLUGIN wellbore_density_calc : public calc_perf_density_base
    {
    public:

      typedef calc_perf_density_base      base_t;
      typedef wellbore_density_calc       this_t;

      typedef base_t::item_t              item_t;
      typedef base_t::sp_well_t           sp_well_t;
      typedef base_t::sp_calc_model_t     sp_calc_model_t;

    public:

      /**
       * \brief  For each well perforation calculates density value
       * \param  well well to calculate perforation density value
       * \param  calc_model
       * */
      void
      calculate (sp_well_t &well, const sp_calc_model_t &calc_model) const;

      //! blue-sky type declaration
      BLUE_SKY_TYPE_DECL (wellbore_density_calc);

    protected:

      /**
       * \brief  For each well perforation calculates density value
       * \param  well well to calculate perforation density value
       * \param  calc_model
       * \param  bhp Precalculated BHP value
       * \return True if density calculated
       * */
      bool
      density_calc (sp_well_t &well, const sp_calc_model_t &calc_model, item_t bhp) const;

    };

} // namespace blue_sky


#endif  // #ifndef BS_WELLBORE_DENSITY_CALC_H_

