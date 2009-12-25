/**
 *       \file  apply_wefac.h
 *      \brief  applies wefac factor to passed value
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  22.12.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_MAIN_LOOP_CALC_APPLY_WEFAC_H_
#define BS_MAIN_LOOP_CALC_APPLY_WEFAC_H_

namespace blue_sky
  {

  /**
   * \brief  applies wefac factor to passed value
   * \param  item value to wefac applied
   * \param  wefac wefac factor, applies only if wefac > 0.0
   * \return return wefac * item
   * */
  template <typename item_t, typename wefac_t>
  inline item_t
  apply_wefac (item_t item, wefac_t wefac)
  {
    return wefac > 0. ? (item_t)(wefac * item) : (item);
  }

} // namespace blue_sky


#endif  // #ifndef BS_MAIN_LOOP_CALC_APPLY_WEFAC_H_
