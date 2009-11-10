/**
 *
 * */
#ifndef BS_MAIN_LOOP_CALC_APPLY_WEFAC_H_
#define BS_MAIN_LOOP_CALC_APPLY_WEFAC_H_

namespace blue_sky
  {

  template <typename item_t, typename wefac_t>
  inline item_t
  apply_wefac (item_t item, wefac_t wefac)
  {
    return wefac > 0. ? (item_t)(wefac * item) : (item);
  }

} // namespace blue_sky


#endif  // #ifndef BS_MAIN_LOOP_CALC_APPLY_WEFAC_H_
