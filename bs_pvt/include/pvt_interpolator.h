/**
 * \file pvt_interpolator.h
 * \brief
 * \author Miryanov Sergey
 * \date 07.05.2008
 */
#ifndef BS_PVT_INTRPOLATOR_H_
#define BS_PVT_INTRPOLATOR_H_

namespace blue_sky
  {
  namespace pvt
    {

    template <typename item_t>
    BS_FORCE_INLINE item_t
    interpolate (item_t y1, item_t y2, item_t diff)
    {
      return y2 + (y2 - y1) * diff;
    }

    template <typename item_t>
    BS_FORCE_INLINE item_t
    interpolate (item_t x, item_t x1, item_t x2, item_t y1, item_t y2, item_t &diff)
    {
      diff = (x - x2) / (x2 - x1);
      return interpolate (y1, y2, diff);
    }

    template <typename item_t>
    BS_FORCE_INLINE item_t
    wtf (item_t a, item_t b)
    {
      return (b - a) / (b + a);
    }

    template <typename item_t>
    BS_FORCE_INLINE item_t
    interpolate_x (item_t y1, item_t y2, item_t diff)
    {
      return 1.0 / y1 + ((1.0 / y2) - (1.0 / y1)) * diff;
    }

    template <typename item_t>
    BS_FORCE_INLINE item_t
    interpolate_x (item_t x, item_t x1, item_t x2, item_t y1, item_t y2, item_t &diff)
    {
      diff = (x - x1) / (x2 - x1);
      return interpolate_x (y1, y2, diff);
    }

    template <typename item_t>
    BS_FORCE_INLINE item_t
    interpolate_x2 (item_t y1, item_t y2, item_t diff)
    {
      return 1.0 / y1 + ((1.0 / y1) - (1.0 / y2)) * diff;
    }

    template <typename item_t>
    BS_FORCE_INLINE item_t
    interpolate_x2 (item_t x, item_t x1, item_t x2, item_t y1, item_t y2, item_t &diff)
    {
      diff = (x - x1) / (x1 - x2);

      return interpolate_x2 (y1, y2, diff);
    }

    template <typename item_t>
    BS_FORCE_INLINE item_t
    interpolate_x3 (item_t y1, item_t y2, item_t diff)
    {
      return y1 + (y1 - y2) * diff;
    }

    template <typename item_t>
    BS_FORCE_INLINE item_t
    interpolate_x3 (item_t x, item_t x1, item_t x2, item_t y1, item_t y2, item_t &diff)
    {
      diff = (x - x2) / (x2 - x1);

      return interpolate_x3 (y1, y2, diff);
    }

  } // namespace pvt
} // namespace blue_sky


#endif  // #ifndef BS_PVT_INTRPOLATOR_H_
