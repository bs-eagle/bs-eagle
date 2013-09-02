/**
 * \file scal_interpolate.h
 * \brief interpolate, binary_search and scale function for SCAL
 * \author Sergey Miryanov
 * \date 30.12.2008
 * */
#ifndef BS_SCAL_SCAL_INTERPOLATE_H_
#define BS_SCAL_SCAL_INTERPOLATE_H_

#include "interpolation_macro.h"

namespace blue_sky {

  #define  interpolation_threshold 1.0e-3

  /**
   * \brief binary_search
   *
   * \param[in] v value for search
   * \param[in] x vector of values search for
   * \return index of found value or x.size () if not found
   * */
  template <class vector_t, class comparator_t>
  BS_FORCE_INLINE size_t
  binary_search (typename vector_t::value_type v, const vector_t &c, const comparator_t &comparator)
  {
    size_t n  = c.size ();
    size_t il = 0;
    size_t iu = n;
    size_t im = 0;
    while (il != iu)
      {
        im = (il + iu) / 2;
        if (comparator (c[im], v))
          il = im + 1;
        else
          iu = im;
      }

    return il;
  }

  template <class vector_t, class dst_vector_t, class comparator_t>
  BS_FORCE_INLINE void
  interpolate (const vector_t &src, const dst_vector_t &dst, const typename vector_t::value_type &s, typename vector_t::value_type &d, const comparator_t &comparator)
  {
    size_t i = binary_search (s, src, comparator);

    if (i > 0 && i < (size_t)src.size ())
      {
        d = dst[i - 1] + (s - src[i - 1]) * (dst[i] - dst[i - 1]) / (src[i] - src[i - 1]);
      }
    else if (i == 0)
      {
        d = dst[0];
      }
    else if (i == (size_t)src.size ())
      {
        d = dst[i - 1];
      }
    else 
      {
        bs_throw_exception (boost::format ("i (%lu) out of range") % i);
      }
  }

  /**
    * \brief interpolate
    *
    * \param[in] s value for interpolate
    * \param[in] x vector of values X
    * \param[in] y vector of values Y
    * \param[out] cap (capillary) return value
    * \param[out] d_cap (d_cap/d_s) derivative
    * */
  template <typename item_t, typename data_vector_t, typename comparator_t>
  BS_FORCE_INLINE void
  interpolate (item_t s, const data_vector_t &x, const data_vector_t &y, item_t &cap, item_t *d_cap, const comparator_t &comparator)
  {
    size_t i = binary_search (s, x, comparator);

    if (i == (size_t)x.size ())
      {
        cap     = y[i - 1];
        if (d_cap) *d_cap = 0;
      }
    else if (i == 0)
      {
        cap     = y[0];
        if (d_cap) *d_cap = 0;
      }
    else
      {
        item_t d_cap_ = (y[i] - y[i-1]) / (x[i] - x[i - 1]);
        cap = y[i - 1] + (s - x[i - 1]) * d_cap_;
        if (d_cap) *d_cap = d_cap_;
      }
  }

  template <typename item_t, typename data_vector_t, typename comparator_t>
  BS_FORCE_INLINE void
  interpolate (item_t s, const data_vector_t &x, const data_vector_t &y, item_t spr_, int j_, item_t &kr, item_t &d_kr, const comparator_t &comparator)
  {
    size_t i = binary_search (s, x, comparator);

    if (i == (size_t)x.size ())
      {
        kr    = y[i - 1];
        d_kr  = 0;
      }
    else if (i == 0)
      {
        kr    = y[0];
        d_kr  = 0;
      }
    else
      {
        d_kr  = (y[i] - y[i-1]) / (x[i] - x[i - 1]);
        kr    = y[i - 1] + (s - x[i - 1]) * d_kr;
      }

    // so we remove loop
    if (fabs (s - spr_) < interpolation_threshold)
      {
        BS_ASSERT (y[j_] > 0.0);

        d_kr = (y[j_] - y[j_ - 1]) / (x[j_] - x[j_ - 1]);
      }
  }

  /**
   * scale (..., true)
   * */
  template <typename item_t>
  BS_FORCE_INLINE item_t
  scale_table (item_t a, item_t b, item_t c, item_t d, item_t s)
  {
    BS_ASSERT (c != a) (c) (a);
    return b + (s - a) * (d - b) / (c - a);
  }

  /**
   * scale (..., false)
   * */
  template <typename item_t>
  BS_FORCE_INLINE item_t
  scale_not_table (item_t a, item_t b, item_t c, item_t d, item_t s, bool /*is_table_value_calc = true*/)
  {
    BS_ASSERT (d != b) (d) (b);
    return a + (s - b) * (c - a) / (d - b);
  }

} 


#endif // #ifndef BS_SCAL_SCAL_INTERPOLATE_H_
