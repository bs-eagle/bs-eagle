/** 
 * \file matrix_sum_spec.h
 * \brief specialization of sum and copy for small density mx
 * \author Sergey Miryanov
 * \date 13.01.2009
 * */
#ifndef BS_MATRIX_SUM_SPEC_H_
#define BS_MATRIX_SUM_SPEC_H_

namespace blue_sky {

  template <size_t block_size>
  struct sum_impl
  {
  };

  template <>
  struct sum_impl <1>
  {
    template <typename item_t>
    static void
    sum (const item_t *src, item_t *dst)
    {
      dst[0] += src[0];
    }

    template <typename item_t>
    static void
    copy (const item_t *src, item_t *dst)
    {
      dst[0] = src[0];
    }
  };

  template <>
  struct sum_impl <2>
  {
    template <typename item_t>
    static void
    sum (const item_t *src, item_t *dst)
    {
      dst[0] += src[0];
      dst[1] += src[1];
      dst[2] += src[2];
      dst[3] += src[3];
    }
    template <typename item_t>
    static void
    copy (const item_t *src, item_t *dst)
    {
      dst[0] = src[0];
      dst[1] = src[1];
      dst[2] = src[2];
      dst[3] = src[3];
    }
  };

  template <>
  struct sum_impl <3>
  {
    template <typename item_t>
    static void
    sum (const item_t *src, item_t *dst)
    {
      dst[0] += src[0];
      dst[1] += src[1];
      dst[2] += src[2];
      dst[3] += src[3];
      dst[4] += src[4];
      dst[5] += src[5];
      dst[6] += src[6];
      dst[7] += src[7];
      dst[8] += src[8];
    }
    template <typename item_t>
    static void
    copy (const item_t *src, item_t *dst)
    {
      dst[0] = src[0];
      dst[1] = src[1];
      dst[2] = src[2];
      dst[3] = src[3];
      dst[4] = src[4];
      dst[5] = src[5];
      dst[6] = src[6];
      dst[7] = src[7];
      dst[8] = src[8];
    }
  };
}


#endif  // #ifndef BS_MATRIX_SUM_SPEC_H_
