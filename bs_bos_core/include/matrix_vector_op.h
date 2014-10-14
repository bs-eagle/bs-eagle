/**
 *       \file  matrix_vector_op.h
 *      \brief  Matrix-vector operations for derivs
 *              calculations
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  08.10.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_BOS_CORE_MATRIX_VECTOR_OP_H_
#define BS_BOS_CORE_MATRIX_VECTOR_OP_H_

#include "force_inline.h"

namespace blue_sky {

  /**
   * \class v_minus_vs_prod
   * \brief Incapsulates A = A - B * C operation, 
   *        A, B, C vectors
   * */
  template <size_t size>
  struct v_minus_vs_prod
  {
  };

  template <>
  struct v_minus_vs_prod <1>
  {
    template <typename a_t, typename b_t, typename c_t>
    static BS_FORCE_INLINE void
    eliminate (const b_t &b, const c_t &c, a_t &a)
    {
      a[0] -= b[0] * c;
    }
  };

  template <>
  struct v_minus_vs_prod <2>
  {
    template <typename a_t, typename b_t, typename c_t>
    static BS_FORCE_INLINE void
    eliminate (const b_t &b, const c_t &c, a_t &a)
    {
      a[0] -= b[0] * c;
      a[1] -= b[1] * c;
    }
  };

  template <>
  struct v_minus_vs_prod <3>
  {
    template <typename a_t, typename b_t, typename c_t>
    static BS_FORCE_INLINE void
    eliminate (const b_t &b, const c_t &c, a_t &a)
    {
      a[0] -= b[0] * c;
      a[1] -= b[1] * c;
      a[2] -= b[2] * c;
    }
  };

  /**
   * \class m_minus_vv_prod
   * \brief Incapsulates A = A - B * C operation,
   *        A - matrix, B, C - vectors
   * */
  template <size_t size>
  struct m_minus_vv_prod
  {
  };

  template <>
  struct m_minus_vv_prod <1>
  {
    template <typename a_t, typename b_t, typename c_t>
    static BS_FORCE_INLINE void
    eliminate (const b_t &b, const c_t &c, a_t &a)
    {
      a[0] -= b[0] * c[0];
    }
  };
  template <>
  struct m_minus_vv_prod <2>
  {
    template <typename a_t, typename b_t, typename c_t>
    static BS_FORCE_INLINE void
    eliminate (const b_t &b, const c_t &c, a_t &a)
    {
      a[0] -= b[0] * c[0];
      a[1] -= b[0] * c[1];
      a[2] -= b[1] * c[0];
      a[3] -= b[1] * c[1];
    }
  };
  template <>
  struct m_minus_vv_prod <3>
  {
    template <typename a_t, typename b_t, typename c_t>
    static BS_FORCE_INLINE void
    eliminate (const b_t &b, const c_t &c, a_t &a)
    {
      a[0] -= b[0] * c[0];
      a[1] -= b[0] * c[1];
      a[2] -= b[0] * c[2];
      a[3] -= b[1] * c[0];
      a[4] -= b[1] * c[1];
      a[5] -= b[1] * c[2];
      a[6] -= b[2] * c[0];
      a[7] -= b[2] * c[1];
      a[8] -= b[2] * c[2];
    }
  };

} // namespace blue_sky

#endif // #ifndef BS_BOS_CORE_MATRIX_VECTOR_OP_H_
