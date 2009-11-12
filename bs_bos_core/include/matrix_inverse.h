/**
 *       \file  matrix_inverse.h
 *      \brief  Calculates inverse matrix
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  15.09.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_CALC_INVERSE_MATRIX_H_
#define BS_CALC_INVERSE_MATRIX_H_

namespace blue_sky
  {

  /**
   * \brief  Calculates inverse matrix
   * \param  n Size of matrix
   * \param  matrix Original matrix
   * \return 0 on success otherwise -1
   * */
  template <typename matrix_t>
  inline void
  calc_inverse_matrix (int n, matrix_t &matrix)
  {
    BS_ASSERT (false && "NOT IMPL YET");

    int i, j, k, l1, l2, l3, l;
    typename matrix_t::value_type temp;
    int *is = 0, *js = 0;

    if (!(is = new int[n]) || !(js = new int[n]))
      {
        throw bs_exception ("calc_inverse_matrix", "can't allocate memory");
      }

    //first step
    for (k = 0; k < n; k++)
      {
        temp = 0.;

        for (i = k; i < n; i++)
          for (j = k, l1 = i * n + k; j < n; j++, l1++)
            if (fabs (matrix[l1]) > temp)
              {
                temp = fabs (matrix[l1]);
                is[k] = i;
                js[k] = j;
              }

        if (temp < EPS_DIFF)
          {
            delete [] is;
            delete [] js;

            throw bs_exception ("calc_inverse_matrix", "can't find inverse matrix");
          }

        if (is[k] != k)
          for (j = 0, l1 = k * n, l2 = is[k] * n; j < n; j++, l1++, l2++)
            {
              temp = matrix[l1];
              matrix[l1] = matrix[l2];
              matrix[l2] = temp;
            }

        if (js[k] != k)
          for (i = 0, l1 = k, l2 = js[k]; i < n; i++, l1 += n, l2 += n)
            {
              temp = matrix[l1];
              matrix[l1] = matrix[l2];
              matrix[l2] = temp;
            }

        l = k * n + k;
        matrix[l] = 1. / matrix[l];

        for (j = 0, l1 = k * n; j < n; j++, l1++)
          if (j != k)
            matrix[l1] *= matrix[l];

        for (i = 0, l2 = k; i < n; i++, l2 += n)
          if (i != k)
            for (j = 0, l1 = i * n, l3 = k * n; j < n; j++, l1++, l3++)
              if (j != k)
                matrix[l1] -= matrix[l2] * matrix[l3];

        for (i = 0, l1 = k; i < n; i++, l1 += n)
          if (i != k)
            matrix[l1] *= - matrix[l];
      }

    //second step
    for (k = n - 1; k > -1; k--)
      {
        if (js[k] != k)
          for (j = 0, l1 = k * n, l2 = js[k] * n; j < n; j++, l1++, l2++)
            {
              temp = matrix[l1];
              matrix[l1] = matrix[l2];
              matrix[l2] = temp;
            }

        if (is[k] != k)
          for (i = 0, l1 = k, l2 = is[k]; i < n; i++, l1 += n, l2 += n)
            {
              temp = matrix[l1];
              matrix[l1] = matrix[l2];
              matrix[l2] = temp;
            }
      }

    delete [] is;
    delete [] js;
  }

} // namespace blue_sky

#endif  // #ifndef BS_CALC_INVERSE_MATRIX_H_


