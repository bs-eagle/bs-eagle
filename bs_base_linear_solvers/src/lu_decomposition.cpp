/*!
  \file lu_decomposition.cpp
  \author Borschuk Oleg
  \version 0.1
  \date 3.10.2004
  \brief file include functions to build LU decomposition and to solve matrix
*/
#include "bs_base_linear_solvers_stdafx.h"

#include <math.h>

#include "lu_decomposition.h"

#define LU_2x2(m_str_0, m_str_1)                \
  m_str_0[1] /= m_str_0[0];                     \
  m_str_1[1] -= m_str_0[1] * m_str_1[0];

#define LU_3x3(m_str_0, m_str_1, m_str_2)       \
  LU_2x2(m_str_0, m_str_1)                      \
  m_str_0[2] /= m_str_0[0];                     \
  m_str_1[2] -= m_str_0[1] * m_str_1[0];        \
  m_str_1[2] /= m_str_1[1];                     \
  m_str_2[1] -= m_str_0[1] * m_str_2[0];        \
  m_str_2[2] -= m_str_0[2] * m_str_2[0];        \
  m_str_2[2] -= m_str_2[1] * m_str_1[2];




/*!
  \fn int lu_block_decomposition (int matrix_size, double *block, int block_size)
  \brief build lu decomposition for matrix lock
  \param matrix_size -- size of matrix
  \param block -- poiter to the first element of block, if need to decompose full
         matrix send pointer to the matrix and set block_size to zero,
         matrix must be stored by strings
  \param block_size -- size of block
  \return 0 if success
  \return < 0 if error occur
*/
template <class fp_type, class i_type> i_type
lu_tools<fp_type, i_type>::lu_block_decomposition (i_type matrix_size, fp_type *block, i_type block_size)
{
  using namespace blue_sky;

  BS_ASSERT (block);
  BS_ASSERT (matrix_size >= 1) (matrix_size);
  BS_ASSERT (block_size == matrix_size) (block_size) (matrix_size);

  // declaration
  i_type i, j, k;

  // main loop
  fp_type *uprow = &block[0];
  for (k = 0; k < block_size; ++k, uprow += matrix_size)
    {
      // check diagonal element
      if (fabs (uprow[k]) < MIN_DIV)
        return -2;

      fp_type diag = (fp_type)(1.0 / uprow[k]);
      for (j = k + 1; j < block_size; ++j)
        {
          uprow[j] *= diag;
        }

      // upgrade other elements
      fp_type *row = &block[(k + 1) * matrix_size];
      for (i = k + 1; i < block_size; ++i, row += matrix_size)
        {
          fp_type row_k = row[k];
          for (j = k + 1; j < block_size; ++j)
            {
              row[j] -= uprow[j] * row_k;
            }
        }
    }
  return 0;
}


/*!
  \fn int lu_block_find_L_roots (int matrix_size, double *block,
                                 double *r_side, int block_size = 0)
  \brief find y in Ly = b
  \param matrix_size -- size of matrix
  \param block -- poiter to the first element of block, if need to find roots in full
         matrix send pointer to the matrix and set block_size to zero,
         matrix must be stored by strings
  \param r_side -- pointer to the first element of block in right side array
  \param block_size -- size of block
  \return 0 if success
  \return < 0 if error occur
*/
template <class fp_type, class i_type> i_type
lu_tools<fp_type, i_type>::lu_block_find_L_roots (i_type matrix_size, fp_type *block, fp_type *r_side, i_type block_size)
{
  using namespace blue_sky;
  BS_ASSERT (block);
  BS_ASSERT (r_side);
  BS_ASSERT (matrix_size >= 1) (matrix_size);
  BS_ASSERT (block_size >= 1) (block_size);
  BS_ASSERT (block_size == matrix_size) (block_size) (matrix_size);

  for (i_type i = 0; i < block_size; ++i)
    {
      BS_ASSERT (fabs (block [i * matrix_size + i]) >= MIN_DIV) (i) (matrix_size) (block [i * matrix_size + i]) (MIN_DIV);
      //if (fabs (block[i * matrix_size + i]) < MIN_DIV)
      //  {
      //    return -2;
      //  }

      r_side[i] /= block[i * matrix_size + i];
      fp_type r_side_i = r_side[i];
      fp_type *b = &block[(i + 1) * matrix_size + i];
      for (i_type j = i + 1; j < block_size; ++j, b += matrix_size)
        {
          r_side[j] -= r_side_i * b[0];
        }
    }
  return 0;
}

/*!
  \fn int lu_block_find_U_roots (int matrix_size, double *block,
                                 double *r_side, int block_size = 0)
  \brief find x in Lx = y
  \param matrix_size -- size of matrix
  \param block -- poiter to the first element of block, if need to find roots in full
         matrix send pointer to the matrix and set block_size to zero,
         matrix must be stored by strings
  \param r_side -- pointer to the first element of block in right side array
  \param block_size -- size of block
  \return 0 if success
  \return < 0 if error occur
*/
template <class fp_type, class i_type> i_type
lu_tools<fp_type, i_type>::lu_block_find_U_roots (i_type matrix_size, fp_type *block, fp_type *r_side, i_type block_size)
{
  using namespace blue_sky;
  BS_ASSERT (block);
  BS_ASSERT (r_side);
  BS_ASSERT (matrix_size >= 1) (matrix_size);
  BS_ASSERT (block_size >= 1) (block_size);
  BS_ASSERT (block_size == matrix_size) (block_size) (matrix_size);

  for (i_type i = block_size - 1; i >= 0; --i)
    {
      fp_type r_side_i = r_side[i];
      fp_type *b = &block[(i - 1) * matrix_size + i];
      for (i_type j = i - 1; j >= 0; --j, b -= matrix_size)
        {
          r_side[j] -= r_side_i * b[0];
        }
    }
  return 0;
}

/*!
  \fn int lu_block_find_U (int matrix_size, double *block_L,
                           double *block_U, int nrow, int ncol)
  \brief find U in LU = A using known L and A
  \param matrix_size -- size of matrix
  \param block_L -- poiter to the first element of block L
  \param block_U -- pointer to the first element of block A and U
  \param block_size -- size of block
  \param ncol -- number of columns in U block
  \return 0 if success
  \return < 0 if error occur
*/
template <class fp_type, class i_type> i_type
lu_tools<fp_type, i_type>::lu_block_find_U (i_type matrix_size, fp_type *block_L, fp_type *block_U, i_type block_size, i_type ncol)
{
  i_type i, j, k;
  // check input variables
  if (matrix_size < 1 || !block_L || !block_U)
    return -1;
  for (k = 0; k < ncol; ++k)
    {
      for (i = 0; i < block_size; ++i)
        {
          //todo: not need to check values every time
          if (fabs (block_L[i * matrix_size + i]) < MIN_DIV)
            return -2;
          block_U[i * matrix_size + k] /= block_L[i * matrix_size + i];
          for (j = i + 1; j < block_size; ++j)
            block_U[j * matrix_size + k] -= block_U[i * matrix_size + k]
                                            * block_L[j * matrix_size + i];
        }

    }
  return 0;
}

/*!
  \fn int lu_block_find_L (int matrix_size, double *block_L,
                           double *block_U, int nrow, int ncol)
  \brief find L in LU = A using known U and A
  \param matrix_size -- size of matrix
  \param block_L -- poiter to the first element of block L and A
  \param block_U -- pointer to the first element of block U
  \param nrow -- number of rows in L block
  \param block_size -- size of block
  \return 0 if success
  \return < 0 if error occur
*/
template <class fp_type, class i_type> i_type
lu_tools<fp_type, i_type>::lu_block_find_L (i_type matrix_size, fp_type *block_L, fp_type *block_U, i_type block_size, i_type nrow)
{
  i_type i, j, k;
  // check input variables
  if (matrix_size < 1 || !block_L || !block_U)
    return -1;
  for (k = 0; k < nrow; ++k)
    {
      for (i = 0; i < block_size; ++i)
        {
          for (j = i + 1; j < block_size; ++j)
            block_L[j + k * matrix_size] -= block_L[i + k * matrix_size]
                                            * block_U[j + i * matrix_size];
        }

    }
  return 0;
}

/*!
  \fn int lu_block_upgrade (int matrix_size, double *block_A, double *block_L,
                            double *block_U, int block_size, int nrow, int ncol)
  \brief update block_A using given block_L and block_U
  \param matrix_size -- size of matrix
  \param block_A -- pointer to the first element of block A
  \param block_L -- poiter to the first element of block L
  \param block_U -- pointer to the first element of block U
  \param block_size -- size of block
  \param nrow -- number of rows in L block
  \param ncol -- number of columns in U block
  \return 0 if success
  \return < 0 if error occur
*/
template <class fp_type, class i_type> i_type
lu_tools<fp_type, i_type>::lu_block_upgrade (i_type matrix_size, fp_type *block_A, fp_type *block_L,
    fp_type *block_U, i_type block_size, i_type nrow, i_type ncol)
{
  i_type i, j, k;
  fp_type update_value;

  // check input data
  if (!block_A || !block_L || !block_U || matrix_size < 1
      || block_size < 1 || nrow < 1 || ncol < 1)
    return -1;

  // main loop
  for (i = 0; i < nrow; ++i)
    {
      for (j = 0; j < ncol; ++j)
        {
          update_value = 0;
          for (k = 0; k < block_size; ++k)
            update_value += block_L[i * matrix_size + k] * block_U[k * matrix_size +j];
          block_A[i * matrix_size + j] -= update_value;
        }
    }
  return 0;
}

/*!
  \fn int lu_decomposition (int matrix_size, double *matrix, int block_size)
  \brief build lu decomposition using block method on one thread
  \param matrix_size -- size of matrix
  \param matrix -- matrix to decompose
  \param block_size -- size of block (default value is 60)
  \return 0 if success
  \return < 0 if error occur
*/
template <class fp_type, class i_type> i_type
lu_tools<fp_type, i_type>::lu_decomposition (i_type matrix_size, fp_type *matrix, i_type block_size)
{
  i_type i, j, k;
  i_type global_i;
  i_type ret_code;
  // last block size
  i_type last_block_size;
  // number of blocks in matrix string
  i_type number_of_blocks;
  // block on the corner
  fp_type *corner_block = 0;
  // L block
  fp_type *block_L = 0;
  // U block
  fp_type *block_U = 0;
  // A block
  fp_type *block_A = 0;


  // check input data
  if (matrix_size < 1 || !matrix)
    return -1;
  if (block_size >= matrix_size || block_size < 2)
    return lu_block_decomposition (matrix_size, matrix);

  // calculate number of blocks and size of last block
  number_of_blocks = matrix_size / block_size;
  last_block_size = matrix_size % block_size;
  if (last_block_size > 0)
    ++number_of_blocks;
  else
    last_block_size = block_size;

  // main loop
  for (i = 0; i < number_of_blocks; ++i)
    {
      // find LU decomposition for block[i, i]
      global_i = i * block_size;
      corner_block = &(matrix[global_i * matrix_size + global_i]);
      ret_code = lu_block_decomposition (matrix_size,
                                         corner_block,
                                         ((i + 1) < number_of_blocks ? block_size : last_block_size));
      if (ret_code) return ret_code;

      // find all U blocks in block string
      for (k = i + 1; k < number_of_blocks; ++k)
        {
          ret_code = lu_block_find_U (matrix_size,
                                      corner_block,
                                      &(matrix[global_i * matrix_size + k * block_size]),
                                      block_size,
                                      ((k + 1) < number_of_blocks ? block_size : last_block_size));
          if (ret_code) return ret_code;
        }
      // find all L blocks in block column
      for (k = i + 1; k < number_of_blocks; ++k)
        {
          ret_code = lu_block_find_L (matrix_size,
                                      &(matrix[k * block_size * matrix_size + global_i]),
                                      corner_block,
                                      block_size,
                                      ((k + 1) < number_of_blocks ? block_size : last_block_size));
          if (ret_code) return ret_code;
        }
      // update other blocks
      for (k = i + 1; k < number_of_blocks; ++k)
        {
          block_L = &(matrix[k * block_size * matrix_size + global_i]);
          for (j = i + 1; j < number_of_blocks; ++j)
            {
              block_U = &(matrix[global_i * matrix_size + j * block_size]);
              block_A = &(matrix[k * block_size * matrix_size + j * block_size]);
              ret_code = lu_block_upgrade (matrix_size, block_A, block_L, block_U, block_size,
                                           ((k + 1) < number_of_blocks ? block_size : last_block_size),
                                           ((j + 1) < number_of_blocks ? block_size : last_block_size));
              if (ret_code) return ret_code;
            }
        }
    }
  return 0;
}

/*!
  \fn int lu_upgrade_right_side (int matrix_size, double *block, int nrow, int ncol,
                                 double *roots, double *r_side)
  \brief update right side
  \param matrix_size -- size of matrix
  \param block -- pointer to the first element of block
  \param nrow -- number of rows in upgrade block
  \param ncol -- number of columns in upgrade block
  \param roots -- known vector of roots
  \param r_side -- vector of upgrade roots
  \return 0 if success
  \return < 0 if error occur
*/
template <class fp_type, class i_type> i_type
lu_tools<fp_type, i_type>::lu_upgrade_right_side (i_type matrix_size, fp_type *block, i_type nrow, i_type ncol,
    fp_type *roots, fp_type *r_side)
{
  i_type i, j;
  fp_type upgrade_value = 0;
  fp_type *block_row = 0;

  // check input data
  if (matrix_size < 1 || !block || nrow < 1 || ncol < 1 || !roots || !r_side)
    return -1;

  for (i = 0; i < nrow; ++i)
    {
      upgrade_value = 0;
      block_row = &(block[i * matrix_size]);
      for (j = 0; j < ncol; ++j)
        upgrade_value += block_row[j] * roots[j];
      r_side[i] -= upgrade_value;
    }
  return 0;
}


/*!
  \fn int lu_find_L_roots (int matrix_size, double *matrix, double *r_side,
                           int block_size = DEF_BLOCK_SIZE)
  \brief find y in Ly = b
  \param matrix_size -- size of matrix
  \param matrix -- pointer to the lu matrix
  \param r_side -- given right side
  \param block_size -- block size, default value DEF_BLOCK_SIZE
  \return 0 if success
  \return < 0 if error occur
*/
template <class fp_type, class i_type> i_type
lu_tools<fp_type, i_type>::lu_find_L_roots (i_type matrix_size, fp_type *matrix, fp_type *r_side, i_type block_size)
{
  i_type i;
  i_type k;
  i_type ret_code;
  // last block size
  i_type last_block_size;
  // number of blocks in matrix string
  i_type number_of_blocks;

  // check input data
  if (matrix_size < 1 || !matrix || !r_side || block_size < 1)
    return -1;

  // calculate number of blocks and size of last block
  number_of_blocks = matrix_size / block_size;
  last_block_size = matrix_size % block_size;
  if (last_block_size > 0)
    ++number_of_blocks;
  else
    last_block_size = block_size;

  // main loop
  for (i = 0; i < number_of_blocks; ++i)
    {//libconfig::ParseException
      ret_code = lu_block_find_L_roots (matrix_size,
                                        &(matrix[i * block_size * (matrix_size + 1)]),
                                        &(r_side[i * block_size]),
                                        ((i + 1) < number_of_blocks ? block_size : last_block_size));
      if (ret_code) return ret_code;

      // update all blocks in block column
      for (k = i + 1; k < number_of_blocks; ++k)
        {
          ret_code = lu_upgrade_right_side (matrix_size,
                                            &(matrix[(k * matrix_size + i) * block_size]),
                                            ((k + 1) < number_of_blocks ? block_size :last_block_size),
                                            ((i + 1) < number_of_blocks ? block_size :last_block_size),
                                            &(r_side[i * block_size]),
                                            &(r_side[k * block_size]));
          if (ret_code) return ret_code;

        }
    }
  return 0;
}

/*!
  \fn int lu_find_U_roots (int matrix_size, double *matrix, double *r_side,
                           int block_size = DEF_BLOCK_SIZE)
  \brief find x in Ux = y
  \param matrix_size -- size of matrix
  \param matrix -- pointer to the lu matrix
  \param r_side -- given right side
  \param block_size -- block size, default value DEF_BLOCK_SIZE
  \return 0 if success
  \return < 0 if error occur
*/
template <class fp_type, class i_type> i_type
lu_tools<fp_type, i_type>::lu_find_U_roots (i_type matrix_size, fp_type *matrix, fp_type *r_side, i_type block_size)
{
  i_type i;
  i_type k;
  i_type ret_code;
  // last block size
  i_type last_block_size;
  // number of blocks in matrix string
  i_type number_of_blocks;

  // check input datalibconfig::ParseException
  if (matrix_size < 1 || !matrix || !r_side || block_size < 1)
    return -1;

  // calculate number of blocks and size of last block
  number_of_blocks = matrix_size / block_size;
  last_block_size = matrix_size % block_size;
  if (last_block_size > 0)
    ++number_of_blocks;
  else
    last_block_size = block_size;

  // main loop
  for (i = number_of_blocks - 1; i >= 0; --i)
    {
      ret_code = lu_block_find_U_roots (matrix_size,
                                        &(matrix[i * block_size * (matrix_size + 1)]),
                                        &(r_side[i * block_size]),
                                        ((i + 1) < number_of_blocks ? block_size : last_block_size));
      if (ret_code) return ret_code;

      // update all blocks in block column
      for (k = i - 1; k >= 0; --k)
        {
          ret_code = lu_upgrade_right_side (matrix_size,
                                            &(matrix[(k * matrix_size + i) * block_size]),
                                            ((k + 1) < number_of_blocks ? block_size :last_block_size),
                                            ((i + 1) < number_of_blocks ? block_size :last_block_size),
                                            &(r_side[i * block_size]),
                                            &(r_side[k * block_size]));
          if (ret_code) return ret_code;

        }
    }
  return 0;
}

template struct lu_tools<double, int>;
template struct lu_tools<float, int>;
