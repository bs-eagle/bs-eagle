#ifndef __LU_DECOMPOSITION_H__
#define __LU_DECOMPOSITION_H__
/*!
  \file lu_decomposition.h
  \version 0.1
  \date 3.10.2004
  \brief file include declarations of functions to build LU decomposition and to solve matrix 
  
         +---+----------------+               matrix is an array of double where elements stored   
         |***|****************|               by strings 
         |***|block_size******|                 |a(0,0) a(0,1) a(0,2)|
         +---+----------------+ matrix_size     |a(1,0) a(1,1) a(1,2)| => 
         |**block_string******|                 |a(2,0) a(2,1) a(2,2)|
         |********************|
         |********************|                 => {a(0,0) a(0,1) a(0,2) a(1,0) a(1,1) a(1,2) ...} 
         +--------------------+
*/

#include <math.h>

//! minimum for devision
#define MIN_DIV 1.0e-24

//! default block size 
#define DEF_BLOCK_SIZE 60

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
  
  
template <class fp_type_t, class i_type_t, class fp_storage_type_t> 
class blu_solver_impl
{
  public:
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
    inline int 
    lu_block_decomposition (const int matrix_size, fp_storage_type_t *block, const int block_size) const
    {
      // declaration
      i_type_t i, j, k, nb;
      fp_storage_type_t *row = 0;
      fp_storage_type_t *uprow = 0;
      fp_storage_type_t diag;
      
      // check input variables
      if (matrix_size < 1 || !block)
        return -1;
      if (block_size < 1 || block_size > matrix_size)
        nb = matrix_size;
      else
        nb = block_size;
      
      // main loop
      for (k = 0; k < nb; ++k)
        {
          // calculate elements in string
          uprow = &(block[k * matrix_size]);
          // check diagonal element
          if (fabs (uprow[k]) < MIN_DIV)
            return -2;
          diag = 1.0 / uprow[k];
          for (j = k + 1; j < nb; ++j)
            uprow[j] *= diag;
          // upgrade other elements
          for (i = k + 1; i < nb; ++i)
            {
              row = &(block[i * matrix_size]);
              for (j = k + 1; j < nb; ++j)
                row[j] -= uprow[j] * row[k];

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
    inline int 
    lu_block_find_L_roots (const i_type_t matrix_size, const fp_storage_type_t *block, fp_type_t *r_side, const i_type_t block_size) const
    {
      i_type_t i, j, nb;
      // check input variables
      if (matrix_size < 1 || !block || !r_side)
        return -1;
      if (block_size < 1 || block_size > matrix_size)
        nb = matrix_size;
      else
        nb = block_size;
      
      for (i = 0; i < nb; ++i)
        {
          if (fabs (block[i * matrix_size + i]) < MIN_DIV)
            return -2;
          r_side[i] /= block[i * matrix_size + i];
          for (j = i + 1; j < nb; ++j)
            r_side[j] -= r_side[i] * block[j * matrix_size + i];
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
    inline int 
    lu_block_find_U_roots (const i_type_t matrix_size, const fp_storage_type_t *block, fp_type_t *r_side, const i_type_t block_size) const
    {
      i_type_t i, j, nb;

      // check input variables
      if (matrix_size < 1 || !block || !r_side)
        return -1;
      if (block_size < 1 || block_size > matrix_size)
        nb = matrix_size;
      else
        nb = block_size;
      
      for (i = nb - 1; i >= 0; --i)
        for (j = i - 1; j >= 0; --j)
          r_side[j] -= r_side[i] * block[j * matrix_size + i];
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
    inline int 
    lu_block_find_U (const i_type_t matrix_size, const fp_storage_type_t *block_L, fp_storage_type_t *block_U, 
                     const i_type_t block_size, const i_type_t ncol) const
    {
      i_type_t i, j, k;
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
    inline int 
    lu_block_find_L (const i_type_t matrix_size, fp_storage_type_t *block_L, 
                     const fp_storage_type_t *block_U, const i_type_t block_size, 
                     const i_type_t nrow) const
    {
      i_type_t i, j, k;
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
    inline int 
    lu_block_upgrade (const i_type_t matrix_size, fp_storage_type_t *block_A, 
                      const fp_storage_type_t *block_L,
                      const fp_storage_type_t *block_U, 
                      const i_type_t block_size, 
                      const i_type_t nrow, 
                      const i_type_t ncol) const
    {
      i_type_t i, j, k;
      fp_storage_type_t update_value;
      
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
    int 
    lu_decomposition (const i_type_t matrix_size, fp_storage_type_t *matrix, const i_type_t block_size) const
    {
      i_type_t i, j, k;
      i_type_t global_i;
      int ret_code;
      // last block size
      i_type_t last_block_size;
      // number of blocks in matrix string
      i_type_t number_of_blocks;
      // block on the corner
      fp_storage_type_t *corner_block = 0;
      // L block
      fp_storage_type_t *block_L = 0;
      // U block
      fp_storage_type_t *block_U = 0;
      // A block
      fp_storage_type_t *block_A = 0;
      
      
      // check input data
  if (matrix_size < 1 || !matrix)
    return -1;
  if (block_size >= matrix_size || block_size < 2)
    return lu_block_decomposition (matrix_size, matrix, matrix_size);
  
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
    inline int 
    lu_upgrade_right_side (const i_type_t matrix_size,  const fp_storage_type_t *block, 
                           const i_type_t nrow, const i_type_t ncol, 
                           const fp_type_t *roots, fp_type_t *r_side) const 
    {
      i_type_t i, j;
      fp_type_t upgrade_value = 0;
      const fp_storage_type_t *block_row = 0;
      
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
    int 
    lu_find_L_roots (const i_type_t matrix_size, const fp_storage_type_t *matrix, 
                     fp_type_t *r_side, const i_type_t block_size) const
    {
      i_type_t i;
      i_type_t k;
      int ret_code;
      // last block size
      i_type_t last_block_size;
      // number of blocks in matrix string
      i_type_t number_of_blocks;
      
      // check input data
      if (matrix_size < 1 || !matrix || !r_side)
        return -1;
      if (block_size >= matrix_size || block_size < 2)
        return lu_block_find_L_roots (matrix_size, matrix, r_side, matrix_size);
      
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
    int 
    lu_find_U_roots (const i_type_t matrix_size, const fp_storage_type_t *matrix, 
                     fp_type_t *r_side, const i_type_t block_size) const
    {
      i_type_t i;
      i_type_t k;
      int ret_code;
      // last block size
      i_type_t last_block_size;
      // number of blocks in matrix string
      i_type_t number_of_blocks;
      
      // check input data
      if (matrix_size < 1 || !matrix || !r_side)
        return -1;
      if (block_size >= matrix_size || block_size < 2)
        return lu_block_find_U_roots (matrix_size, matrix, r_side, matrix_size);
      
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
};
                     
#endif //__LU_DECOMPOSITION_H__
