#ifndef __LU_DECOMPOSITION_H__
#define __LU_DECOMPOSITION_H__


/*!
  \file lu_decomposition.h
  \author Borschuk Oleg
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

//! minimum for devision
#define MIN_DIV 1.0e-24

//! default block size
#define DEF_BLOCK_SIZE 60

template <class fp_type, class i_type>
struct BS_API_PLUGIN lu_tools
  {
    // build LU decomposition for matrix block
    static i_type lu_block_decomposition (i_type matrix_size, fp_type *block, i_type block_size = 0);

    // find y in Ly = b
    static i_type lu_block_find_L_roots (i_type matrix_size, fp_type *block, fp_type *r_side, i_type block_size = 0);

    // find x in Ux = y
    static i_type lu_block_find_U_roots (i_type matrix_size, fp_type *block, fp_type *r_side, i_type block_size = 0);

    // find U in LU = A using known L
    static i_type lu_block_find_U (i_type matrix_size, fp_type *block_L, fp_type *block_U, i_type block_size, i_type ncol);

    // find L in LU = A using known U
    static i_type lu_block_find_L (i_type matrix_size, fp_type *block_L, fp_type *block_U, i_type block_size, i_type nrow);

    // upgrade block
    static i_type lu_block_upgrade (i_type matrix_size, fp_type *block_A, fp_type *block_L,
                                    fp_type *block_U, i_type block_size, i_type nrow, i_type ncol);

    // build lu decomposition using block method on one thread
    static i_type lu_decomposition (i_type matrix_size, fp_type *matrix, i_type block_size = DEF_BLOCK_SIZE);

    // update right side
    static i_type lu_upgrade_right_side (i_type matrix_size, fp_type *block, i_type nrow, i_type ncol,
                                         fp_type *roots, fp_type *r_side);
    // find y in Ly = b
    static i_type lu_find_L_roots (i_type matrix_size, fp_type *matrix, fp_type *r_side,
                                   i_type block_size = DEF_BLOCK_SIZE);

    // find x in Ux = y
    static i_type lu_find_U_roots (i_type matrix_size, fp_type *matrix, fp_type *r_side,
                                   i_type block_size = DEF_BLOCK_SIZE);
  };

#endif //__LU_DECOMPOSITION_H__
