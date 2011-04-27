#ifndef __DENS_MATRIX_IFACE_H
#define __DENS_MATRIX_IFACE_H
/**
 * @file dens_matrix_iface.h
 * @brief Interface class for dense matrix
 * @date 2009-12-07
 */

#include "matrix_iface.h"

namespace blue_sky
{
  /**
   * @brief interface class for density matrix
   */
  class dens_matrix_iface: public matrix_iface
    {

    public:
      typedef matrix_iface                                      base_t;
      typedef dens_matrix_iface                                 this_t;

      typedef smart_ptr <base_t, true>                          sp_matrix_t;
      typedef smart_ptr<this_t, true>                           sp_dens_matrix_t;

      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
    public:
      //! destructor
      virtual ~dens_matrix_iface ()
        {};

      // ----------------------------------------------
      // methods for memory initialization
      // ----------------------------------------------

      /**
       * @brief initialize matrix
       *
       * @param matrix  -- <INPUT> referense to the matrix class
       *
       * @return 0 if success
       */
      virtual int init_by_matrix (sp_matrix_t matrix) = 0;

      /**
       * @brief initialize matrix
       *
       * @param new_n_rows              -- <INPUT> number of rows in matrix
       * @param new_n_cols              -- <INPUT> number of columns in matrix
       * @param calc_block_size         -- <INPUT> size of block used in matrix vector product (60 is a good choice)
       *
       * @return 0 if success
       */
      virtual int init (const t_long new_n_rows, const t_long new_n_cols, const t_long calc_block_size) = 0;

      /**
       * @brief make a copy of given matrix
       *
       * @param matrix  -- <INPUT> given matrix
       *
       * @return 0 if success
       */
      virtual int copy (sp_dens_matrix_t matrix) = 0;

      /**
       * @brief return reference to the vector of values
       */
      virtual spv_float get_values () = 0;

      /**
       * @brief return const reference to the vector of values
       */
      //virtual const fp_storage_vector_type_t &get_values_const () const = 0;


      /**
       * @brief return size of the calculation block in matrix vector product
       *        -1 -- use sequential algorithm
       */
      virtual t_long get_calc_block_size () const = 0;

      /**
       * @brief set size of the calculation block in matrix vector product
       *        if block_size * block_size can be put into the cache 1 memory
       *        this speedup calculation.
       *        60-80 is a good choice.
       *
       * @param block_size      -- <INPUT> new calculation block size (-1 for sequential calculation)
       */
      virtual void set_calc_block_size (const t_long block_size) = 0;
      // ---------------------------------------------
      // INTERNAL checking
      // ---------------------------------------------

      //! check for correctness the structure of matrix (rows_ptr and cols_ind)
      virtual int internal_check () const = 0;
    };
}//namespace blue_sky

#endif //__DENS_MATRIX_IFACE_H
