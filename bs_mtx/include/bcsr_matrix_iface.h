#ifndef __BCSR_MATRIX_IFACE_H
#define __BCSR_MATRIX_IFACE_H
/** 
 * @file bcsr_matrix_iface.h
 * @brief BS interface for Block CSR matrix storage and manipulation
 * @author Oleg Borschuk
 * @date 2009-07-28
 */

#include "matrix_iface.h"

namespace blue_sky
{
  /** 
   * @brief interface class for block CSR matrix storage and manipulation
   */
  class BS_API_PLUGIN bcsr_matrix_iface: public matrix_iface
    {

    public:
      typedef matrix_iface                                      base_t;
      typedef bcsr_matrix_iface                                 this_t;

      typedef smart_ptr <base_t, true>                          sp_matrix_t;
      typedef smart_ptr <this_t, true>                          sp_bcsr_matrix_t;

      
      //-----------------------------------------
      //  METHODS
      //-----------------------------------------

    public:
      //! destructor
      virtual ~bcsr_matrix_iface ()
        {};

      // ----------------------------------------------
      // methods for memory initialization 
      // ----------------------------------------------
      
      //! allocate memory n_rows, n_cols, n_block_size, n_non_zeros, cols_ind, rows_ptr, values
      virtual int init_by_matrix (sp_bcsr_matrix_t matrix) = 0;

      //! allocate memory n_rows, n_cols, n_block_size, n_non_zeros, cols_ind, rows_ptr, values
      virtual int init (const t_long new_n_rows, const t_long new_n_cols, const t_long new_n_blok_size,
                        const t_long new_n_non_zeros) = 0;

      //! allocate memory n_rows, n_cols, n_block_size, n_non_zeros, cols_ind, rows_ptr
      virtual int init_struct (const t_long new_n_rows, const t_long new_n_cols, const t_long new_n_non_zeros) = 0;

      /** 
       * @brief initialize matrix by numpy arrays
       * 
       * @param n_rows_         -- <INPUT> number of rows
       * @param n_cols_         -- <INPUT> number of columns
       * @param n_block_size_   -- <INPUT> sizeof matrix block 
       * @param rows_           -- <INPUT> numpy array (rows_ptr information)
       * @param cols_           -- <INPUT> numpy array (cols_ind information)
       * @param values_         -- <INPUT> numpy array (values information)
       * @param make_copy       -- <INPUT> if true -- copy memory of rows_, cols_, values_; if false use the same memory
       * 
       * @return 0 if success
       */
      virtual int init_by_arrays (const t_long n_rows_, const t_long n_cols_, const t_long n_block_size_,
                                  spv_long rows_, spv_long cols_, spv_float values_, bool make_copy) = 0;

      /** 
       * @brief initialize matrix by raw pointers
       * 
       * @param n_rows_         -- <INPUT> number of rows
       * @param n_cols_         -- <INPUT> number of columns
       * @param n_block_size_   -- <INPUT> sizeof matrix block 
       * @param rows_           -- <INPUT> rows_ptr information
       * @param cols_           -- <INPUT> cols_ind information
       * @param values_         -- <INPUT> values information
       * @param make_copy       -- <INPUT> if true -- copy memory of rows_, cols_, values_; if false use the same memory
       * 
       * @return 0 if success
       */
      virtual int init_by_raw_ptr (const t_long n_rows_, const t_long n_cols_, const t_long n_block_size_,
                                  t_long *rows_, t_long *cols_, t_float *values_, bool make_copy) = 0;


      //! setup n_rows, and allocate memory for rows_ptr
      virtual int alloc_rows_ptr (const t_long new_n_rows) = 0;


      //! allocate memory for column indexes and set it with default values
      virtual int alloc_cols_ind (const t_long new_n_non_zeros) = 0;

      //! allocate memory for values and set it with default values
      virtual int alloc_values (const t_long new_n_non_zeros, const t_long new_n_blok_size) = 0;

      //! allocate memory for cols_ind and values
      virtual int alloc_cols_ind_and_values (const t_long new_n_non_zeros, const t_long new_n_block_size) = 0;

      //! set number of columns
      virtual void set_n_cols (const t_long new_n_cols) = 0;

      // ------------------------------------------------
      // Data initialization methods
      // ------------------------------------------------

      //! copy all data from given matrix
      virtual int copy (sp_bcsr_matrix_t matrix) = 0;

      // ------------------------------------------------
      // GET matrix private DATA methods
      // ------------------------------------------------
      //! return number of nonzeros elements
      virtual t_long get_n_non_zeros () const = 0;


      /** 
       * @return smart pointer to the rows_ptr
       */
      virtual spv_long get_rows_ptr () = 0;

      //virtual const i_vector_type_t &get_rows_ptr_const () const = 0;

      virtual spv_long get_cols_ind () = 0;

      //virtual const i_vector_type_t &get_cols_ind_const () const = 0;

      virtual spv_float get_values () = 0;

      //virtual const fp_storage_vector_type_t &get_values_const () const = 0;
      
      // ---------------------------------------------
      // INTERNAL checking
      // ---------------------------------------------

      //! check for correctness the structure of matrix (rows_ptr and cols_ind)
      virtual int internal_check () const = 0;

    };


}//namespace blue_sky


#endif //__BCSR_MATRIX_IFACE_H

