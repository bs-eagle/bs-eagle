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
  template <class strat_t>
  class bcsr_matrix_iface: public matrix_iface<strat_t>
    {

    public:
      typedef matrix_iface <strat_t>                            base_t;
      //typedef typename strat_t::fp_vector_type                  fp_vector_type_t;
      //typedef typename strat_t::i_vector_type                   i_vector_type_t;
      //typedef typename strat_t::fp_storage_vector_type          fp_storage_vector_type_t;
      typedef typename strat_t::fp_type_t                       fp_type_t;
      typedef typename strat_t::i_type_t                        i_type_t;
      typedef typename strat_t::fp_storage_type_t               fp_storage_type_t;

      typedef smart_ptr<bs_array<fp_type_t>, true>              sp_fp_array_t;
      typedef smart_ptr<bs_array<i_type_t>, true>               sp_i_array_t;
      typedef smart_ptr<bs_array<fp_storage_type_t>, true>      sp_fp_storage_array_t;

      //typedef ndarray<fp_type_t>                                fp_numpy_t;
      //typedef ndarray<fp_storage_type_t>                        fp_storage_numpy_t;
      //typedef ndarray<i_type_t>                                 i_numpy_t;
      //typedef smart_ptr<fp_numpy_t, true>                       sp_fp_array_t;
      //typedef smart_ptr<i_numpy_t, true>                        sp_i_numpy_t;
      //typedef smart_ptr<fp_storage_numpy_t, true>               sp_fp_storage_numpy_t;

      typedef bcsr_matrix_iface<strat_t>                        this_t;

      typedef smart_ptr <base_t, true>                          sp_matrix_t;
      typedef smart_ptr <this_t, true>                          sp_bcsr_matrix_t;

      
      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
      //blue-sky class declaration
      //BLUE_SKY_TYPE_DECL_T(bcsr_matrix_iface);


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
      virtual int init (const i_type_t new_n_rows, const i_type_t new_n_cols, const i_type_t new_n_blok_size,
                        const i_type_t new_n_non_zeros) = 0;

      //! allocate memory n_rows, n_cols, n_block_size, n_non_zeros, cols_ind, rows_ptr
      virtual int init_struct (const i_type_t new_n_rows, const i_type_t new_n_cols, const i_type_t new_n_non_zeros) = 0;

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
      virtual int init_by_arrays (const i_type_t n_rows_, const i_type_t n_cols_, const i_type_t n_block_size_,
                                  sp_i_array_t rows_, sp_i_array_t cols_, sp_fp_storage_array_t values_, bool make_copy) = 0;

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
      virtual int init_by_raw_ptr (const i_type_t n_rows_, const i_type_t n_cols_, const i_type_t n_block_size_,
                                  i_type_t *rows_, i_type_t *cols_, fp_storage_type_t *values_, bool make_copy) = 0;


      //! setup n_rows, and allocate memory for rows_ptr
      virtual int alloc_rows_ptr (const i_type_t new_n_rows) = 0;


      //! allocate memory for column indexes and set it with default values
      virtual int alloc_cols_ind (const i_type_t new_n_non_zeros) = 0;

      //! allocate memory for values and set it with default values
      virtual int alloc_values (const i_type_t new_n_non_zeros, const i_type_t new_n_blok_size) = 0;

      //! allocate memory for cols_ind and values
      virtual int alloc_cols_ind_and_values (const i_type_t new_n_non_zeros, const i_type_t new_n_block_size) = 0;

      //! set number of columns
      virtual void set_n_cols (const i_type_t new_n_cols) = 0;

      // ------------------------------------------------
      // Data initialization methods
      // ------------------------------------------------

      //! copy all data from given matrix
      virtual int copy (sp_bcsr_matrix_t matrix) = 0;

      // ------------------------------------------------
      // GET matrix private DATA methods
      // ------------------------------------------------
      //! return number of nonzeros elements
      virtual i_type_t get_n_non_zeros () const = 0;


      /** 
       * @return smart pointer to the rows_ptr
       */
      virtual sp_i_array_t get_rows_ptr () = 0;

      //virtual const i_vector_type_t &get_rows_ptr_const () const = 0;

      virtual sp_i_array_t get_cols_ind () = 0;

      //virtual const i_vector_type_t &get_cols_ind_const () const = 0;

      virtual sp_fp_storage_array_t get_values () = 0;

      //virtual const fp_storage_vector_type_t &get_values_const () const = 0;
      
      // ---------------------------------------------
      // INTERNAL checking
      // ---------------------------------------------

      //! check for correctness the structure of matrix (rows_ptr and cols_ind)
      virtual int internal_check () const = 0;

    };


}//namespace blue_sky


#endif //__BCSR_MATRIX_IFACE_H

