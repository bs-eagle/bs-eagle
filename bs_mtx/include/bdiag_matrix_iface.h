#ifndef __DIAGONAL_MATRIX_IFACE_H
#define __DIAGONAL_MATRIX_IFACE_H
/** 
 * @file bdiag_matrix_iface.h
 * @brief BS interface for manipulation with block diagonal matrix
 * @author Oleg Borschuk
 * @date 2009-07-28
 */
#include "matrix_iface.h"

namespace blue_sky
{
  /** 
   * @brief interface class for block diagonal matrix storage and manipulation
   * <float, int>
   * <float, long>
   * <double, int>
   * <double, long>
   */
  template <class strat_t>
  class bdiag_matrix_iface: public matrix_iface<strat_t>
    {

    public:
      typedef matrix_iface <strat_t>                            base_t;
      typedef typename strat_t::fp_vector_type                  fp_vector_type_t;
      typedef typename strat_t::i_vector_type                   i_vector_type_t;
      typedef typename strat_t::fp_storage_vector_type          fp_storage_vector_type_t;
      typedef typename strat_t::fp_type_t                       fp_type_t;
      typedef typename strat_t::i_type_t                        i_type_t;
      typedef typename strat_t::fp_storage_type_t               fp_storage_type_t;

      typedef bs_array<fp_type_t>                               fp_array_t;
      typedef bs_array<i_type_t>                                i_array_t;
      typedef bs_array<fp_storage_type_t>                       fp_storage_array_t;

      typedef smart_ptr<fp_array_t, true>                       sp_fp_array_t;
      typedef smart_ptr<i_array_t, true>                        sp_i_array_t;
      typedef smart_ptr<fp_storage_array_t, true>               sp_fp_storage_array_t;

      typedef bdiag_matrix_iface<strat_t>                       this_t;

      typedef smart_ptr <base_t, true>                          sp_matrix_t;
      typedef smart_ptr<this_t, true>                           sp_bdiag_matrix_t;
      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
      //blue-sky class declaration
      //BLUE_SKY_TYPE_DECL_T(bdiag_matrix_iface);


    public:
      //! destructor
      virtual ~bdiag_matrix_iface () 
        {};

      // ----------------------------------------------
      // methods for memory initialization 
      // ----------------------------------------------
      
      //! allocate memory from given matrix
      virtual int init_by_matrix (sp_bdiag_matrix_t matrix) = 0;

      //! allocate memory using given new_n_rows and new_n_blok_size 
      virtual int init (const i_type_t new_n_rows, const i_type_t new_n_blok_size) = 0;

      // ------------------------------------------------
      // Data initialization methods
      // ------------------------------------------------

      //! copy all data from given matrix
      virtual int copy (sp_bdiag_matrix_t matrix) = 0;
      
      //! return pointer to the storage
      virtual sp_fp_storage_array_t get_diag () = 0;

    };

}//namespace blue_sky

#endif //__DIAGONAL_MATRIX_IFACE_H

