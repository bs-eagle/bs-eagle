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
  class bdiag_matrix_iface: public matrix_iface
    {

    public:
      typedef matrix_iface                                      base_t;
      typedef bdiag_matrix_iface                                this_t;

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
      virtual int init (const t_long new_n_rows, const t_long new_n_blok_size) = 0;

      // ------------------------------------------------
      // Data initialization methods
      // ------------------------------------------------

      //! copy all data from given matrix
      virtual int copy (sp_bdiag_matrix_t matrix) = 0;
      
      //! return pointer to the storage
      virtual spv_float get_diag () = 0;

    };

}//namespace blue_sky

#endif //__DIAGONAL_MATRIX_IFACE_H

