#ifndef __BCSR_AMG_MATRIX_IFACE_H
#define __BCSR_AMG_MATRIX_IFACE_H
/** 
 * @file bcsr_amg_matrix_iface.h
 * @brief BS interface for Block CSR matrix storage and manipulation
 * @date 2009-07-28
 */

#include "bcsr_matrix_iface.h"

namespace blue_sky
{
  /** 
   * @brief interface class for block CSR matrix storage and manipulation
   */
  class bcsr_amg_matrix_iface: public bcsr_matrix_iface
    {

    public:
      typedef matrix_iface                                      base_t;
      typedef bcsr_matrix_iface                                 bcsr_matrix_t;

      typedef smart_ptr <base_t, true>                          sp_matrix_t;
      typedef smart_ptr<bcsr_matrix_t, true>                    sp_bcsr_matrix_t;
      
      //-----------------------------------------
      //  METHODS
      //-----------------------------------------


    public:
      //! destructor
      virtual ~bcsr_amg_matrix_iface ()
        {};

      //! build B = A^T from given csr matrix, return 0 if success,
      //! rows number and offset can be set manually with new_n_rows and rows_offset
      virtual int build_transpose (sp_bcsr_matrix_t matrix, const t_long rows_offset/* = 0*/,
                           const t_long cols_offset/* = 0*/, const t_long new_n_rows/* = 0*/) = 0;

      //! build B = A^T from given csr matrix using matrix structure only, without values, return 0 if success,
      //! rows number and offset can be set manually with new_n_rows and rows_offset
      virtual int build_transpose_struct (sp_bcsr_matrix_t matrix,
                                  const t_long rows_offset/* = 0*/,
                                  const t_long cols_offset/* = 0*/,
                                  const t_long new_n_rows/* = 0*/) = 0;
      /** 
       * @brief realize M = RAP (triple matrix product)
       * 
       * @param r_matrix        -- <INPUT> R matrix
       * @param a_matrix        -- <INPUT> A matrix
       * @param p_matrix        -- <INPUT> P matrix
       * @param update          -- <INPUT> if true than use existing matrix structure and update only values
       * 
       * @return 0 if success
       */
      virtual int triple_matrix_product (sp_bcsr_matrix_t r_matrix, sp_bcsr_matrix_t a_matrix, 
                                         sp_bcsr_matrix_t p_matrix, const bool update) = 0;

    };


}//namespace blue_sky
#endif //__BCSR_AMG_MATRIX_IFACE_H
