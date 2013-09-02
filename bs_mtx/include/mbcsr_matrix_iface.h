/**
 * @file mbcsr_matrix_iface.h
 * @brief Block CSR multi matrix interface class
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-02-27
 */
#ifndef MBCSR_MATRIX_IFACE_DT4R0T39

#define MBCSR_MATRIX_IFACE_DT4R0T39

#include "bcsr_matrix_iface.h"
#include <string>

namespace blue_sky
{
  class mbcsr_matrix_iface : public matrix_iface
    {
      //////////////////////////////////////////////////////////////////////////
      // TYPES
      //////////////////////////////////////////////////////////////////////////
      typedef bcsr_matrix_iface                                 bcsr_iface_t;
      typedef mbcsr_matrix_iface                                this_t;
      typedef matrix_iface                                      matrix_iface_t;

      typedef smart_ptr<bcsr_iface_t, true>                     sp_bcsr_iface_t;
    public:

      //////////////////////////////////////////////////////////////////////////
      // METHODS
      //////////////////////////////////////////////////////////////////////////
    public:

      //! return smart pointer to the matrix with name
      virtual sp_bcsr_iface_t get_matrix (const std::string &name) = 0;

      //! add new BCSR matrix to the list
      virtual void add_matrix (const std::string &name, sp_bcsr_iface_t m) = 0;

      //! clear all
      virtual void clear () = 0;


      /**
       * @brief merge all matrixies together and apply filter
       *        filter is a array length n_rows with elements 0 or 1
       *        if 0 in the i position than merged matrix i-th row should include
       *        only diagonal element, the same for the i-th column
       *
       * @param filter -- <INPUT> array length n_rows, elements should be 0 or 1
       *
       * @return merged BCSR matrix
       */
      virtual sp_bcsr_iface_t merge (spv_int filter) = 0;
    };

} // namespace blue_sky


#endif /* end of include guard: MBCSR_MATRIX_IFACE_DT4R0T39 */
