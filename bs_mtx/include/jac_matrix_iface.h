/** 
 * @file jac_matrix_iface.h
 * @brief BS interface for jacobian matrix (consist from acc_matrix (bdiag_matrix), flux_matrix (bcsr_matrix_iface))
 * @author Oleg Borschuk
 * @date 2009-07-30
 */
#ifndef __JAC_MATRIX_IFACE_H
#define __JAC_MATRIX_IFACE_H

#include "bcsr_matrix_iface.h"
#include "bdiag_matrix_iface.h"

namespace blue_sky
{
  class jac_matrix_iface : public matrix_iface
    {
      //////////////////////////////////////////////////////////////////////////
      // TYPES
      //////////////////////////////////////////////////////////////////////////
      typedef bcsr_matrix_iface                                 bcsr_matrix_iface_t;
      typedef bdiag_matrix_iface                                bdiag_matrix_iface_t;
      typedef jac_matrix_iface                                  this_t;
      typedef matrix_iface                                      matrix_iface_t;

      typedef smart_ptr<bcsr_matrix_iface_t, true>              sp_bcsr_matrix_iface_t;
      typedef smart_ptr<bdiag_matrix_iface_t, true>             sp_bdiag_matrix_iface_t;
    public:

      //////////////////////////////////////////////////////////////////////////
      // METHODS
      //////////////////////////////////////////////////////////////////////////
    public:

      //! return flux (regular) matrix
      virtual sp_bcsr_matrix_iface_t get_flux_matrix () = 0;

      //! return facility (irregular) matrix
      virtual sp_bcsr_matrix_iface_t get_facility_matrix () = 0;

      //! return accumulative matrix
      virtual sp_bdiag_matrix_iface_t get_accum_matrix () = 0;

    };

} // namespace blue_sky

#endif //__JAC_MATRIX_IFACE_H
