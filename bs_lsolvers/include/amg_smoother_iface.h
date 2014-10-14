#ifndef __AMG_SMOOTHER_IFACE_H
#define __AMG_SMOOTHER_IFACE_H
/** 
 * @file amg_smoother_iface.h
 * @brief 
 * @date 2009-12-12
 */

#include "bs_tree.h"
#include "bs_assert.h"

#include "prop_iface.h"
#include "lsolver_iface.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "bcsr_matrix_iface.h"
#include "matrix_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include <string>

namespace blue_sky
{
  /** 
   * @brief interface class for matrix storage and manipulation
   */
  
  class BS_API_PLUGIN amg_smoother_iface: public lsolver_iface
    {
    public:
      //! matrix interface type
      typedef matrix_iface                                      matrix_t;
      //! this_t
      typedef lsolver_iface                                     base_t;

      //! bcsr_matrix
      typedef bcsr_matrix_iface                                 bcsr_t;
      //! prop 
      typedef prop_iface                                        prop_t;
      //! short name to smart pointer to this class
      typedef smart_ptr<base_t, true>                           sp_base_t;              
      //! short name to smart pointer to properties holder class
      typedef smart_ptr<prop_t, true>                           sp_prop_t;              
      typedef smart_ptr<matrix_t, true>                         sp_matrix_t;
      typedef smart_ptr<bcsr_t, true>                           sp_bcsr_t;

      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
      
    public:
      /** 
       * @brief smooth solution (make onlu one iteration)
       * 
       * @param matrix  -- <INPUT> Block CSR matrix
       * @param coarse  -- <INPUT> markers vector
       * @param rhs     -- <INPUT> Right hand side
       * @param sol     -- <INPUT/OUTPUT> solution
       * 
       * @return 0 if suceess
       */
      virtual int smooth (sp_bcsr_t matrix, spv_long coarse, const t_long iter_number, 
                          spv_double rhs, spv_double sol) = 0;
      
    public:
      //! destructor
      virtual ~amg_smoother_iface ()
        {}

    };

}//namespace blue_sky
#endif //__AMG_SMOOTHER_IFACE_H
