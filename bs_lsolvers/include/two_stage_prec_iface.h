/** 
 * @file two_stage_prec_iface.h
 * @brief 
 * @date 2009-09-01
 */
#ifndef TWO_STAGE_PREC_IFACE_H__
#define TWO_STAGE_PREC_IFACE_H__

#include "bs_assert.h"
#include "bs_tree.h"

#include "prop_iface.h"
#include "lsolver_iface.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "matrix_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include <string>

namespace blue_sky
{
  /** 
   * @brief interface class for matrix storage and manipulation
   */
  
  class BS_API_PLUGIN two_stage_prec_iface: public lsolver_iface
    {
    public:
      //! matrix interface type
      typedef matrix_iface                                      matrix_t;
      //! this_t
      typedef lsolver_iface                                     base_t;
      //! prop 
      typedef prop_iface                                        prop_t;
      //! short name to smart pointer to this class
      typedef smart_ptr<base_t, true>                           sp_base_t;              
      //! short name to smart pointer to properties holder class
      typedef smart_ptr<prop_t, true>                           sp_prop_t;              

      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
      
    public:
      //! set preconditioner
      virtual void set_prec_2 (sp_base_t prec) = 0;
      
    public:
      //! destructor
      virtual ~two_stage_prec_iface ()
        {}

    };

}//namespace blue_sky

#endif //TWO_STAGE_PREC_IFACE_H__
