/** 
 * @file two_stage_prec_iface.h
 * @brief 
 * @date 2009-09-01
 */
#ifndef TWO_STAGE_PREC_IFACE_H__
#define TWO_STAGE_PREC_IFACE_H__

#include <string>

#include "bs_assert.h"
#include "bs_tree.h"

#include "prop_iface.h"
#include "lsolver_iface.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "matrix_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()


namespace blue_sky
{
  /** 
   * @brief interface class for matrix storage and manipulation
   */
  template <class strategy_t>
  class two_stage_prec_iface: public lsolver_iface<strategy_t>
    {
    public:
      //! matrix interface type
      typedef matrix_iface<strategy_t>                  matrix_t;
      //! internal fp type
      typedef typename strategy_t::fp_type_t            fp_type_t;
      //! internal integer type
      typedef typename strategy_t::i_type_t             i_type_t;
      //! this_t
      typedef lsolver_iface <strategy_t>                base_t;
      //! prop 
      typedef prop_iface<float, int, std::string, bool> prop_t;
      //! short name to smart pointer to this class
      typedef smart_ptr<base_t, true>                   sp_base_t;              
      //! short name to smart pointer to properties holder class
      typedef smart_ptr<prop_t, true>               sp_prop_t;              

      typedef bs_array<fp_type_t>                               fp_array_t;
      typedef bs_array<i_type_t>                                i_array_t;
      //typedef bs_array<fp_storage_type_t>                       fp_storage_array_t;

      typedef smart_ptr<fp_array_t, true>                       sp_fp_array_t;
      typedef smart_ptr<i_array_t, true>                        sp_i_array_t;
      //typedef smart_ptr<fp_storage_array_t, true>               sp_fp_storage_array_t;

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
