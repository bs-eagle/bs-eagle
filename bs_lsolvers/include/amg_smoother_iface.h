#ifndef __AMG_SMOOTHER_IFACE_H
#define __AMG_SMOOTHER_IFACE_H
/** 
 * @file amg_smoother_iface.h
 * @brief 
 * @date 2009-12-12
 */
#include <string>

#include "bs_assert.h"
#include "bs_tree.h"

#include "prop_iface.h"
#include "lsolver_iface.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "bcsr_matrix_iface.h"
#include "matrix_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()


namespace blue_sky
{
  /** 
   * @brief interface class for matrix storage and manipulation
   */
  template <class strategy_t>
  class amg_smoother_iface: public lsolver_iface<strategy_t>
    {
    public:
      //! matrix interface type
      typedef matrix_iface<strategy_t>                  matrix_t;
      //! internal fp type
      typedef typename strategy_t::fp_type_t            fp_type_t;
      //! internal fp type for matrix storage
      typedef typename strategy_t::fp_storage_type_t    fp_storage_type_t;
      //! internal integer type
      typedef typename strategy_t::i_type_t             i_type_t;
      //! this_t
      typedef lsolver_iface <strategy_t>                base_t;

      typedef bs_array<fp_type_t>                               fp_array_t;
      typedef bs_array<i_type_t>                                i_array_t;
      typedef bs_array<fp_storage_type_t>                       fp_storage_array_t;

      typedef smart_ptr<fp_array_t, true>                       sp_fp_array_t;
      typedef smart_ptr<i_array_t, true>                        sp_i_array_t;
      typedef smart_ptr<fp_storage_array_t, true>               sp_fp_storage_array_t;
      //
      //! bcsr_matrix
      typedef bcsr_matrix_iface<strategy_t> bcsr_t;
      //! prop 
      typedef prop_iface<float, int, std::string, bool> prop_t;
      //! short name to smart pointer to this class
      typedef smart_ptr<base_t, true>                   sp_base_t;              
      //! short name to smart pointer to properties holder class
      typedef smart_ptr<prop_t, true>                   sp_prop_t;              
      typedef smart_ptr<matrix_t, true>                 sp_matrix_t;
      typedef smart_ptr<bcsr_t, true>                   sp_bcsr_t;

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
      virtual int smooth (sp_bcsr_t matrix, sp_i_array_t coarse, const i_type_t iter_number, 
                          sp_fp_array_t rhs, sp_fp_array_t sol) = 0;
      
    public:
      //! destructor
      virtual ~amg_smoother_iface ()
        {}

    };

}//namespace blue_sky
#endif //__AMG_SMOOTHER_IFACE_H
