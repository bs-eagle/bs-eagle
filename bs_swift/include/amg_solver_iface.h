/** 
 * @file amg_solver_iface.h
 * @brief Algebraic Multi Grid linear system solver interface file
 * @author 
 * @date 2009-12-07
 */
#ifndef __AMG_SOLVER_IFACE_H
#define __AMG_SOLVER_IFACE_H


#include "bs_assert.h"
#include "bs_tree.h"
#include "conf.h"
#include BS_FORCE_PLUGIN_IMPORT ()
#include "lsolver_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
{
  /** 
   * @brief interface class for matrix storage and manipulation
   */
  class amg_solver_iface: public lsolver_iface
    {
    public:
      typedef amg_solver_iface                                  this_t;
      //! short name to smart pointer to this class
      typedef smart_ptr<this_t, true>                           sp_this_t;              

      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
      
      //virtual const bcsr_t &get_A (int layer) const = 0;

      //virtual const bcsr_t &get_P (int layer) const = 0;
      
      //virtual const fp_vector_type_t &get_rhs (int layer)

      
    public:

    public:
      //! destructor
      virtual ~amg_solver_iface ()
        {}

    };
};
#endif //__AMG_SOLVER_IFACE_H
