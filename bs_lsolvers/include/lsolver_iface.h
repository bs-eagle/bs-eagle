/** 
 * @file lsolver_iface.h
 * @brief 
 * @author Oleg Borschuk
 * @date 2009-09-01
 */
#ifndef __LSOLVER_IFACE_H
#define __LSOLVER_IFACE_H

#include <string>

#include "bs_assert.h"
#include "bs_tree.h"
#include "bs_array.h"

#include "prop_iface.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "matrix_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()



namespace blue_sky
{
  /** 
   * @brief interface class for linear system solver and preconditioners 
   *    (Ax = b), where A -- matrix, x and b -- vectors,
   *            Calling order:
   *                    1) setup properties of the solver using get_prop () and set_prop ()
   *                    2) setup preconditioner by calling set_prec ()
   *                    3) call setup () -- for prepare solver and preconditioner
   *                    4) call solve () -- for solving linear system (for sequence of linear systems A(x_i) = b_i, 1--3 should be call only once)  
   */
  template <class strategy_t>
  class lsolver_iface: public bs_node
    {
    public:
      //! matrix interface type
      typedef matrix_iface<strategy_t>                  matrix_t;
      //! internal fp type
      typedef typename strategy_t::fp_type_t            fp_type_t;
      //! internal integer type
      typedef typename strategy_t::i_type_t             i_type_t;
      
      typedef bs_array<fp_type_t>                               fp_array_t;
      typedef bs_array<i_type_t>                                i_array_t;
      //typedef bs_array<fp_storage_type_t>                       fp_storage_array_t;

      typedef smart_ptr<fp_array_t, true>                       sp_fp_array_t;
      typedef smart_ptr<i_array_t, true>                        sp_i_array_t;
      //typedef smart_ptr<fp_storage_array_t, true>               sp_fp_storage_array_t;

      //! this_t
      typedef lsolver_iface <strategy_t>                this_t;
      //! prop 
      typedef prop_iface<float, int, std::string, bool> prop_t;
      //! short name to smart pointer to this class
      typedef smart_ptr<this_t, true>                   sp_this_t;              
      //! short name to smart pointer to properties holder class
      typedef smart_ptr<prop_t, true>                   sp_prop_t;              
      typedef smart_ptr<matrix_t, true>                 sp_matrix_t;

      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
      
    public:
      /** 
       * @brief solve linear system (Ax=b), where A -- <matrix>, x -- <sol>, b -- <rhs> 
       * 
       * @param matrix  -- <INPUT> matrix to solve
       * @param rhs     -- <INPUT> right hand side
       * @param sol     -- <INPUT/OUTPUT> solution
       * 
       * @return 0 if success
       */
      virtual int solve (sp_matrix_t matrix, sp_fp_array_t rhs, sp_fp_array_t sol) = 0;

      /** 
       * @brief do only one solver iteration (Ax=b), where A -- <matrix>, x -- <sol>, b -- <rhs> 
       *                this method called if linear solver used as a preconditioner 
       * 
       * @param matrix  -- <INPUT> matrix to solve
       * @param rhs     -- <INPUT> right hand side
       * @param sol     -- <INPUT/OUTPUT> solution
       * 
       * @return 0 if success
       */
      virtual int solve_prec (sp_matrix_t matrix, sp_fp_array_t rhs, sp_fp_array_t sol) = 0;

      
      /** 
       * @brief prepare linear solver and preconditioner for solution
       * 
       * @param matrix -- <INPUT/OUTPUT> 
       * 
       * @return 0 if success
       */
      virtual int setup (sp_matrix_t matrix) = 0;

      /** 
       * @brief set preconditioner
       * 
       * @param prec -- <INPUT/OUTPUT> -- smart pointer to the preconditioner
       */
      virtual void set_prec (sp_this_t prec) = 0;
      
      /** 
       * @brief set solver's properties
       * 
       * @param prop -- <INPUT/OUTPUT> smart pointer to the properties class
       */
      virtual void set_prop(sp_prop_t prop) = 0;

      /** 
       * @brief return smart pointer to the linear solver properties
       * 
       * @return smart pointer 
       */
      virtual sp_prop_t get_prop() = 0;

      /** 
       * @brief for iterative linear solvers return final residual, for exact methods 0
       * 
       * @return final residual
       */
      virtual fp_type_t get_final_residual () const = 0;

      //! return number of used iterations
      /** 
       * @brief for iterative linear solvers return number of iteration used to get solution 
       *              or maximum allowed number of iterations, for exact methods 1
       * 
       * @return number of iterations
       */
      virtual int get_niters () const = 0;

#ifdef BSPY_EXPORTING_PLUGIN
      /** 
       * @brief used for python wraper (method __str__ ())
       * 
       * @return string
       */
      virtual std::string py_str () const = 0;
#endif //BSPY_EXPORTING_PLUGIN

    public:
      //! destructor
      virtual ~lsolver_iface ()
        {}

    };

}//namespace blue_sky

#endif //__LSOLVER_IFACE_H
