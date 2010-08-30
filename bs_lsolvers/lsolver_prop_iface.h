/** 
 * @file l_solver_prop_iface.h
 * @brief 
 * @author Oleg Borschuk
 * @date 2009-08-21
 */
#ifndef __LSOLVER_PROP_IFACE_H
#define __LSOLVER_PROP_IFACE_H

#include BS_FORCE_PLUGIN_IMPORT ()
#include "matrix_iface.h"
#include "named_pbase_access.h"
#include "bos_report.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {

  /**
  * \brief properties for linear solvers
  */
  class BS_API_PLUGIN lsolver_prop_iface 
    {
    public:
      typedef double    fp_type_t;
      typedef int       i_type_t;

      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
    public:
      //! destructor
      virtual ~linear_solver_prop ()
        {}

      //! set maximum number of iterations
      virtual int set_max_iters (int /*n_iters*/) = 0;

      //! set tolerance
      virtual void set_tolerance (fp_type_t /*new_tol*/) = 0;

      //! set tolerance
      virtual void set_matbal_tolerance (fp_type_t /*new_tol*/) = 0;

      //! set number of iters
      virtual void set_iters (int /*n_iters*/) = 0;

      //! set successively converged flag
      virtual void set_success (int /*success*/) = 0;

      //! set final resid
      virtual void set_final_resid (fp_type_t /*final_resid*/) = 0;

      //! set relative factor
      virtual void set_relative_factor (fp_type_t /*relative_factor*/) = 0;

      //! return != 0 if method successfully converged
      virtual int check_convergence () const = 0;

      //! return number of iteration
      virtual int get_iters () const = 0;

      //! return relative residual denominator
      virtual fp_type_t get_relative_factor () const = 0;

      //! return maximum allowed number of iterations
      virtual int get_max_iters () const = 0;

      //! return tolerance
      virtual fp_type_t get_tolerance () const = 0;

      //! get successively converged flag
      virtual int get_success () const = 0;

      //! get final resid
      virtual fp_type_t get_final_resid () const = 0;

      //-----------------------------------------
      //  VARIABLES
      //-----------------------------------------
    public:

    public:
      virtual const std::string &get_params_name (i_type_t idx) = 0;
    };
  }
#endif //__LSOLVER_PROP_IFACE_H

