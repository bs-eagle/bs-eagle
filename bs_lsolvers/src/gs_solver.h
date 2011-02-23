#ifndef __GS_SOLVER_H
#define __GS_SOLVER_H
/** 
 * @file gs_solver.h
 * @brief GS linear iterative solver
 * @date 2009-12-16
 */
#include <string>
#include <sstream>

#include <amg_smoother_iface.h>
#include BS_FORCE_PLUGIN_IMPORT ()
#include "strategies.h"
#include "matrix_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {
  /**
  * @brief GS linear solver
  */
  
  class BS_API_PLUGIN gs_solver : public amg_smoother_iface
    {
      //-----------------------------------------
      // TYPES
      //-----------------------------------------
    public:
      //! matrix interface type
      typedef matrix_iface                              matrix_t;
      typedef lsolver_iface                             base_t;

      //! bcsr_matrix
      typedef bcsr_matrix_iface                         bcsr_t;
      typedef smart_ptr<bcsr_t, true>                   sp_bcsr_t;
      //! prop 
      typedef prop_iface<float, int, std::string, bool> prop_t;
      //! short name to smart pointer to this class
      typedef smart_ptr<base_t, true>                   sp_base_t;              
      //! short name to smart pointer to properties holder class
      typedef smart_ptr<prop_t, true>                   sp_prop_t;              

      typedef smart_ptr<matrix_t, true>                 sp_matrix_t;    ///< short name to smart pointer on matrix_t


      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
    public:
      // destructor
      virtual ~gs_solver ();

      // solve
      virtual int solve (sp_matrix_t matrix, spv_double rhs, spv_double sol);

      virtual int solve_prec (sp_matrix_t matrix, spv_double rhs, spv_double sol);

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
                          spv_double rhs, spv_double sol);
      // setup
      virtual int setup (sp_matrix_t matrix);

      //! set preconditioner
      virtual void set_prec (sp_base_t /*prec_*/)
        {
        }
      
      virtual void set_prop (sp_prop_t prop_);

      //! get properties
      virtual sp_prop_t get_prop() 
        {
          return prop;
        }

      //! return final residual
      virtual t_double get_final_residual () const
        {
          return prop->get_f (final_res_idx);
        }

      //! return number of used iterations
      virtual int get_niters () const
        {
          return prop->get_i (iters_idx);
        }

#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const
        {
          std::stringstream s;

          s << "GS linear solver.\n";
          s << "Properties:\n";
          s << prop->py_str ();

          return s.str ();
        }
#endif //BSPY_EXPORTING_PLUGIN

    protected:
      void init_prop ();
      //-----------------------------------------
      //  VARIABLES
      //-----------------------------------------
    public:
    protected:
      sp_prop_t         prop;         //!< properties for solvers
      spv_double        sp_r;

      int               tol_idx;
      int               max_iters_idx;
      int               final_res_idx;
      int               iters_idx;
      int               success_idx;
      int               invers_idx;
      int               cf_type_idx;

    public:
      BLUE_SKY_TYPE_DECL (gs_solver);
    };

  }	// namespace blue_sky

#endif //__GS_SOLVER_H

