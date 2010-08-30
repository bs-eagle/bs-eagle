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
  template <class strategy_t>
  class BS_API_PLUGIN gs_solver : public amg_smoother_iface<strategy_t>
    {
      //-----------------------------------------
      // TYPES
      //-----------------------------------------
    public:
      //! matrix interface type
      typedef matrix_iface<strategy_t>                  matrix_t;
      //! internal fp type
      typedef typename strategy_t::fp_type_t            fp_type_t;
      typedef typename strategy_t::fp_storage_type_t    fp_storage_type_t;
      //! internal integer type
      typedef typename strategy_t::i_type_t             i_type_t;
      //! this_t
      typedef lsolver_iface <strategy_t>                        base_t;

      typedef bs_array<fp_type_t>                               fp_array_t;
      typedef bs_array<i_type_t>                                i_array_t;
      typedef bs_array<fp_storage_type_t>                       fp_storage_array_t;

      typedef smart_ptr<fp_array_t, true>                       sp_fp_array_t;
      typedef smart_ptr<i_array_t, true>                        sp_i_array_t;
      typedef smart_ptr<fp_storage_array_t, true>               sp_fp_storage_array_t;

      //! bcsr_matrix
      typedef bcsr_matrix_iface<strategy_t>                     bcsr_t;
      typedef smart_ptr<bcsr_t, true>                           sp_bcsr_t;
      //! prop 
      typedef prop_iface<float, int, std::string, bool> prop_t;
      //! short name to smart pointer to this class
      typedef smart_ptr<base_t, true>                   sp_base_t;              
      //! short name to smart pointer to properties holder class
      typedef smart_ptr<prop_t, true>               sp_prop_t;              

      typedef smart_ptr<matrix_t, true>                 sp_matrix_t;    ///< short name to smart pointer on matrix_t


      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
    public:
      // destructor
      virtual ~gs_solver ();

      // solve
      virtual int solve (sp_matrix_t matrix, sp_fp_array_t rhs, sp_fp_array_t sol);

      virtual int solve_prec (sp_matrix_t matrix, sp_fp_array_t rhs, sp_fp_array_t sol);

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
                          sp_fp_array_t rhs, sp_fp_array_t sol);
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
      virtual fp_type_t get_final_residual () const
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
      sp_fp_array_t     sp_r;

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

