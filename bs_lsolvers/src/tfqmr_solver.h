/**
 * \file tfqmr_solver.h
 * \brief TFQMR linear solver
 * \author SalimgareevaEM
 * \date
 * */
#ifndef BS_TFQMR_LINEAR_SOLVER_H_
#define BS_TFQMR_LINEAR_SOLVER_H_

#include <string>
#include <sstream>

#include "lsolver_iface.h"
#include BS_FORCE_PLUGIN_IMPORT ()
#include "strategies.h"
#include "matrix_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()


namespace blue_sky
  {
  /**
  * @brief CGS linear solver
  */
  template <class strategy_t>
  class BS_API_PLUGIN tfqmr_solver : public lsolver_iface<strategy_t>
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

      typedef bs_array<fp_type_t>                               fp_array_t;
      typedef bs_array<i_type_t>                                i_array_t;
      typedef bs_array<fp_storage_type_t>                       fp_storage_array_t;

      typedef smart_ptr<fp_array_t, true>                       sp_fp_array_t;
      typedef smart_ptr<i_array_t, true>                        sp_i_array_t;
      typedef smart_ptr<fp_storage_array_t, true>               sp_fp_storage_array_t;
      //! prop 
      typedef prop_iface<float, int, std::string, bool> prop_t;
      //typedef linear_solver_base<strategy_t>      this_t;         ///< typedef to this type
      typedef lsolver_iface<strategy_t>                 base_t;         ///< typedef to this type. in child classes used as a short name of base class
      typedef smart_ptr<base_t, true>                   sp_base_t;      ///< short name to smart pointer to this class
      typedef smart_ptr<prop_t, true>               sp_prop_t;      ///< short name to smart pointer to properties holder class

      typedef smart_ptr<matrix_t, true>           sp_matrix_t;    ///< short name to smart pointer on matrix_t


      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
    public:
      // destructor
      virtual ~tfqmr_solver ();

      // solve
      virtual int solve (sp_matrix_t matrix, sp_fp_array_t rhs, sp_fp_array_t sol);

      virtual int solve_prec (sp_matrix_t matrix, sp_fp_array_t rhs, sp_fp_array_t sol);

      // setup
      virtual int setup (sp_matrix_t matrix);

      //! set preconditioner
      virtual void set_prec (sp_base_t prec_)
        {
          prec = prec_;
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

          s << "TFQMR linear solver.\n";
          if (prec)
            {
              s << "preconditioned by:\n" << prec->py_str ();
            }
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
      sp_base_t         prec;         //!< pointer to the preconditioner
      sp_prop_t         prop;         //!< properties for solvers
      sp_fp_array_t     sp_p;
      sp_fp_array_t     sp_v;
      sp_fp_array_t     sp_w;
      sp_fp_array_t     sp_u;
      sp_fp_array_t     sp_q;
      sp_fp_array_t     sp_d;
      sp_fp_array_t     sp_res;
      sp_fp_array_t     sp_r;
      sp_fp_array_t     sp_rtilde;
      sp_fp_array_t     sp_tmp;
      sp_fp_array_t     sp_rhat;
      sp_fp_array_t     sp_y;
      int               tol_idx;
      int               max_iters_idx;
      int               final_res_idx;
      int               iters_idx;
      int               success_idx;

    public:
      BLUE_SKY_TYPE_DECL (tfqmr_solver);
    };

  }	// namespace blue_sky

#endif // #ifndef BS_CGS_LINEAR_SOLVER_H_

