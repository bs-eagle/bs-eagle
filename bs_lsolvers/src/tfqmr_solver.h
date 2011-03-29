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
  
  class BS_API_PLUGIN tfqmr_solver : public lsolver_iface
    {

      //-----------------------------------------
      // TYPES
      //-----------------------------------------
    public:
      //! matrix interface type
      typedef matrix_iface                              matrix_t;
      //! prop 
      typedef prop_iface                                prop_t;
      typedef lsolver_iface                             base_t;         ///< typedef to this type. in child classes used as a short name of base class
      typedef smart_ptr<base_t, true>                   sp_base_t;      ///< short name to smart pointer to this class
      typedef smart_ptr<prop_t, true>                   sp_prop_t;      ///< short name to smart pointer to properties holder class

      typedef smart_ptr<matrix_t, true>                 sp_matrix_t;    ///< short name to smart pointer on matrix_t


      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
    public:
      // destructor
      virtual ~tfqmr_solver ();

      // solve
      virtual int solve (sp_matrix_t matrix, spv_double rhs, spv_double sol);

      virtual int solve_prec (sp_matrix_t matrix, spv_double rhs, spv_double sol);

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
      spv_double        sp_p;
      spv_double        sp_v;
      spv_double        sp_w;
      spv_double        sp_u;
      spv_double        sp_q;
      spv_double        sp_d;
      spv_double        sp_res;
      spv_double        sp_r;
      spv_double        sp_rtilde;
      spv_double        sp_tmp;
      spv_double        sp_rhat;
      spv_double        sp_y;

    public:
      BLUE_SKY_TYPE_DECL (tfqmr_solver);
    };

  }	// namespace blue_sky

#endif // #ifndef BS_CGS_LINEAR_SOLVER_H_

