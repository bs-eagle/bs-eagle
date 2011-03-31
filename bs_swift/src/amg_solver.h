#ifndef __AMG_SOLVER_H
#define __AMG_SOLVER_H

/**
 * @file amg_solver.h
 * @brief class for implementing AMG solver
 * @author
 * @version
 * @date 2011-03-30
 */

#include "bcsr_matrix_iface.h"
#include "amg_solver_iface.h"

namespace blue_sky
  {

  class BS_API_PLUGIN amg_solver : public amg_solver_iface
    {

      //-----------------------------------------
      // TYPES
      //-----------------------------------------

    public:
      //! matrix interface type
      typedef matrix_iface                              matrix_t;
      typedef bcsr_matrix_iface                         bcsr_matrix_t;

      //! prop
      typedef prop_iface                                prop_t;
      typedef lsolver_iface                             base_t;         ///< typedef to this type. in child classes used as a short name of base class
      typedef smart_ptr<base_t, true>                   sp_base_t;      ///< short name to smart pointer to this class
      typedef smart_ptr<prop_t, true>                   sp_prop_t;      ///< short name to smart pointer to properties holder class

      typedef smart_ptr<matrix_t, true>                 sp_matrix_t;         ///< short name to smart pointer on matrix_t
      typedef smart_ptr<bcsr_matrix_t, true>            sp_bcsr_matrix_t;    ///< short name to smart pointer on matrix_t
      typedef std::vector<sp_bcsr_matrix_t>             sp_bcsr_matrix_t_list;
      //-----------------------------------------
      //  METHODS
      //-----------------------------------------

    public:
      // destructor
      virtual ~amg_solver ();

      // solve
      virtual int solve (sp_matrix_t matrix, spv_double rhs, spv_double sol);

      virtual int solve_prec (sp_matrix_t matrix, spv_double rhs, spv_double sol);

      // setup
      virtual int setup (sp_matrix_t matrix);

      //! set preconditioner
      virtual void set_prec (sp_base_t prec_)
        {
          //prec = prec_;
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

          s << "AMG\n";
          //if (prec)
          //  {
          //    s << "preconditioned by:\n" << prec->py_str ();
          //  }
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
      spv_long          aver_cop;     //!< average operator complexity

      sp_bcsr_matrix_t_list A;        //!< coarse level matrices vector
      sp_bcsr_matrix_t_list S;        //!< coarse level strength matrices vector
      sp_bcsr_matrix_t_list P;        //!< coarse level prolongation matrices vector

      BLUE_SKY_TYPE_DECL (amg_solver);
    };

  }	// namespace blue_sky

#endif //__AMG_SOLVER_H

