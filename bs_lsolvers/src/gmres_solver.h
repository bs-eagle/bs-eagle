/**
 * \file gmres_solver.h
 * \brief GMRES linear solver
 * \author SalimgareevaEM
 * \date
 * */
#ifndef BS_GMRES_LINEAR_SOLVER_H_
#define BS_GMRES_LINEAR_SOLVER_H_

#include <string>
#include <sstream>

#include "lsolver_iface.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "matrix_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {
  const std::string m_idx = "m_size";
  /**
  * @brief GMRES linear solver
  */
  class BS_API_PLUGIN gmres_solver : public lsolver_iface
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
      virtual ~gmres_solver ();


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

          s << "GMRES linear solver.\n";
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
      std::vector<spv_double> vec_p;
      spv_double        sp_w;
      spv_double        sp_r;
      spv_double        sp_s;
      spv_double        sp_c;
      spv_double        sp_rs;
      spv_double        sp_hh;


    public:
      BLUE_SKY_TYPE_DECL (gmres_solver);
    };

  }	// namespace blue_sky

#endif // #ifndef BS_GMRES_LINEAR_SOLVER_H_

