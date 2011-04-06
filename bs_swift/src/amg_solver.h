#ifndef __AMG_SOLVER_H
#define __AMG_SOLVER_H

/**
 * @file amg_solver.h
 * @brief class for implementing AMG solver
 * @author
 * @version
 * @date 2011-03-30
 */

#include "bcsr_amg_matrix_iface.h"
#include "amg_solver_iface.h"
#include "amg_smbuilder_iface.h"
#include "amg_coarse_iface.h"
#include "amg_pbuild_iface.h"

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
      typedef bcsr_amg_matrix_iface                     bcsr_t;
      typedef std::vector<spv_long>                     spv_long_vec;   ///< vector of smart pointers to vector<long>

      //! prop
      typedef prop_iface                                prop_t;
      typedef lsolver_iface                             base_t;         ///< typedef to this type. in child classes used as a short name of base class
      typedef smart_ptr<base_t, true>                   sp_base_t;      ///< short name to smart pointer to this class
      typedef smart_ptr<prop_t, true>                   sp_prop_t;      ///< short name to smart pointer to properties holder class

      typedef amg_smbuilder_iface                       smbuild_t;
      typedef amg_coarse_iface                          coarse_t;
      typedef amg_pbuild_iface                          pbuild_t;
      typedef smart_ptr<smbuild_t, true>                sp_smbuild_t;
      typedef smart_ptr<coarse_t, true>                 sp_coarse_t;
      typedef smart_ptr<pbuild_t, true>                 sp_pbuild_t;

      //!
      typedef smart_ptr<matrix_t, true>                 sp_matrix_t;    ///< short name to smart pointer on matrix_t
      typedef smart_ptr<bcsr_t, true>                   sp_bcsr_t;      ///< short name to smart pointer on matrix_t
      typedef std::vector<sp_bcsr_t>                    sp_bcsr_t_vec;  ///< vector of smart pointers to bcsr_amg_matrix_iface
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
      virtual void set_prec (sp_base_t /*prec_*/)
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

      void set_strength_type (int level, int type)
        {
          strength_type->resize (level);
          (*strength_type)[level] = type;
        }

      void set_coarse_type (int level, int type)
        {
          coarse_type->resize (level);
          (*coarse_type)[level] = type;
        }

      void set_interp_type (int level, int type)
        {
          interp_type->resize (level);
          (*interp_type)[level] = type;
        }

      int get_strength_type (int level)
        {
          if (strength_type->empty ())
            return 0;
          return (level < strength_type->size ()) ? (*strength_type)[level] : (*strength_type)[strength_type->size () - 1];
        }

      int get_coarse_type (int level)
        {
          if (coarse_type->empty ())
            return 0;
          return (level < coarse_type->size ()) ? (*coarse_type)[level] : (*coarse_type)[coarse_type->size () - 1];
        }

      int get_interp_type (int level)
        {
          if (interp_type->empty ())
            return 0;
          return (level < interp_type->size ()) ? (*interp_type)[level] : (*interp_type)[interp_type->size () - 1];
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
      sp_base_t             prec;         //!< pointer to the preconditioner
      sp_prop_t             prop;         //!< properties for solvers
      spv_long_vec          s;            //!< S markers
      spv_long_vec          cf;           //!< CF markers

      spv_int               strength_type;//!<
      spv_int               coarse_type;  //!<
      spv_int               interp_type;  //!<

      sp_bcsr_t_vec A;        //!< coarse level matrices vector
      sp_bcsr_t_vec P;        //!< coarse level prolongation matrices vector

      BLUE_SKY_TYPE_DECL (amg_solver);
    };

  }	// namespace blue_sky

#endif //__AMG_SOLVER_H

