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
  const std::string strength_threshold_idx = "strength_threshold_idx";
  const std::string max_row_sum_idx = "max_row_sum";
  const std::string n_last_level_points_idx = "n_last_level_points";

  class BS_API_PLUGIN amg_solver : public amg_solver_iface
    {

      //-----------------------------------------
      // TYPES
      //-----------------------------------------

    public:
      //! matrix interface type
      typedef matrix_iface                              matrix_t;
      typedef bcsr_amg_matrix_iface                     bcsr_t;
      typedef std::vector<spv_long>                     vec_spv_long;   ///< vector of smart pointers to vector<long>

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
      typedef std::vector<sp_smbuild_t>                 vec_sp_smbuild_t;
      typedef std::vector<sp_coarse_t>                  vec_sp_coarse_t;
      typedef std::vector<sp_pbuild_t>                  vec_sp_pbuild_t;

      //!
      typedef smart_ptr<matrix_t, true>                 sp_matrix_t;    ///< short name to smart pointer on matrix_t
      typedef smart_ptr<bcsr_t, true>                   sp_bcsr_t;      ///< short name to smart pointer on matrix_t
      typedef std::vector<sp_bcsr_t>                    vec_sp_bcsr_t;  ///< vector of smart pointers to bcsr_amg_matrix_iface
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

      //! return strength_threshold
      virtual t_double get_strength_threshold () const
        {
          return prop->get_f (strength_threshold_idx);
        }

      //! return max_row_sum
      virtual t_double get_max_row_sum () const
        {
          return prop->get_f (max_row_sum_idx);
        }

      //! return n_last_level_points
      virtual t_long get_n_last_level_points() const
        {
          return prop->get_i (n_last_level_points_idx);
        }

      //! set strength matrix builder for amg level
      void set_smbuilder (unsigned int level, sp_smbuild_t sp_smbuilder_iface)
        {
          unsigned int size = smbuilder.size ();
          //use last smbuilder [size-1]) for [size],..,[level-1]
          if (size <= level)
            {
              smbuilder.resize (level + 1);
              for (unsigned int i = size; i < level; i++)
                {
                  smbuilder[i] = smbuilder[size - 1];
                }
            }
          smbuilder[level] = sp_smbuilder_iface;
        }

      //! set coarse type for amg level
      void set_coarser (unsigned int level, sp_coarse_t sp_coarse_iface)
        {
          unsigned int size = coarser.size ();
          //use last coarser [size-1]) for [size],..,[level-1]
          if (size <= level)
            {
              coarser.resize (level + 1);
              for (unsigned int i = size; i < level; i++)
                {
                  coarser[i] = coarser[size - 1];
                }
            }
          coarser[level] = sp_coarse_iface;
        }

      //! set interpolation matrix builder for amg level
      void set_pbuilder (unsigned int level, sp_pbuild_t sp_pbuilder_iface)
        {
          unsigned int size = pbuilder.size ();
          //use last pbuilder [size-1]) for [size],..,[level-1]
          if (size <= level)
            {
              pbuilder.resize (level + 1);
              for (unsigned int i = size; i < level; i++)
                {
                  pbuilder[i] = pbuilder[size - 1];
                }
            }
          pbuilder[level] = sp_pbuilder_iface;
        }

      sp_smbuild_t get_smbuilder (unsigned int level)
        {
          return (level < smbuilder.size ()) ? smbuilder[level] : smbuilder[smbuilder.size () - 1];
        }

      sp_coarse_t get_coarser (unsigned int level)
        {
          return (level < coarser.size ()) ? coarser[level] : coarser[coarser.size () - 1];
        }

      sp_pbuild_t get_pbuilder (unsigned int level)
        {
          return (level < pbuilder.size ()) ? pbuilder[level] : pbuilder[pbuilder.size () - 1];
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
      sp_base_t             lu_solver;    //!< LU solver for last level

      vec_spv_long          s;            //!< S markers
      vec_spv_long          cf;           //!< CF markers

      vec_sp_smbuild_t      smbuilder;    //!<
      vec_sp_coarse_t       coarser;      //!<
      vec_sp_pbuild_t       pbuilder;     //!<

      vec_sp_bcsr_t a;                    //!< coarse level matrices vector
      vec_sp_bcsr_t p;                    //!< coarse level prolongation matrices vector

      BLUE_SKY_TYPE_DECL (amg_solver);
    };

  }	// namespace blue_sky

#endif //__AMG_SOLVER_H

