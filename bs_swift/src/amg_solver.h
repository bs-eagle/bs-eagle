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
#include "dens_matrix_iface.h"
#include "amg_solver_iface.h"
#include "amg_smbuilder_iface.h"
#include "amg_coarse_iface.h"
#include "amg_pbuild_iface.h"
#include "amg_smoother_iface.h"

namespace blue_sky
  {
  const std::string strength_threshold_idx = "strength_threshold";
  const std::string max_row_sum_idx = "max_row_sum";
  const std::string n_last_level_points_idx = "n_last_level_points";
  const std::string n_levels_idx = "n_levels";
  const std::string n_levels_max_idx = "n_levels_max";
  const std::string cop_idx = "Cop";
  const std::string n_pre_smooth_iters_idx = "n_pre_smooth_iters";
  const std::string n_post_smooth_iters_idx = "n_post_smooth_iters";

  class BS_API_PLUGIN amg_solver : public amg_solver_iface
    {

      //-----------------------------------------
      // TYPES
      //-----------------------------------------

    public:
      //! matrix interface type
      typedef matrix_iface                              matrix_t;
      typedef bcsr_amg_matrix_iface                     bcsr_t;
      typedef dens_matrix_iface                         dens_matrix_t;
      typedef std::vector<spv_long>                     vec_spv_long;   ///< vector of smart pointers to vector<long>
      typedef std::vector<spv_double>                   vec_spv_double; ///< vector of smart pointers to vector<double>

      //! prop
      typedef prop_iface                                prop_t;
      typedef lsolver_iface                             base_t;         ///< typedef to this type. in child classes used as a short name of base class
      typedef smart_ptr<base_t, true>                   sp_base_t;      ///< short name to smart pointer to this class
      typedef smart_ptr<prop_t, true>                   sp_prop_t;      ///< short name to smart pointer to properties holder class

      typedef smart_ptr<amg_smbuilder_iface, true>      sp_smbuild_t;
      typedef smart_ptr<amg_coarse_iface, true>         sp_coarse_t;
      typedef smart_ptr<amg_pbuild_iface, true>         sp_pbuild_t;
      typedef smart_ptr<amg_smoother_iface, true>       sp_smooth_t;
      typedef std::vector<sp_smbuild_t>                 vec_sp_smbuild_t;
      typedef std::vector<sp_coarse_t>                  vec_sp_coarse_t;
      typedef std::vector<sp_pbuild_t>                  vec_sp_pbuild_t;
      typedef std::vector<sp_smooth_t>                  vec_sp_smooth_t;

      //!
      typedef smart_ptr<matrix_t, true>                 sp_matrix_t;    ///< short name to smart pointer on matrix_t
      typedef smart_ptr<bcsr_t, true>                   sp_bcsr_t;      ///< short name to smart pointer on matrix_t
      typedef std::vector<sp_bcsr_t>                    vec_sp_bcsr_t;  ///< vector of smart pointers to bcsr_amg_matrix_iface
      typedef smart_ptr<dens_matrix_t, true>            sp_dens_matrix_t;    ///< short name to smart pointer on matrix_t
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

      //! get preconditioner
      virtual sp_base_t get_prec ()
        {
          return NULL;//prec;
        }

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

      //! return maximal number of AMG levels
      virtual int get_n_levels_max () const
        {
          return prop->get_i (n_levels_max_idx);
        }

      //! set strength matrix builder for amg level
      void set_smbuilder (unsigned int level, sp_smbuild_t sp_smbuilder_iface)
        {
          set_tool (smbuilder_vec, level, sp_smbuilder_iface);
        }

      //! set coarse type for amg level
      void set_coarser (unsigned int level, sp_coarse_t sp_coarse_iface)
        {
          set_tool (coarser_vec, level, sp_coarse_iface);
        }

      //! set interpolation matrix builder for amg level
      void set_pbuilder (unsigned int level, sp_pbuild_t sp_pbuilder_iface)
        {
          set_tool (pbuilder_vec, level, sp_pbuilder_iface);
        }

      //! set pre smoothing method for amg level
      void set_pre_smoother (unsigned int level, sp_smooth_t sp_smooth_iface)
        {
          set_tool (pre_smoother_vec, level, sp_smooth_iface);
        }

      //! set post smoothing method for amg level
      void set_post_smoother (unsigned int level, sp_smooth_t sp_smooth_iface)
        {
          set_tool (post_smoother_vec, level, sp_smooth_iface);
        }

      const vec_sp_bcsr_t get_matrices () const
        {
          return a;
        }

      const vec_sp_bcsr_t get_p_matrices () const
        {
          return p;
        }

      const vec_spv_long get_cf_markers () const
        {
          return cf;
        }

      const vec_spv_long get_s_markers () const
        {
          return s;
        }

      const vec_spv_double get_sol () const
        {
          return sol;
        }

      const vec_spv_double get_rhs () const
        {
          return rhs;
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

    private:
      //! store number of AMG setup levels to prop
      //! store maximal number of AMG levels to prop
      virtual void set_n_levels (t_long n_levels)
        {
          t_long n_levels_max = prop->get_i (n_levels_max_idx);
          if (n_levels > n_levels_max)
            prop->set_i (n_levels_max_idx, n_levels);
          prop->set_i (n_levels_idx, n_levels);
        }

      //! set tool for amg level
      //! set tools_vec[level] = new_tool, add (and set =last) elements if need
      template <class sp_tool_iface> inline
      void set_tool (std::vector<sp_tool_iface> tools_vec,
                     const unsigned int level,
                     const sp_tool_iface new_tool)
        {
          unsigned int size = tools_vec.size ();
          //use last tools_vec [size-1]) for [size],..,[level-1]
          if (size <= level)
            {
              tools_vec.resize (level + 1);
              for (unsigned int i = size; i < level; i++)
                {
                  tools_vec[i] = tools_vec[size - 1];
                }
            }
          tools_vec[level] = new_tool;
        }

      //! get tool for amg level
      //! return tools_vec[level] if exist,  else return last: tools_vec[size-1]
      template <class sp_tool_iface> inline
      sp_tool_iface get_tool (const std::vector<sp_tool_iface> tools_vec,
                              const unsigned int level) const
        {
          return (level < tools_vec.size ()) ? tools_vec[level] :
                  tools_vec[tools_vec.size () - 1];
        }

      //-----------------------------------------
      //  VARIABLES
      //-----------------------------------------
    public:
    protected:
      sp_base_t             prec;         //!< pointer to the preconditioner

    private:
      sp_prop_t             prop;         //!< properties for solvers
      sp_base_t             lu_solver;    //!< LU solver for last level
      sp_dens_matrix_t      lu_fact;      //!< dense matrix on last level
      spv_double            wksp;         //!< temporary vector
      sp_bcsr_t             r_matrix;     //!< temporary matrix for transposed p_matrix

      vec_sp_smbuild_t      smbuilder_vec;    //!< vector of smbuilders on each level
      vec_sp_coarse_t       coarser_vec;      //!< vector of coarsers on each level
      vec_sp_pbuild_t       pbuilder_vec;     //!< vector of pbuilders on each level
      vec_sp_smooth_t       pre_smoother_vec; //!< vector of pre_smoothers on each level
      vec_sp_smooth_t       post_smoother_vec;//!< vector of post_smoothers on each level

      vec_sp_bcsr_t         a;            //!< coarse level matrices vector
      vec_sp_bcsr_t         p;            //!< coarse level prolongation matrices vector
      vec_spv_long          s;            //!< S markers
      vec_spv_long          cf;           //!< CF markers
      vec_spv_double        rhs;          //!< coarse level rhs vector
      vec_spv_double        sol;          //!< coarse level sol vector

      BLUE_SKY_TYPE_DECL (amg_solver);
    };

  }	// namespace blue_sky

#endif //__AMG_SOLVER_H

