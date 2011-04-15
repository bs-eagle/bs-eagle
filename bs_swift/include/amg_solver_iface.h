/**
 * @file amg_solver_iface.h
 * @brief Algebraic Multi Grid linear system solver interface file
 * @author
 * @date 2009-12-07
 */
#ifndef __AMG_SOLVER_IFACE_H
#define __AMG_SOLVER_IFACE_H


#include "bs_assert.h"
#include "bs_tree.h"
#include "conf.h"
#include BS_FORCE_PLUGIN_IMPORT ()
#include "lsolver_iface.h"
#include "amg_smbuilder_iface.h"
#include "amg_coarse_iface.h"
#include "amg_pbuild_iface.h"
#include "amg_smoother_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
{
  /**
   * @brief interface class for matrix storage and manipulation
   */
  class amg_solver_iface: public lsolver_iface
    {
    public:
      typedef amg_solver_iface                                  this_t;
      //! short name to smart pointer to this class
      typedef smart_ptr<this_t, true>                           sp_this_t;

      typedef smart_ptr<amg_smbuilder_iface, true>              sp_smbuild_t;
      typedef smart_ptr<amg_coarse_iface, true>                 sp_coarse_t;
      typedef smart_ptr<amg_pbuild_iface, true>                 sp_pbuild_t;
      typedef smart_ptr<amg_smoother_iface, true>               sp_smooth_t;

      typedef smart_ptr<bcsr_amg_matrix_iface, true>            sp_bcsr_t;
      typedef std::vector<sp_bcsr_t>                            vec_sp_bcsr_t;
      //-----------------------------------------
      //  METHODS
      //-----------------------------------------

      //virtual const bcsr_t &get_A (int layer) const = 0;

      //virtual const bcsr_t &get_P (int layer) const = 0;

      //virtual const fp_vector_type_t &get_rhs (int layer)

      virtual void set_smbuilder (unsigned int level, sp_smbuild_t sp_smbuild_iface) = 0;
      virtual void set_coarser (unsigned int level, sp_coarse_t sp_coarse_iface) = 0;
      virtual void set_pbuilder (unsigned int level, sp_pbuild_t sp_pbuild_iface) = 0;
      virtual void set_pre_smoother (unsigned int level, sp_smooth_t sp_smooth_iface) = 0;
      virtual void set_post_smoother (unsigned int level, sp_smooth_t sp_smooth_iface) = 0;
      virtual const vec_sp_bcsr_t get_level_matrices ()  const = 0;
    public:

    public:
      //! destructor
      virtual ~amg_solver_iface ()
        {}

    };
};
#endif //__AMG_SOLVER_IFACE_H
