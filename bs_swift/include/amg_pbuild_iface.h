/** 
 * @file amg_pbuild_iface.h
 * @brief interface for AMG Prolangation matrix builder
 * @author 
 * @version 
 * @date 2010-03-01
 */
#ifndef __AMG_PBUILD_H
#define __AMG_PBUILD_H

#include "bs_object_base.h"
#include "bs_assert.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "strategies.h"
#include "bcsr_amg_matrix_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {
  /**
   * \brief Interpolation base
   */
  class amg_pbuild_iface: public objbase
    {
    public:
      //! bcsr matrix
      typedef bcsr_amg_matrix_iface                             bcsr_t;
      typedef smart_ptr<bcsr_t, true>                           sp_bcsr_t;

      //! this_t
      typedef amg_pbuild_iface                                  this_t;
      //! short name to smart pointer to this class
      typedef smart_ptr<this_t, true>                           sp_this_t;              


    //-----------------------------------------
    //  METHODS
    //-----------------------------------------
    public:
      //! destructor
      virtual ~amg_pbuild_iface ()
        {}

      /** 
       * @brief build prolangation matrix
       * 
       * @param a_matrix        -- <INPUT> A matrix
       * @param n_coarse_size   -- <INPUT> number of points in coarse grid
       * @param max_connections -- <INPUT> maximum allowed connections 
       * @param cf_markers      -- <INPUT/OUTPUT> CF markers
       * @param s_markers       -- <INPUT> non zero elements of Strength matrix
       * @param p_matrix        -- <OUTPUT> prolangation matrix
       * 
       * @return 
       */
      virtual int build (sp_bcsr_t a_matrix, 
                         const t_long n_coarse_size,
                         const t_long max_connections,
                         spv_long cf_markers,
                         spv_long s_markers,
                         sp_bcsr_t p_matrix) = 0;

#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const = 0;
#endif //BSPY_EXPORTING_PLUGIN
    };
  };
#endif //__AMG_PBUILD_H
