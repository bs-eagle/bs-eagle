/** 
 * @file amg_smbuilder_iface.h
 * @brief 
 * @author 
 * @version 
 * @date 2010-03-09
 */
#ifndef __AMG_SMBUILDER_IFACE_H
#define __AMG_SMBUILDER_IFACE_H

#include <string>

#include "bs_object_base.h"
#include "bs_assert.h"

#include "conf.h"
#include BS_FORCE_PLUGIN_IMPORT ()
#include "bcsr_amg_matrix_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
{

  class amg_smbuilder_iface: public objbase
    {
    public:
      typedef bcsr_amg_matrix_iface                             bcsr_t;
      typedef smart_ptr<bcsr_t, true>                           sp_bcsr_t;

      //! this_t
      typedef amg_smbuilder_iface                               this_t;
      //! short name to smart pointer to this class
      typedef smart_ptr<this_t, true>                           sp_this_t;              

      
    //-----------------------------------------
    //  METHODS
    //-----------------------------------------
    public:

      /** 
       * @brief build strength matrix Strength matrix
       * 
       * @param a_matrix                -- <INPUT> matrix
       * @param strength_threshold      -- <INPUT> strength threshold
       * @param max_row_sum             -- <INPUT>
       * @param s_markers               -- <OUTPUT> strength matrix 
       * 
       * @return 0 if success
       */
      virtual int build (sp_bcsr_t a_matrix, 
                         t_double strength_threshold, 
                         t_double max_row_sum,
                         spv_long s_markers) const = 0;
      //! destructor
      virtual ~amg_smbuilder_iface ()
        {}

#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const = 0;
#endif //BSPY_EXPORTING_PLUGIN

    };
};
#endif //__AMG_SMBUILDER_IFACE_H
