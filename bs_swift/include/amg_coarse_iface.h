/** 
 * @file amg_coarse_iface.h
 * @brief interface for coarse grid builder
 * @author 
 * @version 
 * @date 2010-03-01
 */
#ifndef __AMG_CORSE_IFACE_H
#define __AMG_CORSE_IFACE_H
#include <string>

#include "conf.h"
#include "bs_object_base.h"
#include "bs_assert.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "bcsr_amg_matrix_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
{

  class amg_coarse_iface: public objbase
    {
    public:
      //! bcsr matrix
      typedef bcsr_amg_matrix_iface                             bcsr_t;
      typedef smart_ptr<bcsr_t, true>                           sp_bcsr_t;

      //! this_t
      typedef amg_coarse_iface                                  this_t;
      //! short name to smart pointer to this class
      typedef smart_ptr<this_t, true>                           sp_this_t;              

      
    //-----------------------------------------
    //  METHODS
    //-----------------------------------------
    public:

      /** 
       * @brief build CF markers and Strength matrix
       * 
       * @param a_matrix        -- <INPUT> A matrix
       * @param wksp            -- <INPUT/OUTPUT> workspace 
       * @param cf_markers      -- <OUTPUT> CF markers 
       * @param s_markers       -- <OUTPUT> non zero elements of Strength matrix
       * 
       * @return number of coarse points (n_coarse_size) or < 0 if error occur
       */
      virtual t_long build (sp_bcsr_t a_matrix, 
                            spv_double wksp, 
                            spv_long cf_markers,
                            spv_long s_markers) = 0;
      //! destructor
      virtual ~amg_coarse_iface ()
        {}

#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const = 0;
#endif //BSPY_EXPORTING_PLUGIN

    };
};

#endif //__AMG_CORSE_IFACE_H
