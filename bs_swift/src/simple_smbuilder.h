/**
 * @file simple_smbuilder.h
 * @brief
 * @author
 * @version
 * @date 2010-03-09
 */
#ifndef __SIMPLESMBUILDER_H
#define __SIMPLESMBUILDER_H

#include "amg_smbuilder_iface.h"

namespace blue_sky
{
  class BS_API_PLUGIN simple_smbuilder: public amg_smbuilder_iface
    {
      //! bcsr matrix
      typedef bcsr_amg_matrix_iface                             bcsr_t;
      typedef smart_ptr<bcsr_t, true>                           sp_bcsr_t;
      //! this_t
      typedef simple_smbuilder                                  this_t;

      BLUE_SKY_TYPE_DECL (simple_smbuilder);
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
       * @return max_connections
       */
      virtual t_long build (sp_bcsr_t a_matrix,
                            t_double strength_threshold,
                            t_double max_row_sum,
                            spv_long s_markers) const;

      //! destructor
      virtual ~simple_smbuilder ()
        {}

#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const
        {
          return std::string ("Simple strength matrix builder.");
        }
#endif //BSPY_EXPORTING_PLUGIN

    //-----------------------------------------
    //  VARIABLES
    //-----------------------------------------
    protected:
    };
}
#endif //__SIMPLESMBUILDER_H

