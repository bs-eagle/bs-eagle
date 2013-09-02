/**
 * @file direct_pbuild.h
 * @brief
 * @author
 * @version
 * @date 2010-03-16
 */
#ifndef __DIRECT_PBUILD_H
#define __DIRECT_PBUILD_H

#include "amg_pbuild_iface.h"

namespace blue_sky
{
  /**
   * @brief Prolangation matrix builder (Direct algorithm)
   */
  class BS_API_PLUGIN direct_pbuild: public amg_pbuild_iface
    {
      //! bcsr matrix
      typedef bcsr_amg_matrix_iface                             bcsr_t;
      typedef smart_ptr<bcsr_t, true>                           sp_bcsr_t;
      //! this_t
      typedef direct_pbuild                                     this_t;

      BLUE_SKY_TYPE_DECL (direct_pbuild);
    //-----------------------------------------
    //  METHODS
    //-----------------------------------------
    public:
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
                         sp_bcsr_t p_matrix);
      //! destructor
      virtual ~direct_pbuild ()
        {}

#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const
        {
          return std::string ("Prolongation matrix builder (Direct algorithm).");
        }
#endif //BSPY_EXPORTING_PLUGIN

    //-----------------------------------------
    //  VARIABLES
    //-----------------------------------------
    protected:

    };
}
#endif //__DIRECT_PBUILD_H

