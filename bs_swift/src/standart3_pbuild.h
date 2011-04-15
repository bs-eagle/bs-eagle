/**
 * @file standart3_pbuild.h
 * @brief
 * @author
 * @version
 * @date 2010-03-16
 */
#ifndef __STANDART3_PBUILD_H
#define __STANDART3_PBUILD_H

#include "amg_pbuild_iface.h"

namespace blue_sky
{
  /**
   * @brief Prolangation matrix builder (Standart algorithm)
   */
  class BS_API_PLUGIN standart3_pbuild: public amg_pbuild_iface
    {
      typedef bcsr_amg_matrix_iface                      bcsr_t;
      typedef smart_ptr<bcsr_t, true>                    sp_bcsr_t;
      //! this_t
      typedef standart3_pbuild                           this_t;

      BLUE_SKY_TYPE_DECL (standart3_pbuild);
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
                         spv_long sp_cf_markers,
                         spv_long sp_s_markers,
                         sp_bcsr_t p_matrix);
      //! destructor
      virtual ~standart3_pbuild ()
        {}

#ifdef BSPY_EXPORTING_PLUGIN
      virtual std::string py_str () const
        {
          return std::string ("Prolongation matrix builder (Standart3 algorithm).");
        }
#endif //BSPY_EXPORTING_PLUGIN

    //-----------------------------------------
    //  VARIABLES
    //-----------------------------------------
    protected:
      spv_long sp_cur_connections;
      spv_long sp_cur_connections_markers;

    };
}
#endif //__STANDART3_PBUILD_H

