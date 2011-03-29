/** 
 * @file coarse_tools.h
 * @brief 
 * @author 
 * @version 
 * @date 2010-03-03
 */
#ifndef __COARSE_TOOLS_H
#define __COARSE_TOOLS_H

#include "bs_object_base.h"
#include "bs_assert.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "strategies.h"
#include "bcsr_amg_matrix_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {
    class BS_API_PLUGIN coarse_tools: public objbase
    {
      //! bcsr matrix
      typedef bcsr_amg_matrix_iface                             bcsr_t;
      typedef smart_ptr<bcsr_t, true>                           sp_bcsr_t;
      
      typedef coarse_tools                                      this_t;
      BLUE_SKY_TYPE_DECL_T(coarse_tools);
    //-----------------------------------------
    //  METHODS
    //-----------------------------------------
    public:

      /** 
       * @brief build transposed strength matrix
       * 
       * @param matrix          -- <INPUT> matrix
       * @param s_markers       -- <INPUT> strength values
       * @param st_markers      -- <OUTPUT> transposed strength values
       * @param st_cols         -- <OUTPUT> CSR structure for transposed strength matrix
       * @param st_rows         -- <OUTPUT> CSR structure for transposed strength matrix
       * 
       * @return 0 if success 
       */
      static int build_transposed_strength_matrix (sp_bcsr_t matrix, 
                                                   spv_long s_markers,
                                                   spv_long st_markers,
                                                   spv_long st_cols,
                                                   spv_long st_rows
                                                  );
      
      /** 
       * @brief Build independant set
       * 
       * @param matrix          -- <INPUT> matrix
       * @param graph           -- <INPUT> current graph
       * @param last_in_graph   -- <INPUT> 
       * @param meassure_array  -- <INPUT> meassure array for grid points
       * @param ff              -- <INPUT> current ff value
       * @param s_markers       -- <INPUT> strength matrix
       * @param cf_markers      -- <INPUT/OUTPUT> CF markers
       * 
       * @return 0 if success 
       */
      static int build_independent_set (sp_bcsr_t matrix,
                                        spv_long graph,
                                        const t_long last_in_graph,
                                        spv_double meassure_array,
                                        const t_long ff,
                                        spv_long s_markers,
                                        spv_long cf_markers
                                       );

      //! destructor
      virtual ~coarse_tools ()
        {}
    };
  };
#endif //__COARSE_TOOLS_H

