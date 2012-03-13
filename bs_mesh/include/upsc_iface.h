/** 
 * @file upsc_iface.h
 * @brief interface for the upscaling
 * @author Alina Yapparova
 * @version 
 * @date 2012-20-02
 */
#ifndef UPSC_IFACE_H

#define UPSC_IFACE_H


#include <string>

#include "bs_object_base.h"
#include "conf.h"

#include "bs_mesh_stdafx.h"
#include "bs_mesh_grdecl.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include BS_STOP_PLUGIN_IMPORT ()

#include "prop_iface.h"

namespace blue_sky
{
class upsc_iface : public objbase
  {
    public:
      typedef BS_SP (prop_iface)                      sp_prop_t;

    public:
      /** 
       * @brief destructor
       */
      virtual ~upsc_iface ()
        {}

      /** 
       * @brief return SP to the property 
       */
      virtual sp_prop_t get_prop () = 0;
      virtual boost::python::tuple king_method ( t_long nx, t_long ny, t_long nz, t_long nz_upsc,
                            spv_float vol, spv_float ntg, spv_float poro, spv_float perm, 
                            spv_float swat ) = 0;
      
      virtual spv_float upscale_grid ( t_long Nx, t_long Ny, t_long Nz, spv_float zcorn, spv_uint layers ) = 0;

#ifdef BSPY_EXPORTING_PLUGIN
      /** 
       * @brief python print wrapper
       * 
       * @return return table description
       */
      virtual std::string py_str () const = 0;

#endif //BSPY_EXPORTING_PLUGIN
};

}  // end of bluesky name space

#endif /* end of include guard: UPSC_IFACE_H */
