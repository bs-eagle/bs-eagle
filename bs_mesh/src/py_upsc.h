/** 
 * @file py_upsc.h
 * @brief python wraper for hdm upscaling
 * @author Alina Yapparova
 * @version 
 * @date 2012-20-02
 */
#ifndef PY_UPSC_H

#define PY_UPSC_H


#include <string>
#include "upsc_iface.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "dummy_base.h"
#include "construct_python_object.h"
#include "make_me_happy.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include "export_python_wrapper.h"

#ifdef BSPY_EXPORTING_PLUGIN
namespace blue_sky
  {
  namespace python
    {

  PY_EXPORTER (py_upsc_exporter, default_exporter)
    .def ("get_prop",                         &T::get_prop, 
        args (""), "Return upsc properties")
    .def ("__str__",                            &T::py_str)
    .def ("upscale_grid_zcolumn", &T::upscale_grid_zcolumn)
    .def ("upscale_grid", &T::upscale_grid)
    .def ("upscale_cubes_xy", &T::upscale_cubes_xy)
    .def ("upscale_sat_cube_xy", &T::upscale_sat_cube_xy)
    .def ("upscale_permz_zcolumn", &T::upscale_permz_zcolumn)
    .def ("upscale_perm_block", &T::upscale_perm_block)
    .def ("upscale_saturation_cube", &T::upscale_saturation_cube)
    .def ("king_method", &T::king_method)
  PY_EXPORTER_END;                               

  //! export matrices to python                  
  void py_export_upsc ();                    

    
  } // namespace python
} // namespace blue_sky

#endif // #ifdef BSPY_EXPORTING_PLUGIN
#endif /* end of include guard: PY_UPSC_H */

