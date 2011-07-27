/** 
 * @file py_gis.h
 * @brief python interface to WELL GIS 
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-26
 */
#ifndef PY_GIS_M2WDGY4K

#define PY_GIS_M2WDGY4K

#include <string>
#include "table_iface.h"

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

  PY_EXPORTER (py_gis_exporter, default_exporter)
    .def ("get_table",                       &T::get_table, 
        args (""), "Return table with depth and GIS curves")
    .def ("get_prop",                         &T::get_prop, 
        args (""), "Return GIS properties")
    .def ("read_from_las_file",               &T::read_from_las_file, 
        args ("fname"), "Read gis from LAS file format")
    .def ("__str__",                            &T::py_str)
  PY_EXPORTER_END;                               

  //PY_EXPORTER (py_dummy_exporter, default_exporter)
                                                 
  //PY_EXPORTER_END;                               
                                                 
  //! export matrices to python                  
  void py_export_gis ();                    

    
  } // namespace python
} // namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN



#endif /* end of include guard: PY_GIS_M2WDGY4K */
