/** 
 * @file py_pool.h
 * @brief python wrapper for #h5_pool_iface and #h5_pool
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-03-12
 */
#ifndef PY_POOL_UHSZMDSY

#define PY_POOL_UHSZMDSY


#include <string>
#include "pool_iface.h"

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

  PY_EXPORTER (py_pool_exporter, default_exporter)
    .def ("open_file",                          &T::open_file, 
        args ("fname", "path"), "Open or create h5 file")
    .def ("close_file",                         &T::close_file, 
        args (""), "Close h5 file")
    .def ("flush",                              &T::flush, 
        args (""), "Flush buffers to the file")
    .def ("get_fp_data",                        &T::get_fp_data, 
        args ("name"), "Return array with given name")
    .def ("get_i_data",                         &T::get_i_data, 
        args ("name"), "Return array with given name")
    .def ("set_fp_data",                        &T::set_fp_data, 
        args ("name", "array"), "Store array in the pool")
    .def ("set_i_data",                         &T::set_i_data, 
        args ("name", "array"), "Store array in the pool")
    .def ("get_data_type",                      &T::get_data_type, 
        args ("name"), "Get array data type")
    .def ("list_data",                          &T::py_list_data)
    .def ("py_set_pool_dims",                   &T::py_set_pool_dims, 
        args ("dims"), "Set pool dimensions")
    .def ("py_get_pool_dims",                   &T::py_get_pool_dims, 
        args (""), "Get pool dimensions")
    .def ("__str__",                            &T::py_str)
  PY_EXPORTER_END;                               

  //! export matrices to python                  
  void py_export_pool ();                    

    
  } // namespace python
} // namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN

#endif /* end of include guard: PY_POOL_UHSZMDSY */
