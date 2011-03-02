/** 
 * @file py_mbcsr_matrix.h
 * @brief python wrapper for #mbcsr_matrix
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-03-02
 */
#ifndef PY_MBCSR_MATRIX_774X6T7T

#define PY_MBCSR_MATRIX_774X6T7T

#include "mbcsr_matrix_iface.h"
#include "mbcsr_matrix.h"
#include "py_matrix_iface.h"

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

  PY_EXPORTER (py_mbcsr_matrix_iface_exporter, py_matrix_iface_exporter)
    .def ("get_matrix",                     
          &T::get_matrix,
          args ("name"), "Return matrix by name")
    .def ("add_matrix",                               
          &T::add_matrix,
          args ("name", "matrix"), "Add new matrix to the list (map)")
    .def ("clear",                        
          &T::clear,
          args (""), "Clear all")
    .def ("merge",                     
          &T::merge,
          args (""), "Merge all matrix")
  PY_EXPORTER_END;

  //! export matrices to python                  
  void py_export_mbcsr_matrices ();                    

    
  } // namespace python
} // namespace blue_sky
#endif //


#endif /* end of include guard: PY_MBCSR_MATRIX_774X6T7T */
