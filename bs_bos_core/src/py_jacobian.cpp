/**
 *       \file  py_jacobian.cpp
 *      \brief  Exports python wrapper for jacobian,
 *              see jacobian.h
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  12.05.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "py_jacobian.h"
#include "jacobian.h"

#ifdef BSPY_EXPORTING_PLUGIN
#include "construct_python_object.h"
namespace blue_sky {
namespace python {

  void
  export_jacobian (const char *name)
  {
    using namespace boost::python;

    typedef BS_SP (mbcsr_matrix_iface) (jacobian::*get_matrix_t) () const;
    get_matrix_t get_matrix = &jacobian::get_matrix;

    class_ <jacobian, boost::noncopyable> (name, no_init)
      .def ("__cons__",         make_constructor (construct_python_object <jacobian>))
      .def ("__init__",         make_function (init_python_object <jacobian>))
      .add_property ("solver",  make_function (&jacobian::get_solver, return_value_policy <copy_const_reference> ()),  make_function (&jacobian::set_solver))
      .add_property ("prec",    make_function (&jacobian::get_prec, return_value_policy <copy_const_reference> ()),    make_function (&jacobian::set_prec))
      .def ("matrix",           get_matrix)
      ;
  }

  void
  py_export_jacobian ()
  {
    export_jacobian ("jacobian");
    boost::python::register_ptr_to_python <smart_ptr <jacobian, true> > ();
  }


} // namespace python
} // namespace blue_sky 
#endif
