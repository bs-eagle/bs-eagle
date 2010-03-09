/**
 *       \file  py_jacobian.cpp
 *      \brief  Exports python wrapper for jacobian,
 *              see jacobian.h
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  12.05.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#include "stdafx.h"

#include "py_jacobian.h"
#include "jacobian.h"

#ifdef BSPY_EXPORTING_PLUGIN
#include "construct_python_object.h"
namespace blue_sky {
namespace python {

  template <typename jacobian_t>
  void
  export_jacobian (const char *name)
  {
    using namespace boost::python;

    class_ <jacobian_t, boost::noncopyable> (name, no_init)
      .def ("__cons__",         make_constructor (construct_python_object <jacobian_t>))
      .def ("__init__",         make_function (init_python_object <jacobian_t>))
      .add_property ("solver",  make_function (&jacobian_t::get_solver, return_value_policy <copy_const_reference> ()),  make_function (&jacobian_t::set_solver))
      .add_property ("prec",    make_function (&jacobian_t::get_prec, return_value_policy <copy_const_reference> ()),    make_function (&jacobian_t::set_prec))
      .add_property ("jmx",     make_function (&jacobian_t::get_jmatrix, return_value_policy <copy_const_reference> ()), make_function (&jacobian_t::set_jmatrix))
      ;
  }

  void
  py_export_jacobian ()
  {
    export_jacobian <jacobian <base_strategy_fi> >   ("jacobian_fi");
    export_jacobian <jacobian <base_strategy_di> >   ("jacobian_di");
    export_jacobian <jacobian <base_strategy_mixi> > ("jacobian_mixi");

    boost::python::register_ptr_to_python <smart_ptr <jacobian <base_strategy_fi>, true> > ();
    boost::python::register_ptr_to_python <smart_ptr <jacobian <base_strategy_di>, true> > ();
    boost::python::register_ptr_to_python <smart_ptr <jacobian <base_strategy_mixi>, true> > ();
  }


} // namespace python
} // namespace blue_sky 
#endif