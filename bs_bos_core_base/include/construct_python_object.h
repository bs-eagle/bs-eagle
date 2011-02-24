/** 
 * \file construct_python_object.h
 * \brief function for create and initialize c++ objects from python
 * \author Sergey Miryanov
 * \date 08.05.2009
 * */
#ifndef BS_CONSTRUCT_PYTHON_OBJECT_H_
#define BS_CONSTRUCT_PYTHON_OBJECT_H_

#include "bs_kernel.h"
#include "throw_exception.h"
#include "py_object_handler.h"

#include <boost/python.hpp>

namespace blue_sky {
namespace python {

  template <typename object_t>
  static boost::python::object 
  init_python_object (boost::python::object py_obj)
  {
    using namespace boost::python;
    object return_value = call_method <object> (py_obj.ptr (), "__cons__");

    object_t *obj = extract <object_t *> (py_obj);
    if (!obj)
      {
        bs_throw_exception ("Can't extract c++ object from python object");
      }

    boost::python::detail::initialize_wrapper (py_obj.ptr (), obj);
    return return_value;
  }

  template <typename object_t>
  static smart_ptr <object_t, true>
  construct_python_object ()
  {
    smart_ptr <object_t, true> obj (BS_KERNEL.create_object (object_t::bs_type ()), bs_dynamic_cast ());
    if (!obj)
      {
        bs_throw_exception ("Can't create object: " + object_t::bs_type ().stype_);
      }

    return obj;
  }

  template <typename object_t>
  static void 
  destroy_python_object (const smart_ptr <object_t, true> & /*obj*/)
  {

  }


} // namespace python
} // namespace blue_sky



#endif // #ifndef BS_CONSTRUCT_PYTHON_OBJECT_H_
