/**
 * \file pass_arg_to_python.h
 * \brief helper for passing arguments (parameters) to python functions from c++
 * \author Sergey Miryanov
 * \date 26.06.2009
 * */
#ifndef BS_BOS_CORE_BASE_PASS_ARG_TO_PYTHON_H_
#define BS_BOS_CORE_BASE_PASS_ARG_TO_PYTHON_H_

#include <boost/type_traits.hpp>
#include <boost/mpl/if.hpp>

namespace blue_sky {
namespace python {
namespace tools {

  template <typename T>
  static boost::python::pointer_wrapper <T *>
  pass_arg_to_python (T *p)
  {
    return boost::python::ptr (p);
  }

  template <typename T>
  static
  typename boost::mpl::if_ <typename boost::is_reference <T>::type, boost::reference_wrapper <T>, T>::type
  pass_arg_to_python (T p)
  {
    return boost::is_reference <T>::value ? boost::ref (p) : p;
  }

} // namespace tools
} // namespace python
} // namespace blue_sky



#endif // #ifndef BS_BOS_CORE_BASE_PASS_ARG_TO_PYTHON_H_

