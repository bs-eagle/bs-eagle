/** 
 * @file conf.cpp
 * @brief conf helper methods
 * @author Mark Khait
 * @version 1.0
 * @date 2011-04-22
*/
#include "bs_bos_core_base_stdafx.h"
#include "conf.h"


namespace blue_sky {

std::string get_long_type ()
{
  if (sizeof (t_long) == sizeof (int))
    return "int32";
  else if (sizeof (t_long) == sizeof (long))
    return "int64";
  else return "unknown type";
}

std::string get_float_type ()
{
  if (sizeof (t_float) == sizeof (float))
    return "float32";
  else if (sizeof (t_float) == sizeof (double))
    return "float64";
  else return "unknown type";
}

std::string get_double_type ()
{
  if (sizeof (t_double) == sizeof (float))
    return "float32";
  else if (sizeof (t_double) == sizeof (double))
    return "float64";
  else return "unknown type";
}


namespace python {

  void
  export_type_helper ()
  {
    using namespace boost::python;
    def("get_long_type", get_long_type);
    def("get_double_type", get_double_type);
    def("get_float_type", get_float_type);
  }

} // namespace python
} // namespace blue_sky


