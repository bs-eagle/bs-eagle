/** 
 * \file numpy_tools.h
 * \brief functions to extract data from numpy objects
 * \author Sergey Miryanov
 * \date 25.05.2009
 * */
#ifndef BS_BOS_CORE_DATA_STORAGE_NUMPY_TOOLS_H_
#define BS_BOS_CORE_DATA_STORAGE_NUMPY_TOOLS_H_

#include <boost/format.hpp>

namespace blue_sky {
namespace numpy_tools {
namespace type_helper
{
  template <typename T>
  struct type_kind
  {
  };

  template <>
  struct type_kind <float>
  {
    enum {
      numpy_type = 'f',
    };
    //static const char numpy_type  = 'f';
    //static const char cpp_type    = 'f';
  };
  template <>
  struct type_kind <int>
  {
    enum {
      numpy_type = 'i',
    };
    //static const char numpy_type  = 'i';
    //static const char cpp_type    = 'i';
  };
  template <>
  struct type_kind <unsigned char>
  {
    enum {
      numpy_type = 'u',
    };
    //static const char numpy_type  = 'u';
    //static const char cpp_type    = 'u';
  };
}

  template <typename item_t>
  item_t *
  get_buffer (const boost::python::object &obj, boost::python::ssize_t size)
  {
    using namespace boost;
    using namespace boost::python;
    namespace bp = boost::python;

    std::string kind = extract <std::string> (obj.attr ("dtype").attr ("kind"));
    BS_ASSERT (kind.length ()) ((char)type_helper::type_kind <item_t>::numpy_type);

    if ((char)type_helper::type_kind <item_t>::numpy_type != kind[0])
      {
        bs_throw_exception((format ("Types mismatch: %s -> %s") % (char)type_helper::type_kind <item_t>::numpy_type % kind).str ());
      }

    bp::ssize_t item_size = extract <bp::ssize_t> (obj.attr ("itemsize"));
    if (sizeof (item_t) != item_size)
      {
        bs_throw_exception ((format ("Invalid itemsize: %d -> %d") % sizeof (item_t) % item_size).str ());
      }

    object buf = obj.attr ("data");
    PyObject *buffer = buf.ptr ();
    BS_ASSERT (buffer);

    PyBufferProcs *pb = buffer->ob_type->tp_as_buffer;
    BS_ASSERT (pb);

    void *result = 0;
    bp::ssize_t buf_size = (*pb->bf_getwritebuffer) (buffer, 0, &result);
    if (buf_size != static_cast <bp::ssize_t> (size * sizeof (item_t)))
      {
        bs_throw_exception ((format ("Invalid buffer size: %d -> %d") % (size * sizeof (item_t)) % buf_size).str ());
      }

    BS_ASSERT (result);
    return (item_t *)result;
  }

  boost::python::ssize_t 
  get_len (const boost::python::object &obj)
  {
    namespace bp = boost::python;
    return bp::len (obj.attr ("data")) / (bp::ssize_t)bp::extract<bp::ssize_t> (obj.attr ("itemsize"));
  }


} // namespace numpy_tools
} // namespace blue_sky

#endif // #ifndef BS_BOS_CORE_DATA_STORAGE_NUMPY_TOOLS_H_
