/** 
 * @file py_pool.cpp
 * @brief python wrapper for #h5_pool_iface and #h5_pool
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-03-12
 */
#ifdef BSPY_EXPORTING_PLUGIN
#include <boost/python.hpp>
#endif

#include "py_pool.h"
#include "h5_pool.hpp"
#include "bs_misc.h"

using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {
  //////////////////////////////////////////////////////////////////////////
  void py_export_h5_pool_serialize();

namespace {

  hid_t get_file_id(const h5_pool& p) {
    return p.file_id;
  }
  std::string get_fname(const h5_pool& p) {
    return p.fname;
  }
  // wide string version of open_file
  void open_file_w(h5_pool& p, const std::wstring& fname) {
    p.open_file(
#ifdef UNIX
        wstr2str(fname)
#else
        wstr2str(fname, "ru_RU.CP1251")
#endif
    );
  }
}

  template< typename T >
  struct h5_pool_exporter_plus {
    template< typename class_t >
    static class_t &
    export_class(class_t& class__) {
      py_pool_exporter< T >::export_class(class__)
        .add_property("file_id", &get_file_id)
        .add_property("fname", &get_fname)
        .def("open_file", &open_file_w);
      return class__;
    }
  };

  //////////////////////////////////////////////////////////////////////////
  //! export matrices to python
  void py_export_pool ()
  {
    using namespace boost::python;

    base_exporter <h5_pool_iface, py_pool_exporter>::export_class ("h5_pool_iface");

    class_exporter <h5_pool, h5_pool_iface, h5_pool_exporter_plus>::export_class ("h5_pool");

    py_export_h5_pool_serialize();

    // register implicit conversion to interface
    implicitly_convertible<
      smart_ptr< h5_pool >,
      smart_ptr< h5_pool_iface >
    >();
  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
