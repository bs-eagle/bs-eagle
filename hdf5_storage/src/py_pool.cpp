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
#include "bs_serialize.h"

using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //! export matrices to python
  void py_export_pool ()
  {
    using namespace boost::python;

    base_exporter <h5_pool_iface, py_pool_exporter>::export_class ("h5_pool_iface");

    class_exporter <h5_pool, h5_pool_iface, py_pool_exporter>::export_class ("h5_pool");

    def("serialize_to_str", &blue_sky::serialize_to_str< h5_pool_iface >);
    def("serialize_from_str", &blue_sky::serialize_from_str< h5_pool_iface >);

    // register implicit conversion to interface
    implicitly_convertible<
      smart_ptr< h5_pool >,
      smart_ptr< h5_pool_iface >
    >();
  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
