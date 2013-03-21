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

BLUE_SKY_TYPE_SERIALIZE_GUID(blue_sky::h5_pool)
BLUE_SKY_CLASS_SRZ_FCN_DECL(serialize, blue_sky::h5_pool)

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

    std::string (*s2s_hpi)(smart_ptr< h5_pool_iface, true >&) = &blue_sky::serialize_to_str< h5_pool_iface >;
    std::string (*s2s_pi)(smart_ptr< h5_pool, true >&) = &blue_sky::serialize_to_str< h5_pool >;
    smart_ptr< h5_pool_iface, true > (*sfs_hpi)(const std::string&) = &blue_sky::serialize_from_str< h5_pool_iface >;
    smart_ptr< h5_pool, true > (*sfs_pi)(const std::string&) = &blue_sky::serialize_from_str< h5_pool >;

    def("serialize_to_str", s2s_hpi);
    def("serialize_to_str", s2s_pi);
    def("serialize_from_str", sfs_hpi);
    def("serialize_from_str", sfs_pi);

    // register implicit conversion to interface
    implicitly_convertible<
      smart_ptr< h5_pool >,
      smart_ptr< h5_pool_iface >
    >();
  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
