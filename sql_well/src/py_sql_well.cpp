/** 
 * @file py_sql_well.cpp
 * @brief python wraper for well fracture storage
 * @author Oleg Borschuk
 * @version 
 * @date 2011-07-29
 */

#include "py_sql_well.h"
#include "sql_well.h"
#include "bs_serialize.h"

using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {
  void py_export_compdat_ident();

  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //! export matrices to python
  void py_export_sql_well ()
  {
    using namespace boost::python;

    base_exporter <well_pool_iface, py_sql_well_exporter>::export_class ("well_pool_iface");

    class_exporter <sql_well, well_pool_iface, py_sql_well_exporter>::export_class ("sql_well");

    py_export_compdat_ident();

    def("serialize_to_str", &blue_sky::serialize_to_str< well_pool_iface >);
    def("serialize_from_str", &blue_sky::serialize_from_str< well_pool_iface >);

    // register implicit conversion to interface
    implicitly_convertible<
      smart_ptr< sql_well >,
      smart_ptr< well_pool_iface >
    >();
  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
