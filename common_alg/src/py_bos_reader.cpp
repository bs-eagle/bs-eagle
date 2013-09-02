/** 
 * @file py_bos_reader.cpp
 * @brief python wraper for BOS ascii file reader
 * @author Oleg Borschuk
 * @version 
 * @date 2012-03-01
 */

#include "py_bos_reader.h"
#include "bos_reader.h"

using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {

  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //! export matrices to python
  void py_export_bos_reader ()
  {
    using namespace boost::python;

    base_exporter <bos_reader_iface, py_bos_reader_exporter>::export_class ("bos_reader_iface");
    //base_exporter <bos_reader, py_bos_reader_exporter>::export_class ("bos_reader");

    class_exporter <bos_reader, bos_reader_iface, py_bos_reader_exporter>::export_class ("bos_reader");


    // register implicit conversion to interface
    implicitly_convertible<
      smart_ptr< bos_reader >,
      smart_ptr< bos_reader_iface >
    >();
  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
