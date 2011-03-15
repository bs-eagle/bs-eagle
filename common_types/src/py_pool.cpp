/** 
 * @file py_pool.cpp
 * @brief python wrapper for #h5_pool_iface and #h5_pool
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-03-12
 */

#include "py_pool.h"
#include "h5_pool.h"

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


  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
