/** 
 * @file py_table.cpp
 * @brief Python #table and #table_iface interface
 * @author Oleg Borschuk
 * @version 1.0
 * @date 2011-02-26
 */

#include "py_table.h"
#include "table.h"

using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //! export matrices to python
  void py_export_table ()
  {
    using namespace boost::python;

    //base_exporter <matrix_iface <seq_vector<float>, seq_vector<int> >,   default_exporter>::export_class ("matrix_iface_fi"); 
    //matrix_exporter::export_base <matrix_iface, py_matrix_iface_exporter> ("matrix_iface");
    base_exporter <table_iface, py_table_exporter>::export_class ("table_iface");

    class_exporter <table, table_iface, py_table_exporter>::export_class ("table");

  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
