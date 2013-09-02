/**
* \file   py_prop.cpp
* \brief  Python wrapper for property
* \author Miryanov Sergey
* \date 2008-04-04
*/
//#include "bs_lsolvers_stdafx.h"
#include "py_list_converter.h"
#include "py_prop.h"
#include "prop.h"

using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //! export matrices to python
  void py_export_prop ()
  {
    using namespace boost::python;

    //base_exporter <matrix_iface <seq_vector<float>, seq_vector<int> >,   default_exporter>::export_class ("matrix_iface_fi"); 
    //matrix_exporter::export_base <matrix_iface, py_matrix_iface_exporter> ("matrix_iface");
    base_exporter <prop_iface, py_prop_exporter>::export_class ("prop_iface");

    class_exporter <prop, prop_iface, py_prop_exporter_2>::export_class ("prop");

  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
