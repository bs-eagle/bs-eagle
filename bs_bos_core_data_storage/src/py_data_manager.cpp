#include "bs_bos_core_data_storage_stdafx.h"

#include "py_data_manager.h"
#include "data_manager.h"
#include "py_data_class.h"

#include "export_python_wrapper.h"

#ifdef BSPY_EXPORTING_PLUGIN
using namespace boost::python;

namespace blue_sky {
namespace python {

 /*
  template <typename T, class strategy_t>
  smart_ptr <idata, true>
  get_data (T *t)
  {
    return t->data;
  }
*/
  PY_EXPORTER (data_manager_exporter, default_exporter)
    //.add_property ("data", get_data <T, strategy_t>)
    .def("init", &T::init)
    .def("read_keyword_file", &T::read_keyword_file)
    .def("get_data", &T::get_data)
    .def("get_reader", &T::get_reader)
    .def("get_keyword_manager", &T::get_keyword_manager)
  PY_EXPORTER_END;
  void py_export_data_manager ()
  {
    //strategy_exporter::export_base <data_manager, data_manager_exporter> ("data_manager");
  }

} // namespace python
} // namespace blue_sky
#endif