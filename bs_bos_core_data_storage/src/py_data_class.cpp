#include "bs_bos_core_data_storage_stdafx.h"

#include "py_data_class.h"
#include "numpy_tools.h"
#include "py_object_handler.h"

#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>

#include "export_python_wrapper.h"

#ifdef BSPY_EXPORTING_PLUGIN
using namespace boost::python;
namespace bp = boost::python;
using namespace blue_sky::tools;

namespace blue_sky {
namespace python {

  PY_EXPORTER (idata_exporter, default_exporter)
    .def("init", &T::init)
    .def("flush_pool", &T::flush_pool)
    .def("set_i_array", &T::set_i_array)
    .def("set_fp_array", &T::set_fp_array)
    .def("get_i_array", &T::get_i_array)
    .def("get_fp_array", &T::get_fp_array)
  PY_EXPORTER_END;

  void py_export_idata()
  {
    base_exporter<idata, idata_exporter>::export_class ("idata");
  }

} // namespace python
} // namespace blue_sky
#endif