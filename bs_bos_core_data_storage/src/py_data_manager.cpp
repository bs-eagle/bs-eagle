#include "bs_bos_core_data_storage_stdafx.h"

#include "py_data_manager.h"
#include "data_manager.h"
#include "py_data_class.h"

#include "export_python_wrapper.h"

using namespace boost::python;

namespace blue_sky {
namespace python {

  template <typename T>
  smart_ptr <idata, true>
  get_data (T *t)
  {
    return t->data;
  }

  PY_EXPORTER (data_manager_exporter, default_exporter)
    .add_property ("data", get_data <T>)
  PY_EXPORTER_END;

  void py_export_data_manager ()
  {
    strategy_exporter::export_base <data_manager, data_manager_exporter> ("data_manager");
  }

} // namespace python
} // namespace blue_sky
