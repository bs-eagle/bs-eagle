#ifndef PY_DATA_CLASS_H
#define PY_DATA_CLASS_H

#include "data_class.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "dummy_base.h"
#include "construct_python_object.h"
#include "make_me_happy.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include "export_python_wrapper.h"

namespace blue_sky {
namespace python {

  /**
   * \brief Exports idata to python
   * */
  void 
  py_export_idata();

} // namespace python
} // namespace blue_sky

#endif // PY_DATA_CLASS_H

