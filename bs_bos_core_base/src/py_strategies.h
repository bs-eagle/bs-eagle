#ifndef __PY_STRATEGIES_H
#define __PY_STRATEGIES_H

#include "strategies.h"
#include BS_FORCE_PLUGIN_IMPORT ()
#include "dummy_base.h"
#include "construct_python_object.h"
#include "make_me_happy.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include "export_python_wrapper.h"

#ifdef BSPY_EXPORTING_PLUGIN
namespace blue_sky
  {
  namespace python
    {

  //! export matrices to python                  
  void py_export_strategies ();                    
    
  } // namespace python
} // namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN

#endif //__PY_STRATEGIES_H
