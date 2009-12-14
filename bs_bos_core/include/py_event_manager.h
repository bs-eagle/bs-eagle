/**
 *       \file  py_event_manager.h
 *      \brief  Python wrapper for event_manager, 
 *              see event_manager.h
 *     \author  Morozov Andrey
 *       \date  07.06.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_PY_EVENT_MANAGER_H
#define BS_PY_EVENT_MANAGER_H

#ifdef BSPY_EXPORTING_PLUGIN

#include "py_event_base.h"
#include "event_manager.h"

namespace blue_sky {
namespace python {

  /**
   * \brief  Exports wrappers to python
   * */
  void
  py_export_event_manager ();

} // namespace python
} // namespace blue_sky

#endif //#ifdef BSPY_EXPORTING_PLUGIN
#endif //#ifndef BS_PY_EVENT_MANAGER_H
