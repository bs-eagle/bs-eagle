/**
 *       \file  py_facility_manager.h
 *      \brief  Python wrappers for facility_manager,
 *              see facility_manager.h
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  17.10.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#ifndef PY_FACILITY_MANAGER_H
#define PY_FACILITY_MANAGER_H

#include "py_calc_well.h"
#include "facility_manager.h"

namespace blue_sky {
namespace python {

  /**
   * \brief  Exports wrappers to python
   * \param  
   * \return 
   * */
  void 
  py_export_facility_manager();

}
}

#endif // PY_FACILITY_MANAGER_H

