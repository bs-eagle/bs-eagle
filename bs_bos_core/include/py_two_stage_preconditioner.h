/**
 *       \file  py_two_stage_preconditioner.h
 *      \brief  Python wrapper for two_stage_preconditioner, 
 *              see two_stage_preconditioner.h
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  16.04.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_PY_TWO_STAGE_PRECONDITIONER_H_
#define BS_PY_TWO_STAGE_PRECONDITIONER_H_

#ifdef BSPY_EXPORTING_PLUGIN
#include "two_stage_preconditioner.h"

namespace blue_sky {
namespace python {

  //! export classes to python
  void py_export_two_stage_prec ();

} // namespace python
} // namespace blue_sky

#endif // #ifdef BSPY_EXPORTING_PLUGIN
#endif // #ifndef BS_PY_TWO_STAGE_PRECONDITIONER_H_
