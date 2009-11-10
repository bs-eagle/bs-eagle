/**
 * \file py_two_stage_preconditioner.h
 * \brief Python wrapper for two_stage_preconditioner
 * \author Miryanov Sergey
 * \date 16.04.2008
 */

#ifndef BS_PY_TWO_STAGE_PRECONDITIONER_H_
#define BS_PY_TWO_STAGE_PRECONDITIONER_H_

#ifdef BSPY_EXPORTING_PLUGIN
#include "two_stage_preconditioner.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "py_bcsr_matrix.h"
#include "py_linear_solvers.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky {
namespace python {

  //! export classes to python
  void py_export_two_stage_prec ();

} // namespace python
} // namespace blue_sky

#endif // #ifdef BSPY_EXPORTING_PLUGIN
#endif // #ifndef BS_PY_TWO_STAGE_PRECONDITIONER_H_
