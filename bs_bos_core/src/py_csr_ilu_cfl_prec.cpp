/**
 * \file py_csr_ilu_cfl_prec.cpp
 * \brief Python wrappers for csr_ilu_cfl_prec
 * \author Salimgareeva Elmira
 * \date 28.09.2009
 */
#include "stdafx.h"
#include "py_csr_ilu_cfl_prec.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "py_linear_solvers.h"
#include BS_STOP_PLUGIN_IMPORT ()

#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {


  //! export wrappers to python
  void
  py_export_csr_ilu_cfl_prec ()
  {
    strategy_exporter::export_class <csr_ilu_cfl_prec, linear_solver_base, default_exporter> ("csr_ilu_cfl_prec_seq");
  }

} // namespace python
} // namespace blue_sky

#endif // #ifdef BSPY_EXPORTING_PLUGIN

