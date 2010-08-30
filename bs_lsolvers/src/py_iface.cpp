/**
* \file   py_matrix_iface.cpp
* \brief  Python wrapper for linear solvers
* \author Miryanov Sergey
* \date 2008-04-04
*/
#include "bs_lsolvers_stdafx.h"
#include "py_iface.h"
#include "matrix_iface.h"
#include "bicgstab_solver.h"
#include "cgs_solver.h"
#include "gmres_solver.h"
#include "bcsr_ilu_prec.h"
#include "two_stage_prec.h"
#include "tfqmr_solver.h"
#include "blu_solver.h"
#include "gs_solver.h"


using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //! export matrices to python
  void py_export_lsolvers ()
  {
    using namespace boost::python;

    strategy_exporter::export_base <lsolver_iface, py_lsolver_iface_exporter> ("lsolver_iface");

    strategy_exporter::export_class <bicgstab_solver,   lsolver_iface, py_lsolver_iface_exporter> ("bicgstab");
    strategy_exporter::export_class <cgs_solver,        lsolver_iface, py_lsolver_iface_exporter> ("cgs");
    strategy_exporter::export_class <gmres_solver,      lsolver_iface, py_lsolver_iface_exporter> ("gmres");
    strategy_exporter::export_class <bcsr_ilu_prec,     lsolver_iface, py_lsolver_iface_exporter> ("bcsr_ilu_prec");
    strategy_exporter::export_class <blu_solver,        lsolver_iface, py_lsolver_iface_exporter> ("blu_solver");
    strategy_exporter::export_class <two_stage_prec,    lsolver_iface, py_lsolver_iface_exporter> ("two_stage_prec");
    strategy_exporter::export_class <tfqmr_solver,      lsolver_iface, py_lsolver_iface_exporter> ("tfqmr");
    strategy_exporter::export_class <amg_smoother_iface,lsolver_iface, py_amg_smoother_iface_exporter> ("amg_smoother_iface");
    strategy_exporter::export_class <gs_solver,         amg_smoother_iface, py_amg_smoother_iface_exporter> ("gs");

  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN

