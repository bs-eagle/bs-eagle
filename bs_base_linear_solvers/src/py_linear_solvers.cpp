/**
* \file   py_linear_solvers.cpp
* \brief  Python wrapper for linear solvers
* \author Miryanov Sergey
* \date 2008-04-04
*/

#include "bs_base_linear_solvers_stdafx.h"
#include "py_linear_solvers.h"

using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {
  //////////////////////////////////////////////////////////////////////////
  //! export linear_solver_prop to python
  void py_export_linear_solver_prop ()
  {
    class_<linear_solver_prop, bases <named_pbase> , boost::noncopyable>("linear_solver_prop", no_init)
      .def ("__cons__",             make_constructor (construct_python_object <linear_solver_prop>))
      .def ("__init__",             make_function (init_python_object <linear_solver_prop>))
      .def ("set_max_iters",        &linear_solver_prop::set_max_iters)
      .def ("set_tolerance",        &linear_solver_prop::set_tolerance)
      .def ("check_convergence",    &linear_solver_prop::check_convergence)
      .def ("get_iters",            &linear_solver_prop::get_iters)
      .def ("set_iters",            &linear_solver_prop::set_iters)
      .def ("get_relative_factor",  &linear_solver_prop::get_relative_factor)
      .def ("get_max_iters",        &linear_solver_prop::get_max_iters)
      .def ("get_tolerance",        &linear_solver_prop::get_tolerance)
      .def ("get_final_resid",      &linear_solver_prop::get_final_resid)
      ;
  }

  //////////////////////////////////////////////////////////////////////////
  //! export linear solvers to python
  void py_export_linear_solvers ()
  {
    using namespace boost::python;

    strategy_exporter::export_base <linear_solver_base, default_exporter> ("linear_solver");

    strategy_exporter::export_class <py_linear_solver,  linear_solver_base, solver_exporter> ("linear_solver_wrapper");
    strategy_exporter::export_class <gmres_solver2,     linear_solver_base, solver_exporter> ("gmres_solver2_seq");
    strategy_exporter::export_class <bicgstab_solver,   linear_solver_base, solver_exporter> ("bicgstab_solver_seq");
    strategy_exporter::export_class <cgs_solver,        linear_solver_base, solver_exporter> ("cgs_solver_seq");
    strategy_exporter::export_class <tfqmr_solver,      linear_solver_base, solver_exporter> ("tfqmr_solver_seq");

#ifdef _MPI
    solver_export <py_linear_solver <gmres_solver2 <mpi_strategy_di>, py_matrix_base <mpi_vector <double>, mpi_vector <int> > > > ("gmres_solver2_mpi_di");
#endif
  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN

