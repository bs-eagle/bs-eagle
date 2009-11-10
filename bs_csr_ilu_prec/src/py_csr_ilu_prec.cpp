/**
 * \file py_csr_ilu_prec.cpp
 * \brief Python wrappers for csr_ilu_prec
 * \author Miryanov Sergey
 * \date 18.04.2008
 */
#include "bs_csr_ilu_prec_stdafx.h"
#include "py_csr_ilu_prec.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "py_linear_solvers.h"
#include BS_STOP_PLUGIN_IMPORT ()

#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {


  template <class solver_t>
  struct csr_ilu_exporter
  {
    template <typename class_t>
    static class_t &
    export_class (class_t &class__)
    {
      using namespace boost::python;

      solver_exporter <solver_t>::export_class (class__)
        .add_property ("ilu_matrix", make_function (&solver_t::get_ilu_matrix))
        ;

      return class__;
    }
  };
  

  //! export wrappers to python
  void
  py_export_csr_ilu_prec ()
  {
    strategy_exporter::export_class <csr_ilu_prec, linear_solver_base, csr_ilu_exporter> ("csr_ilu_prec_seq");
  }

} // namespace python
} // namespace blue_sky

#endif // #ifdef BSPY_EXPORTING_PLUGIN
