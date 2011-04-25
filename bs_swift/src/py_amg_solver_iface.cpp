/**
* \file   py_amg_solver_iface.cpp
* \brief  Python wrapper for amg solver
* \author Sayfullin Ilshat
* \date 2011-04-01
*/

#include "py_amg_solver_iface.h"
#include "amg_solver.h"
#include "py_list_converter.h"

using namespace boost::python;
#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky {
namespace python {
  //////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////
  //! export to python
  void py_export_amg_solver ()
  {
    using namespace boost::python;

    class_exporter<amg_solver_iface, lsolver_iface, py_amg_solver_iface_exporter>::export_class ("amg_solver_iface");
    class_exporter<amg_solver, amg_solver_iface, py_amg_solver_iface_exporter>::export_class ("amg_solver");

	  // register vector of type descriptors <-> Python list converters
	  typedef bspy_converter< list_traits< std::vector<smart_ptr<bcsr_amg_matrix_iface, true> > > > spbcsr_vec_converter;
	  spbcsr_vec_converter::register_from_py ();
	  spbcsr_vec_converter::register_to_py ();

	  // register vector of type descriptors <-> Python list converters
	  typedef bspy_converter< list_traits< std::vector<spv_long> > > spvlong_vec_converter;
	  spvlong_vec_converter::register_from_py ();
	  spvlong_vec_converter::register_to_py ();

	  // register vector of type descriptors <-> Python list converters
	  typedef bspy_converter< list_traits< std::vector<spv_double> > > spvdouble_vec_converter;
	  spvdouble_vec_converter::register_from_py ();
	  spvdouble_vec_converter::register_to_py ();
  }

}	// namespace python
}	// namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN

