#include "stdafx.h"
#include "py_equil_model.h"
#include "equil_model_depth.h"

#ifdef BSPY_EXPORTING_PLUGIN
#include "export_python_wrapper.h"
#include <boost/python/call_method.hpp>
#include "py_list_converter.h"

using namespace boost::python;
namespace blue_sky
{
	namespace python
	{
		PY_EXPORTER(equil_model_depth_exporter, default_exporter)
			.def("calc_equil", &T::py_calc_equil)
			.def("get_pressure", &T::get_pressure)
			.def("get_saturation", &T::get_saturation)
		PY_EXPORTER_END;

		void
		py_export_equil_model ()
		{
			base_exporter <equil_model_depth, equil_model_depth_exporter>::export_class ("equil_model_depth");

			typedef bspy_converter< list_traits< std::vector<double> > > spv_table_list_converter;
	        spv_table_list_converter::register_from_py ();
	        spv_table_list_converter::register_to_py ();
		}
	}
}

#endif