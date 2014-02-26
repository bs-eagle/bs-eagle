#include "stdafx.h"
#include "py_equil_model.h"
#include "equil_model_iface.h"
#include "equil_model_depth.h"
#include "bs_serialize.h"

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
		  .def("init_equil_model", &T::init_equil_model)
			.def("calc_equil", &T::py_calc_equil)
			.def("get_pressure", &T::get_pressure)
			.def("get_saturation", &T::get_saturation)
			.def("get_equil_region_data", &T::get_equil_region_data)
			.def("get_equil_data", &T::get_equil_data)
            .def("get_pbvd_region_data", &T::get_pbvd_region_data)
            .def("get_pbvd_data", &T::get_pbvd_data)
            .def("get_rsvd_region_data", &T::get_rsvd_region_data)
            .def("get_rsvd_data", &T::get_rsvd_data)
			.def("get_n_equil_regions", &T::get_n_equil_regions)
		PY_EXPORTER_END;

		void
		py_export_equil_model ()
		{
			base_exporter <equil_model_iface, equil_model_depth_exporter>::export_class ("equil_model_iface");
      class_exporter <equil_model_depth, equil_model_iface, equil_model_depth_exporter>::export_class ("equil_model_depth");

      // register equil to/from str serialization
      def("serialize_to_str", &blue_sky::serialize_to_str< equil_model_iface >);
      def("serialize_from_str", &blue_sky::serialize_from_str< equil_model_iface >);
		}
	}
}

#endif
