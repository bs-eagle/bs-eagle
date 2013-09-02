/**
 *       \file  py_pvt.cpp
 *      \brief  python wrappers for pvt_
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  08.05.2008
 *  \copyright  This source code is released under the terms of
 *              the BSD License. See LICENSE for more details.
 * */
#include "bs_pvt_stdafx.h"

#include "py_pvt.h"

#ifdef BSPY_EXPORTING_PLUGIN

#include "pvt_base.h"
#include "pvt_dead_oil.h"
#include "pvt_gas.h"
#include "pvt_oil.h"
#include "pvt_water.h"
#include "pvt_dummy.h"
#include "pvt_3p.h"
#include "bs_serialize.h"

#include "export_python_wrapper.h"
#include "py_list_converter.h"

namespace blue_sky  {
namespace python    {

  PY_EXPORTER (pvt_exporter, default_exporter)
    .add_property ("p_step",          &T::get_p_step)
    .add_property ("surface_density", &T::get_surface_density)
    .def ("get_gor_for_pressure",     &T::get_gor_for_pressure)
    .def ("print",                    &T::print)
  PY_EXPORTER_END;


    PY_EXPORTER (pvt_dummy_exporter, default_exporter)
      .def ("get_table", &T::get_table)
      .def ("get_pvt_type", &T::get_pvt_type)
    PY_EXPORTER_END;

	PY_EXPORTER (pvt_3p_exporter, default_exporter)
	  .def ("init_pvt_arrays", &T::init_pvt_arrays)
	  .def ("init_pvt_calc_data", &T::init_pvt_calc_data)
	  .def ("get_table", &T::get_table)
	  .def ("get_n_pvt_regions", &T::get_n_pvt_regions)
	  .def ("get_density_table", &T::get_density_table)
	  .def ("get_tables_list", &T::get_tables_list)
	  .def ("get_density_list", &T::get_density_list)
	PY_EXPORTER_END;

  void
  py_export_pvt ()
  {
    using namespace boost::python;

    base_exporter<pvt_base, pvt_exporter>::export_class ("pvt_base");
    class_exporter<pvt_oil, pvt_base, pvt_exporter>::export_class  ("pvt_oil");
    class_exporter<pvt_gas, pvt_base, pvt_exporter>::export_class ("pvt_gas");
    class_exporter<pvt_water, pvt_base, pvt_exporter>::export_class ("pvt_water");
    class_exporter<pvt_dead_oil, pvt_base, pvt_exporter>::export_class ("pvt_dead_oil");

    base_exporter <pvt_dummy_iface, pvt_dummy_exporter>::export_class ("pvt_dummy_iface");
    class_exporter<pvt_dummy, pvt_dummy_iface, pvt_dummy_exporter>::export_class ("pvt_dummy");

	base_exporter <pvt_3p_iface, pvt_3p_exporter>::export_class ("pvt_3p_iface");
	class_exporter<pvt_3p, pvt_3p_iface, pvt_3p_exporter>::export_class ("pvt_3p");

	  // register vector of type descriptors <-> Python list converters
	  typedef bspy_converter< list_traits< std::list<smart_ptr<table_iface> > > > spv_table_list_converter;
	  spv_table_list_converter::register_from_py ();
	  spv_table_list_converter::register_to_py ();

    // register pvt to/from str serialization
	def("serialize_to_str", &blue_sky::serialize_to_str< pvt_3p_iface >);
	def("serialize_from_str", &blue_sky::serialize_from_str< pvt_3p_iface >);
  }

} // namespace python
} // namespace blue_sky

#endif  // #ifdef BSPY_EXPORTING_PLUGIN

