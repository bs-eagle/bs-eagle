/**
 *       \file  py_scal_wrapper.cpp
 *      \brief  Exports scal_3p (+) to python
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  14.12.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "bs_scal_stdafx.h"

#include "py_scal_wrapper.h"
#include "scal_save_data.h"
#include "scal_2p_dummy.h"

#include "export_python_wrapper.h"

#ifdef BSPY_EXPORTING_PLUGIN
namespace blue_sky
  {
  namespace python
    {

    PY_EXPORTER (scale_array_holder_exporter, default_exporter)
    PY_EXPORTER_END;

    PY_EXPORTER (scal_data_holder_exporter, default_exporter)
    PY_EXPORTER_END;

    PY_EXPORTER (jfunction_exporter, default_exporter)
      .add_property ("st_phase",  &T::st_phase)
      .add_property ("alpha",     &T::alpha)
      .add_property ("beta",      &T::beta)
      .add_property ("valid",     &T::valid)
	  .def ("init", &T::init)
	  .def ("init_ab", &T::init_ab)
    PY_EXPORTER_END;

    PY_EXPORTER (scal_3p_exporter, default_exporter)
      .add_property ("water_scale",      &T::get_water_scale)
      .add_property ("gas_scale",        &T::get_gas_scale)
      .add_property ("water_data",       &T::get_water_data)
      .add_property ("gas_data",         &T::get_gas_data)
      .add_property ("water_jfunction",  &T::get_water_jfunction, &T::set_water_jfunction)
      .add_property ("gas_jfunction",    &T::get_gas_jfunction, &T::set_gas_jfunction)
    .def ("init_tables", &T:: init_scal_input_table_arrays)
	  .def ("init_from_scal", &T::init_from_scal)
	  .def ("get_table", &T::get_tables)
    PY_EXPORTER_END;
    
    PY_EXPORTER (scal_dummy_exporter, default_exporter)
      .def ("get_table", &T::py_get_table)
    PY_EXPORTER_END;

    void py_export_scal ()
    {
      using namespace boost::python;

      enum_ <PHASE_ENUM> ("phase_enum")
        .value ("water", PHASE_WATER)
        .value ("gas", PHASE_GAS)
        .value ("oil", PHASE_OIL)
        ;

      enum_ <RPO_MODEL_ENUM> ("rpo_model_enum")
        .value ("default_model", RPO_DEFAULT_MODEL)
        .value ("stone1_model", STONE1_MODEL)
        .value ("stone2_model", STONE2_MODEL)
        ;

      enum_ <JFUNC_PERM_TYPE_ENUM> ("jfunc_perm_type_enum")
        .value ("xy", JFUNC_PERM_XY)
        .value ("x", JFUNC_PERM_X)
        .value ("y", JFUNC_PERM_Y)
        .value ("z", JFUNC_PERM_Z)
        ;

 
      base_exporter<jfunction, jfunction_exporter>::export_class ("jfunction");
      base_exporter<scale_array_holder, scale_array_holder_exporter>::export_class ("scale_arrays");
      base_exporter<scal_2p_data_holder, scal_data_holder_exporter>::export_class ("scal_data");
      base_exporter<scal_3p_iface, default_exporter>::export_class ("scal_3p_iface");
      class_exporter<scal_3p, scal_3p_iface, scal_3p_exporter>::export_class ("scal_3p");
      
      base_exporter <scal_dummy_iface, scal_dummy_exporter>::export_class ("scal_dummy_iface");
      class_exporter<scal_2p_dummy, scal_dummy_iface, scal_dummy_exporter>::export_class ("scal_2p_dummy");
      class_exporter<scal_3p_dummy, scal_dummy_iface, scal_dummy_exporter>::export_class ("scal_3p_dummy");
      
    }
  }
}
#endif