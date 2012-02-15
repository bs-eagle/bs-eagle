/// @file convert_units_serialize.cpp
/// @brief 
/// @author uentity
/// @version 
/// @date 19.01.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "bs_bos_core_data_storage_stdafx.h"

#include "convert_units_serialize.h"

using namespace blue_sky;

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, physical_constants)
	ar & t.darcy_constant;
	ar & t.gravity_constant;
	ar & t.atmospheric_pressure;
	ar & t.jfunc_constant;
	ar & t.default_injection_bhp_limit;
	ar & t.default_production_bhp_limit;
	ar & t.default_minimal_average_pressure;
	ar & t.default_liquid_rate_limit;
	ar & t.default_gas_rate_limit;
	ar & t.gas_constant;
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, convert_units)
	ar & t.input_constants;
	ar & t.output_constants;
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(save_construct_data, convert_units)
	ar << (const int&)t->get_input_units();
	ar << (const int&)t->get_output_units();
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(load_construct_data, convert_units)
	int units_in, units_out;
	ar >> units_in; ar >> units_out;
	::new(t) convert_units(units_in, units_out);
BLUE_SKY_CLASS_SRZ_FCN_END

BOOST_CLASS_EXPORT_IMPLEMENT(physical_constants)
BLUE_SKY_CLASS_SERIALIZE_INST(physical_constants)
BOOST_CLASS_EXPORT_IMPLEMENT(convert_units)
BLUE_SKY_CLASS_SERIALIZE_INST(convert_units)

