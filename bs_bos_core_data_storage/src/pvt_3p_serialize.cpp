/// @file pvt_3p_serialize.cpp
/// @brief pvt_3p serialization implementation
/// @author uentity
/// @version 
/// @date 26.01.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "bs_bos_core_data_storage_stdafx.h"
#include "bs_pvt_stdafx.h"

#include "pvt_3p_serialize.h"
#include "table_iface.h"

#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

using namespace blue_sky;
namespace boser = boost::serialization;

/*-----------------------------------------------------------------
 * serialize table_iface
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(save, table_iface)
	// dump all info to string and save it
	ar << (const std::string&)t.to_str();
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(load, table_iface)
	// dump all info to string and save it
	std::string prop_data;
	ar >> prop_data;
	t.from_str(prop_data);
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SERIALIZE_SPLIT(table_iface)

// define GUID and implementation of prop_iface serialization
// for GUID use std boost::serialization macro, cause we declare it for interface
// also instantiate code using _BYNAME macro
BOOST_CLASS_EXPORT_KEY2(table_iface, "table")
BLUE_SKY_TYPE_SERIALIZE_DECL_NOGUID(table_iface)
BLUE_SKY_TYPE_SERIALIZE_IMPL_BYNAME(table_iface, "table")

/*-----------------------------------------------------------------
 * serialize pvt_base
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, pvt_base)
	ar & t.p_step;
	ar & t.surface_density;
	ar & t.molar_density;
	ar & t.init_dependent;

	ar & t.pvt_input_props;
	ar & t.pvt_props_table;
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_TYPE_SERIALIZE_DECL(pvt_base)
BLUE_SKY_TYPE_SERIALIZE_IMPL(pvt_base)

/*-----------------------------------------------------------------
 * serialize pvt_dead_oil
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, pvt_dead_oil)
	// just redirect to base class
	boser::bs_base_object< pvt_base >(t);
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_TYPE_SERIALIZE_DECL(pvt_dead_oil)
BLUE_SKY_TYPE_SERIALIZE_IMPL(pvt_dead_oil)

/*-----------------------------------------------------------------
 * serialize pvt_gas
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, pvt_gas)
	// just redirect to base class
	boser::bs_base_object< pvt_base >(t);
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_TYPE_SERIALIZE_DECL(pvt_gas)
BLUE_SKY_TYPE_SERIALIZE_IMPL(pvt_gas)

/*-----------------------------------------------------------------
 * serialize pvt_oil
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, pvt_oil)
	// just redirect to base class
	boser::bs_base_object< pvt_base >(t);
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_TYPE_SERIALIZE_DECL(pvt_oil)
BLUE_SKY_TYPE_SERIALIZE_IMPL(pvt_oil)

/*-----------------------------------------------------------------
 * serialize pvt_water
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, pvt_water)
	// just redirect to base class
	boser::bs_base_object< pvt_base >(t);
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_TYPE_SERIALIZE_DECL(pvt_water)
BLUE_SKY_TYPE_SERIALIZE_IMPL(pvt_water)

/*-----------------------------------------------------------------
 * serialize pvt_3p
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, pvt_3p)
	// register conversion to interface
	boser::bs_void_cast_register(
		static_cast< pvt_3p* >(NULL),
		static_cast< pvt_3p_iface* >(NULL)
	);
	ar & t.n_pvt_regions;
	ar & t.pvt_oil_array;
	ar & t.pvt_water_array;
	ar & t.pvt_gas_array;

	ar & t.density_table;
	ar & t.density;
BLUE_SKY_CLASS_SRZ_FCN_END

BOOST_SERIALIZATION_ASSUME_ABSTRACT(pvt_3p_iface)

BLUE_SKY_TYPE_SERIALIZE_IMPL(pvt_3p)

