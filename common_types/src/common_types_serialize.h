/// @file common_types_serialize.h
/// @brief Various common types serialization
/// @author uentity
/// @version 
/// @date 27.01.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef COMMON_TYPES_SERIALIZE_ZD0U0FWQ
#define COMMON_TYPES_SERIALIZE_ZD0U0FWQ

#include "bs_serialize.h"

#include "table.h"
#include "prop.h"
#include "vartype_table.h"

BLUE_SKY_CLASS_SRZ_FCN_DECL(serialize, blue_sky::table)
BLUE_SKY_CLASS_SRZ_FCN_DECL(serialize, blue_sky::prop)
BLUE_SKY_CLASS_SRZ_FCN_DECL_T(serialize, blue_sky::vartype_table, 1)

// for GUID use std boost::serialization macro, cause we declare it for interface
//BOOST_CLASS_EXPORT_KEY2(blue_sky::table_iface, "table_iface")
//BLUE_SKY_TYPE_SERIALIZE_DECL_NOGUID(blue_sky::table_iface)
BLUE_SKY_TYPE_SERIALIZE_DECL(blue_sky::table)

//BOOST_CLASS_EXPORT_KEY2(blue_sky::prop_iface, "prop_iface")
//BLUE_SKY_TYPE_SERIALIZE_DECL_NOGUID(blue_sky::prop_iface)
BLUE_SKY_TYPE_SERIALIZE_DECL(blue_sky::prop)

// vartype_table_iface defined only for t_float
//BOOST_CLASS_EXPORT_KEY2(blue_sky::vartype_table_iface< t_float >, "vartype_table_iface_float")
//BLUE_SKY_TYPE_SERIALIZE_DECL_T(blue_sky::vartype_table_iface, 1)
BLUE_SKY_TYPE_SERIALIZE_DECL_T(blue_sky::vartype_table, 1)
BLUE_SKY_TYPE_SERIALIZE_GUID_T(blue_sky::vartype_table, t_float)

#endif /* end of include guard: COMMON_TYPES_SERIALIZE_ZD0U0FWQ */

