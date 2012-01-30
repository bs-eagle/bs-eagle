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

#include "table_iface.h"
#include "prop_iface.h"
#include "vartype_table_iface.h"

BLUE_SKY_CLASS_SRZ_FCN_DECL(serialize, blue_sky::table_iface)
BLUE_SKY_CLASS_SRZ_FCN_DECL(serialize, blue_sky::prop_iface)
BLUE_SKY_CLASS_SRZ_FCN_DECL_T(serialize, blue_sky::vartype_table_iface, 1)

// for GUID use std boost::serialization macro, cause we declare it for interface
BOOST_CLASS_EXPORT_KEY2(blue_sky::table_iface, "table")
BLUE_SKY_TYPE_SERIALIZE_DECL_NOGUID(blue_sky::table_iface)

BOOST_CLASS_EXPORT_KEY2(blue_sky::prop_iface, "prop")
BLUE_SKY_TYPE_SERIALIZE_DECL_NOGUID(blue_sky::prop_iface)

// vartype_table_iface defined only for t_float
BOOST_CLASS_EXPORT_KEY2(blue_sky::vartype_table_iface< t_float >, "float_var_table")
BLUE_SKY_TYPE_SERIALIZE_DECL_T(blue_sky::vartype_table_iface, 1)

#endif /* end of include guard: COMMON_TYPES_SERIALIZE_ZD0U0FWQ */

