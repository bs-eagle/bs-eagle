/// @file convert_units_serialize.h
/// @brief Serialization declarations for convert_units & physical_constants
/// @author uentity
/// @version 
/// @date 17.01.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef CONVERT_UNITS_SERIALIZE_D00HJLUP
#define CONVERT_UNITS_SERIALIZE_D00HJLUP

#include "bs_serialize.h"
#include "convert_units.h"

BLUE_SKY_CLASS_SRZ_FCN_DECL(serialize, blue_sky::physical_constants)

BLUE_SKY_CLASS_SRZ_FCN_DECL(serialize, blue_sky::convert_units)
BLUE_SKY_CLASS_SRZ_FCN_DECL(save_construct_data, blue_sky::convert_units)
BLUE_SKY_CLASS_SRZ_FCN_DECL(load_construct_data, blue_sky::convert_units)

BOOST_CLASS_EXPORT_KEY(blue_sky::physical_constants)
BOOST_CLASS_EXPORT_KEY(blue_sky::convert_units)

#endif /* end of include guard: CONVERT_UNITS_SERIALIZE_D00HJLUP */

