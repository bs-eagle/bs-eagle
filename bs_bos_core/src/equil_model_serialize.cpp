/// @file equil_model_serialize.cpp
/// @brief equil_model serialization implementation
/// @author uentity
/// @version 
/// @date 30.01.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "stdafx.h"
#include "bs_serialize.h"
#include "equil_model_depth.h"

#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

using namespace blue_sky;
namespace boser = boost::serialization;

/*-----------------------------------------------------------------
 * serialize equil_model_depth
 *----------------------------------------------------------------*/

BOOST_CLASS_VERSION(equil_model_depth, 1)

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, equil_model_depth)
	// register conversion to base iface
	boser::bs_void_cast_register(
		static_cast< equil_model_depth* >(NULL),
		static_cast< equil_model_iface* >(NULL)
	);
	// serialize data
	ar & t.n_equil_regions;
	ar & t.equil_data;
	if(version > 0) {
        ar & t.pbvd_data;
        ar & t.rsvd_data;
    }
	ar & t.pressure;
	ar & t.saturation;
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_TYPE_SERIALIZE_DECL(equil_model_depth)
BLUE_SKY_TYPE_SERIALIZE_IMPL(equil_model_depth)

