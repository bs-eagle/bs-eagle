/// @file locale_keeper_serialize.cpp
/// @brief Serialization implementation for locale_keeper
/// @author uentity
/// @version 
/// @date 20.01.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#if defined(BSPY_EXPORTING_PLUGIN) && defined(UNIX)
// supress gcc warnings
#include <boost/python/detail/wrap_python.hpp>
#endif

#include "locale_keeper_serialize.h"

#include <boost/archive/polymorphic_iarchive.hpp>
#include <boost/archive/polymorphic_oarchive.hpp>
#include <boost/serialization/string.hpp>
#include <string.h>

using namespace blue_sky;

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(save_construct_data, locale_keeper)
	std::string loc(t->locale_);
	ar << loc;
	ar << t->category_;
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(load_construct_data, locale_keeper)
	std::string loc(t->locale_);
	int category;
	ar >> loc;
	ar >> category;
	::new(t) locale_keeper(loc.c_str(), category);
BLUE_SKY_CLASS_SRZ_FCN_END

// empty serialize
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, locale_keeper)
BLUE_SKY_CLASS_SRZ_FCN_END

BOOST_CLASS_EXPORT_IMPLEMENT(locale_keeper)

