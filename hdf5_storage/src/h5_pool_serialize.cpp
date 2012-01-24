/// @file h5_pool_serialize.cpp
/// @brief h5_pool serialization implementation
/// @author uentity
/// @version 
/// @date 24.01.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "h5_pool_serialize.h"
#include <boost/serialization/array.hpp>

using namespace blue_sky;
namespace boser = boost::serialization;

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(save, h5_pool)
	bool is_open = (t.file_id != 0);
	ar << is_open;
	ar << t.fname;
	ar << t.path;
	ar << boser::make_array(&t.pool_dims, 3);
	ar << t.n_pool_dims;
	// save flag if h5 is open
	// flush buffers
	const_cast< h5_pool& >(t).flush();
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(load, h5_pool)
	bool do_open;
	ar >> do_open;
	ar >> t.fname;
	ar >> t.path;
	if(do_open)
		t.open_file(t.fname, t.path);
	ar >> boser::make_array(&t.pool_dims, 3);
	ar >> t.n_pool_dims;
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SERIALIZE_SPLIT(h5_pool)

BLUE_SKY_TYPE_SERIALIZE_IMPL(h5_pool)

