/// @file h5_pool_serialize.cpp
/// @brief h5_pool serialization implementation
/// @author uentity
/// @version 
/// @date 24.01.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

//#include "bs_bos_core_data_storage_stdafx.h"

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
	if(do_open) {
		t.file_id = H5Fopen(t.fname.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
		if(t.file_id < 0) // try to create file
			t.open_file(t.fname.c_str(), t.path.c_str());
		else {
			// try to open group
			t.group_id = H5Gopen(t.file_id, t.path.c_str());
			if(t.group_id < 0) {
				// try to create group
				t.group_id = H5Gcreate(t.file_id, t.path.c_str (), -1);
				if(t.group_id < 0)
					bs_throw_exception(
						boost::format ("h5_pool_serialize: Can't create group: %s") % t.path
					);
			}
			// init
			t.fill_map();
		}
	}
	ar >> boser::make_array(&t.pool_dims, 3);
	ar >> t.n_pool_dims;
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, h5_pool)
	// register conversion to base interface
	boser::bs_void_cast_register< h5_pool, h5_pool_iface >(
		static_cast< h5_pool* >(NULL),
		static_cast< h5_pool_iface* >(NULL)
	);
	boser::split_free(ar, t, version);
BLUE_SKY_CLASS_SRZ_FCN_END

BOOST_SERIALIZATION_ASSUME_ABSTRACT(h5_pool_iface)
BLUE_SKY_TYPE_SERIALIZE_IMPL(h5_pool)

