/// @file hdm_serialize.h
/// @brief Serialization declaration for hdm
/// @author uentity
/// @version 
/// @date 20.01.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef HDM_SERIALIZE_B9ETGFM9
#define HDM_SERIALIZE_B9ETGFM9

#include "bs_serialize.h"
#include "hdm.h"

BLUE_SKY_CLASS_SRZ_FCN_DECL(serialize, blue_sky::hdm)

BLUE_SKY_TYPE_SERIALIZE_DECL(blue_sky::hdm)

/*-----------------------------------------------------------------
 * declare hdm save and load functions
 * that really will be invoked by clients
 *----------------------------------------------------------------*/

namespace blue_sky {
////////////////////////////////////////////////////////////////////
// note that prj_name should not include extension
//
BS_API_PLUGIN void hdm_serialize_save(
	smart_ptr< hdm > t,
	const std::wstring& prj_path,
	const std::wstring& prj_name,
	const std::wstring& deep_copy_suffix = L""  // h5 & well pool files should be copied using this suffix
);

BS_API_PLUGIN smart_ptr< hdm > hdm_serialize_load(
	const std::wstring& prj_path,
	const std::wstring& prj_name
);

std::string hdm_serialize_to_str(
	smart_ptr< hdm > t,
	const std::wstring& prj_path,
	const std::wstring& prj_name,
	const std::wstring& deep_copy_suffix = L"" // h5 & well pool files should be copied using this suffix
);

smart_ptr< hdm > hdm_serialize_from_str(
	const std::string& hdm_dump,
	const std::wstring& prj_path,
	const std::wstring& prj_name
);

} /* blue_sky */

#endif /* end of include guard: HDM_SERIALIZE_B9ETGFM9 */

