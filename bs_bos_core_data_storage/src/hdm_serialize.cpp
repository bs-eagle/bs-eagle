/// @file hdm_serialize.cpp
/// @brief Serialization implementation for hdm
/// @author uentity
/// @version 
/// @date 20.01.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "bs_bos_core_data_storage_stdafx.h"

#include "hdm_serialize.h"
//#include "sql_well_serialize.h"

using namespace blue_sky;
namespace boser = boost::serialization;

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(save, hdm)
	// save number of pvt regions
	ar << (const t_long&)t.pvt_3p_->get_n_pvt_regions();
	// save number of scal regions
	ar << (const t_long&)t.scal_3p_->get_n_scal_regions();
	// save number of equil regions
	ar << (const t_long&)t.equil_model_->get_n_equil_regions();
	// save sql_well
	ar << t.well_pool_;
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(load, hdm)

BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, hdm)
	// register conversion to base iface
	boser::bs_void_cast_register(
		static_cast< hdm* >(NULL),
		static_cast< hdm_iface* >(NULL)
	);
	// split into save/load
	boser::split_free(ar, t, version);
BLUE_SKY_CLASS_SRZ_FCN_END

// instantiate serialization code
BLUE_SKY_TYPE_SERIALIZE_IMPL(hdm)

