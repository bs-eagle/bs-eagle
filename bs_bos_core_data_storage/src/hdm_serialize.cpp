/// @file hdm_serialize.cpp
/// @brief Serialization implementation for hdm
/// @author uentity
/// @version 
/// @date 20.01.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "bs_bos_core_data_storage_stdafx.h"

#include "hdm_serialize.h"
#include "sql_well_serialize.h"
#include "locale_keeper_serialize.h"
#include "convert_units_serialize.h"
#include "idata_serialize.h"

using namespace blue_sky;
namespace boser = boost::serialization;

/*-----------------------------------------------------------------
 * save hdm
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(save, hdm)
	// save idata
	ar << t.data;
	// save number of pvt regions
	ar << (const t_long&)t.pvt_3p_->get_n_pvt_regions();
	// save number of scal regions
	ar << (const t_long&)t.scal_3p_->get_n_scal_regions();
	// save number of equil regions
	ar << (const t_long&)t.equil_model_->get_n_equil_regions();
	// save sql_well
	ar << t.well_pool_;
	// locale_keeper
	ar << t.lkeeper;
	// constants
	ar << t.ph_const;
BLUE_SKY_CLASS_SRZ_FCN_END

/*-----------------------------------------------------------------
 * load hdm
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(load, hdm)
	// load idata (including h5_pool)
	ar >> t.data;

	// load reg nums
	t_long n_pvt_r, n_scal_r, n_equil_r;
	ar >> n_pvt_r; ar >> n_scal_r; ar >> n_equil_r;
	// init models using numbers read
	t.init_fluids(n_pvt_r, n_scal_r);
	t.init_equil(n_equil_r);

	// 
BLUE_SKY_CLASS_SRZ_FCN_END

/*-----------------------------------------------------------------
 * register conversion to interface + serialize via split
 *----------------------------------------------------------------*/
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

