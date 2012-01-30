/// @file hdm_serialize.cpp
/// @brief Serialization implementation for hdm
/// @author uentity
/// @version 
/// @date 20.01.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "bs_bos_core_data_storage_stdafx.h"

#include "hdm_serialize.h"
#include "common_types_serialize.h"
#include "sql_well_serialize.h"
#include "locale_keeper_serialize.h"
#include "convert_units_serialize.h"
#include "idata_serialize.h"
#include "scal_3p_serialize.h"
#include "pvt_3p_serialize.h"

#include "equil_model_depth.h"

#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

using namespace blue_sky;
namespace boser = boost::serialization;

/*-----------------------------------------------------------------
 * serialize equil_model_depth
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, equil_model_depth)
	// register conversion to base iface
	boser::bs_void_cast_register(
		static_cast< equil_model_depth* >(NULL),
		static_cast< equil_model_iface* >(NULL)
	);
	// serialize data
	ar & t.n_equil_regions;
	ar & t.equil_data;
	ar & t.pressure;
	ar & t.saturation;
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_TYPE_SERIALIZE_DECL(equil_model_depth)
BLUE_SKY_TYPE_SERIALIZE_IMPL(equil_model_depth)

/*-----------------------------------------------------------------
 * save hdm
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(save, hdm)
	// save some variables
	ar << t.data;
	ar << t.scal_3p_;
	ar << t.pvt_3p_;
	ar << t.equil_model_;

	// save sql_well
	ar << t.well_pool_;
	// locale_keeper
	ar << t.lkeeper;
	// constants
	ar << t.ph_const;

	// save number of pvt regions
	ar << (const t_long&)t.pvt_3p_->get_n_pvt_regions();
	// save number of scal regions
	ar << (const t_long&)t.scal_3p_->get_n_scal_regions();
	// save number of equil regions
	ar << (const t_long&)t.equil_model_->get_n_equil_regions();
BLUE_SKY_CLASS_SRZ_FCN_END

/*-----------------------------------------------------------------
 * load hdm
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(load, hdm)
	typedef smart_ptr< type > sp_hdm;
	// load idata (including h5_pool)
	ar >> t.data;
	ar >> t.scal_3p_;
	ar >> t.pvt_3p_;
	ar >> t.equil_model_;

	ar >> t.well_pool_;
	ar >> t.lkeeper;
	ar >> t.ph_const;

	// load reg nums
	t_long n_pvt_r, n_scal_r, n_equil_r;
	ar >> n_pvt_r; ar >> n_scal_r; ar >> n_equil_r;
	// init models using numbers read
	//t.init_fluids(n_pvt_r, n_scal_r);
	//t.init_equil(n_equil_r);

	// init mesh
	t.mesh->init_props(sp_hdm(&t));
	// init keyword_manager
	keyword_params kp;
	kp.hdm = &t;
	t.km->init(sp_hdm(&t));
	switch (t.data->props->get_i("mesh")) {
	case 0:
		t.km->handle_keyword_reactor ("MESH_IJK", kp);
		break;
	case 1:
		t.km->handle_keyword_reactor ("MESH_GRDECL", kp);
		break;
	default:
			bs_throw_exception ("init: wrong mesh choice");
	}
	switch (t.data->props->get_i("init")) {
	case 0:
		t.km->handle_keyword_reactor ("EXPLICIT_MODEL", kp);
		break;
	case 1:
		t.km->handle_keyword_reactor ("EQUIL_MODEL", kp);
		break;
	default:
		bs_throw_exception ("init: wrong mesh choice");
	}
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

