/// @file scal_3p_serialize.cpp
/// @brief scal_3p serialization implementation
/// @author uentity
/// @version 
/// @date 26.01.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "bs_scal_stdafx.h"

#include "conf.h"
#include "scal_3p.h"
#include "py_scal_wrapper.h"
#include "bs_serialize.h"
#include "scal_3p_impl.h"

#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/collection_size_type.hpp>

using namespace blue_sky;
namespace boser = boost::serialization;

/*-----------------------------------------------------------------
 * serialize some PODs
 *----------------------------------------------------------------*/
namespace scal_dp = blue_sky::scal::data_placement;

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, scal_dp::scal_placement_info)
	ar & t.sp_step;
	ar & t.so_step;
	ar & t.krp_step;
	ar & t.krop_step;
	ar & t.pcp_step;

	ar & t.sp_offset;
	ar & t.so_offset;
	ar & t.krp_offset;
	ar & t.krop_offset;
	ar & t.pcp_offset;

	ar & t.type;
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, scal_dp::scale_array_placement_info)
	ar & t.socr_step;
	ar & t.scr_step;
	ar & t.su_step;
	ar & t.sl_step;
	ar & t.pcp_step;

	ar & t.socr_offset;
	ar & t.scr_offset;
	ar & t.su_offset;
	ar & t.sl_offset;
	ar & t.pcp_offset;
	ar & t.size;
BLUE_SKY_CLASS_SRZ_FCN_END

/*-----------------------------------------------------------------
 * serialize scal_region_info
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN_T(serialize, scal_region_info, 1)
	ar & t.So_count;
	ar & t.Sp_count;
	ar & t.Krp_min_greater_zero;
	ar & t.Krop_min_greater_zero;
	ar & t.so_offset;
	ar & t.sp_offset;
	ar & t.spr;
	ar & t.sorp;
	ar & t.kpr;
	ar & t.krorp;
	ar & t.pcp_max;
BLUE_SKY_CLASS_SRZ_FCN_END

/*-----------------------------------------------------------------
 * serialize data_vector
 *----------------------------------------------------------------*/
//BLUE_SKY_CLASS_SRZ_FCN_BEGIN_T(save_construct_data, data_vector, 1)
//	boser::collection_size_type sz(t->size()), step(t->step());
//	ar << sz; ar << step;
//	// data
//	ar << boser::make_array(&(*t)[0], sz);
//BLUE_SKY_CLASS_SRZ_FCN_END
//
//BLUE_SKY_CLASS_SRZ_FCN_BEGIN_T(load_construct_data, data_vector, 1)
//	boser::collection_size_type sz, step;
//	ar >> sz; ar >> step;
//	// data
//	ar << boser::make_array(&(*t)[0], sz);
//BLUE_SKY_CLASS_SRZ_FCN_END

/*-----------------------------------------------------------------
 * serialize scal_region
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, scal_region)
	ar & t.Sp;
	ar & t.So;
	ar & t.Krp;
	ar & t.Krop;
	ar & t.Pcp;
BLUE_SKY_CLASS_SRZ_FCN_END

/*-----------------------------------------------------------------
 * serialize scal_2p_data_holder
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(save, scal_2p_data_holder)
	ar << t.placement_info_;
	ar << t.data_;
	ar << t.region_;
	//a<< & t.region_2_;
	ar << t.scal_table_array;
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(load, scal_2p_data_holder)
	ar >> t.placement_info_;
	ar >> t.data_;
	ar >> t.region_;
	//a>> & t.region_2_;
	ar >> t.scal_table_array;

	// fill region_2_
	t.region_2_.clear();
	t.init_regions();
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, scal_2p_data_holder)
	// register conversion to interface
	boser::bs_void_cast_register(
		static_cast< scal_2p_data_holder* >(NULL),
		static_cast< scal_2p_data_holder_iface* >(NULL)
	);
	// split to save/load
	boser::split_free(ar, t, version);
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_TYPE_SERIALIZE_DECL(scal_2p_data_holder)
BLUE_SKY_TYPE_SERIALIZE_IMPL(scal_2p_data_holder)

/*-----------------------------------------------------------------
 * serialize scale_array_holder
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, scale_array_holder)
	// register conversion to interface
	boser::bs_void_cast_register(
		static_cast< scale_array_holder* >(NULL),
		static_cast< scale_array_holder_iface* >(NULL)
	);
	// serialize data
	ar & t.data_pool;
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_TYPE_SERIALIZE_DECL(scale_array_holder)
BLUE_SKY_TYPE_SERIALIZE_IMPL(scale_array_holder)

/*-----------------------------------------------------------------
 * serialize jfunction
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, jfunction)
	ar & t.st_phase;
	ar & t.alpha;
	ar & t.beta;
	ar & t.is_valid;
	ar & t.plane_a;
	ar & t.plane_b;
	ar & t.perm_type;
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_TYPE_SERIALIZE_DECL(jfunction)
BLUE_SKY_TYPE_SERIALIZE_IMPL(jfunction)

/*-----------------------------------------------------------------
 * serialize scal_3p_impl
 *----------------------------------------------------------------*/
//BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, scal_3p_impl)
//	// register conversion to interface
//	boser::bs_void_cast_register(
//		static_cast< scal_3p_impl* >(NULL),
//		static_cast< scal_3p::scal_3p_impl_base* >(NULL)
//	);
//
//	ar & t.is_w & t.is_g & t.is_o;
//	ar & t.rpo_model & t.n_phases;
//
//	ar & t.water_scale;
//	ar & t.gas_scale;
//	ar & t.water_data;
//	ar & t.gas_data;
//	ar & t.water_jfunc;
//	ar & t.gas_jfunc;
//
//	ar & t.i_w & t.i_g & t.i_o;
//	ar & t.i_w_w & t.i_w_g & t.i_w_o;
//	ar & t.i_g_w & t.i_g_g & t.i_g_o;
//	ar & t.i_o_w & t.i_o_g & t.i_o_o;
//
//	ar & t.i_s_w & t.i_s_g;
//	ar & t.is_scalecrs;
//BLUE_SKY_CLASS_SRZ_FCN_END
//
//BOOST_CLASS_EXPORT(scal_3p_impl)

/*-----------------------------------------------------------------
 * serialize scal_3p
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(save, scal_3p)
	// save data
	ar << t.water_data;
	ar << t.gas_data;
	ar << t.water_scale;
	ar << t.gas_scale;
	ar << t.water_jfunc;
	ar << t.gas_jfunc;

	ar << t.n_scal_regions;
	ar << t.water_input_table;
	ar << t.gas_input_table;
	ar << t.oil_input_table;

	ar << t.is_gas;
	ar << t.is_oil;
	ar << t.is_water;

	// should we save impl_?
	bool save_impl = (t.impl_ != NULL);
	ar << save_impl;
	if(!save_impl) return;

	// save internal impl_ variables
	scal_3p_impl& tt = *static_cast< scal_3p_impl* >(t.impl_);
	ar << tt.is_w << tt.is_g << tt.is_o;
	ar << tt.rpo_model;

	ar << tt.water_scale;
	ar << tt.gas_scale;
	ar << tt.water_data;
	ar << tt.gas_data;
	ar << tt.water_jfunc;
	ar << tt.gas_jfunc;

	ar << tt.i_w << tt.i_g << tt.i_o;
	ar << tt.i_w_w << tt.i_w_g << tt.i_w_o;
	ar << tt.i_g_w << tt.i_g_g << tt.i_g_o;
	ar << tt.i_o_w << tt.i_o_g << tt.i_o_o;

	ar << tt.i_s_w << tt.i_s_g;
	ar << tt.is_scalecrs;
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(load, scal_3p)
	// load data
	ar >> t.water_data;
	ar >> t.gas_data;
	ar >> t.water_scale;
	ar >> t.gas_scale;
	ar >> t.water_jfunc;
	ar >> t.gas_jfunc;

	ar >> t.n_scal_regions;
	ar >> t.water_input_table;
	ar >> t.gas_input_table;
	ar >> t.oil_input_table;

	ar >> t.is_gas;
	ar >> t.is_oil;
	ar >> t.is_water;

	// should we load impl_?
	bool load_impl;
	ar >> load_impl;
	if(!load_impl) return;

	// load internal impl_ variables
	t.impl_ = new scal_3p_impl;
	scal_3p_impl& tt = *static_cast< scal_3p_impl* >(t.impl_);
	ar >> tt.is_w >> tt.is_g >> tt.is_o;
	ar >> tt.rpo_model;

	ar >> tt.water_scale;
	ar >> tt.gas_scale;
	ar >> tt.water_data;
	ar >> tt.gas_data;
	ar >> tt.water_jfunc;
	ar >> tt.gas_jfunc;

	ar >> tt.i_w >> tt.i_g >> tt.i_o;
	ar >> tt.i_w_w >> tt.i_w_g >> tt.i_w_o;
	ar >> tt.i_g_w >> tt.i_g_g >> tt.i_g_o;
	ar >> tt.i_o_w >> tt.i_o_g >> tt.i_o_o;

	ar >> tt.i_s_w >> tt.i_s_g;
	ar >> tt.is_scalecrs;
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, scal_3p)
	// register conversion to interface
	boser::bs_void_cast_register(
		static_cast< scal_3p* >(NULL),
		static_cast< scal_3p_iface* >(NULL)
	);
	// split
	boser::split_free(ar, t, version);
BLUE_SKY_CLASS_SRZ_FCN_END

BOOST_SERIALIZATION_ASSUME_ABSTRACT(scal_3p_iface)
BLUE_SKY_TYPE_SERIALIZE_DECL(scal_3p)
BLUE_SKY_TYPE_SERIALIZE_IMPL(scal_3p)

