/// @file idata_serialize.cpp
/// @brief idata serialization implementation
/// @author uentity
/// @version 
/// @date 26.01.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "bs_bos_core_data_storage_stdafx.h"

#include "auto_value.h"
#include "idata_serialize.h"
#include "h5_pool_serialize.h"
#include "convert_units_serialize.h"
#include "rocktab_table_serialize.h"
//#include "common_types_serialize.h"

#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

using namespace blue_sky;
namespace boser = boost::serialization;

/*-----------------------------------------------------------------
 * serialize pvt_info & scal_info
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, idata::pvt_info)
	ar & t.main_data_;
	ar & (bool&)t.has_density_;
	ar & (t_float&)t.density_;
	ar & (t_float&)t.molar_density_;
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, idata::scal_info)
	ar & t.main_data_;
BLUE_SKY_CLASS_SRZ_FCN_END

BOOST_CLASS_EXPORT_IMPLEMENT(idata::pvt_info)
BOOST_CLASS_EXPORT_IMPLEMENT(idata::scal_info)

/*-----------------------------------------------------------------
 * serialize val_vs_depth
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, val_vs_depth)
	ar & t.tdepth();
	ar & t.tvalues();
BLUE_SKY_CLASS_SRZ_FCN_END

/*-----------------------------------------------------------------
 * serialize idata
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, idata)
	ar & t.props;
	ar & t.h5_pool;
	ar & t.input_units_converter;
	ar & t.output_units_converter;
	ar & t.rocktab;

	ar & t.pvto; ar & t.pvtdo; ar & t.pvtg; ar & t.pvtw;
	ar & t.swof; ar & t.sgof; ar & t.swfn; ar & t.sgfn; ar & t.sof2; ar & t.sof3;

	ar & t.equil_regions;
	ar & t.rock;
	ar & t.equil;
	ar & t.p_ref;

	ar & t.prvd;
	ar & t.rsvd;
	ar & t.pbvd;
BLUE_SKY_CLASS_SRZ_FCN_END

// instantiate serialization code
BLUE_SKY_TYPE_SERIALIZE_IMPL(idata)

