/// @file rocktab_table_serialize.cpp
/// @brief bs_table & rocktab_table serialization implementation
/// @author uentity
/// @version 
/// @date 26.01.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "rocktab_table_serialize.h"

#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

namespace boser = boost::serialization;

/*-----------------------------------------------------------------
 * serialize bs_table
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(save, blue_sky::bs_table)
	// save rows num and col num
	ar << t.n_rows;
	const ulong col_num = ulong(t.columns.size());
	ar << col_num;

	// save data
	for(ulong i = 0; i < col_num; ++i) {
		ar << t.columns[i].name;
		ar << t.columns[i].data;
	}
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(load, blue_sky::bs_table)
	// restore rows num
	ar >> t.n_rows;

	// load data
	ulong col_num;
	ar >> col_num;
	t.columns.resize(col_num);
	for(ulong i = 0; i < t.columns.size(); ++i) {
		ar >> t.columns[i].name;
		ar >> t.columns[i].data;
	}
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SERIALIZE_SPLIT(blue_sky::bs_table)

BOOST_CLASS_EXPORT_IMPLEMENT(blue_sky::bs_table)
BLUE_SKY_CLASS_SERIALIZE_INST(blue_sky::bs_table)

/*-----------------------------------------------------------------
 * serialize rocktab_table
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, blue_sky::rocktab_table)
	// just invoke base class serialization
	ar & boser::base_object< bs_table >(t);
BLUE_SKY_CLASS_SRZ_FCN_END

BOOST_CLASS_EXPORT_IMPLEMENT(blue_sky::rocktab_table)
BLUE_SKY_CLASS_SERIALIZE_INST(blue_sky::rocktab_table)

