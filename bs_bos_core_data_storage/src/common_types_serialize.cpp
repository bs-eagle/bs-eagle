/// @file common_types_serialize.cpp
/// @brief Various common types serialization implementation
/// @author uentity
/// @version 
/// @date 27.01.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "common_types_serialize.h"

#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

using namespace blue_sky;
namespace boser = boost::serialization;

/*-----------------------------------------------------------------
 * serialize table_iface
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(save, table_iface)
	// dump all info to string and save it
	ar << (const std::string&)t.to_str();
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(load, table_iface)
	// dump all info to string and save it
	std::string prop_data;
	ar >> prop_data;
	t.from_str(prop_data);
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SERIALIZE_SPLIT(table_iface)

// instantiate code using _BYNAME macro
BOOST_SERIALIZATION_ASSUME_ABSTRACT(table_iface)
BLUE_SKY_TYPE_SERIALIZE_IMPL_BYNAME(table_iface, "table")

/*-----------------------------------------------------------------
 * serialize prop_iface
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(save, prop_iface)
	// dump all info to string and save it
	ar << (const std::string&)t.to_str();
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(load, prop_iface)
	// dump all info to string and save it
	std::string prop_data;
	ar >> prop_data;
	t.from_str(prop_data);
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SERIALIZE_SPLIT(prop_iface)

// instantiate code using _BYNAME macro
BLUE_SKY_TYPE_SERIALIZE_IMPL_BYNAME(prop_iface, "prop")

/*-----------------------------------------------------------------
 * serialize vartype_table_iface
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN_T(save, vartype_table_iface, 1)
	// save columns number
	t_long num_cols = t.get_n_cols();
	ar << num_cols;
	// save data
	for(t_long i = 0; i < num_cols; ++i) {
		ar << t.get_col_name(i);
		ar << const_cast< type& >(t).get_col_vector(i);
	}
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN_T(load, vartype_table_iface, 1)
	// load columns number
	t_long num_cols;
	ar >> num_cols;
	t.init(num_cols);

	// buffer column vector for adding columns
	typename type::sp_var_type_array_t column = BS_KERNEL.create_object(
		type::var_type_array_t::bs_type()
	);
	// load data
	typename type::vector_t v;
	std::string col_name;
	for(t_long i = 0; i < num_cols; ++i) {
		ar >> col_name;
		ar >> v;

		// column = v
		if(column->size() != v.size())
			column.resize(v.size());
		std::copy(v->begin(), v->end(), column->begin());
		// add column
		t.add_col_vector(i, col_name, column);
	}
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SERIALIZE_SPLIT_T(vartype_table_iface, 1)

BLUE_SKY_TYPE_SERIALIZE_IMPL_T(vartype_table_iface, 1)
// instantiate for t_float
BLUE_SKY_TYPE_SERIALIZE_EXPORT_T(vartype_table_iface, t_float)

