/// @file common_types_serialize.cpp
/// @brief Various common types serialization implementation
/// @author uentity
/// @version 
/// @date 27.01.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "common_types_serialize.h"
#include "table.h"
#include "prop.h"

#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

using namespace blue_sky;
namespace boser = boost::serialization;

//BLUE_SKY_CLASS_SRZ_FCN_DECL(serialize, blue_sky::table)
//BLUE_SKY_CLASS_SRZ_FCN_DECL(serialize, blue_sky::prop)
//BLUE_SKY_CLASS_SRZ_FCN_DECL_T(serialize, blue_sky::vartype_table, 1)

/*-----------------------------------------------------------------
 * serialize table_iface
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(save, table)
	// dump all info to string and save it
	ar << (const std::string&)t.to_str();
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(load, table)
	// dump all info to string and save it
	std::string prop_data;
	ar >> prop_data;
	t.from_str(prop_data);
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, table)
	boser::bs_void_cast_register(
		static_cast< table* >(NULL),
		static_cast< table_iface* >(NULL)
	);
	boser::split_free(ar, t, version);
BLUE_SKY_CLASS_SRZ_FCN_END

//BLUE_SKY_CLASS_SERIALIZE_SPLIT(table_iface)

// instantiate code using _BYNAME macro
BOOST_SERIALIZATION_ASSUME_ABSTRACT(table_iface)
BLUE_SKY_TYPE_SERIALIZE_IMPL(table)
//BLUE_SKY_TYPE_SERIALIZE_IMPL_BYNAME(table_iface, "table")

/*-----------------------------------------------------------------
 * serialize prop_iface
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(save, prop)
	// dump all info to string and save it
	ar << (const std::string&)t.to_str();
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(load, prop)
	// dump all info to string and save it
	std::string prop_data;
	ar >> prop_data;
	t.from_str(prop_data);
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, prop)
	boser::bs_void_cast_register(
		static_cast< prop* >(NULL),
		static_cast< prop_iface* >(NULL)
	);
	boser::split_free(ar, t, version);
BLUE_SKY_CLASS_SRZ_FCN_END

//BLUE_SKY_CLASS_SERIALIZE_SPLIT(prop_iface)

// instantiate code using _BYNAME macro
BOOST_SERIALIZATION_ASSUME_ABSTRACT(prop_iface)
BLUE_SKY_TYPE_SERIALIZE_IMPL(prop)
//BLUE_SKY_TYPE_SERIALIZE_IMPL_BYNAME(prop, "prop")

/*-----------------------------------------------------------------
 * serialize vartype_table_iface
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN_T(save, vartype_table, 1)
	// save columns number
	t_long num_cols = t.get_n_cols();
	ar << num_cols;
	// save data
	for(t_long i = 0; i < num_cols; ++i) {
		ar << (const std::string&)t.get_col_name(i);
		ar << const_cast< type& >(t).get_col_vector(i);
	}
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN_T(load, vartype_table, 1)
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
			column->resize(v.size());
		if(column->size())
			std::copy(v.begin(), v.end(), column->begin());
		// add column
		t.add_col_vector(i, col_name, column);
	}
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN_T(serialize, vartype_table, 1)
	typedef typename type::vector_t::value_type var_type_t;
	boser::bs_void_cast_register(
		static_cast< type* >(NULL),
		static_cast< vartype_table_iface< var_type_t >* >(NULL)
	);
	boser::split_free(ar, t, version);
BLUE_SKY_CLASS_SRZ_FCN_END
//BLUE_SKY_CLASS_SERIALIZE_SPLIT_T(vartype_table_iface, 1)

BOOST_SERIALIZATION_ASSUME_ABSTRACT(vartype_table_iface)
//BLUE_SKY_TYPE_SERIALIZE_IMPL_T(vartype_table, 1)
// instantiate for t_float
BLUE_SKY_TYPE_SERIALIZE_IMPL_T(vartype_table, t_float)

