/// @file common_types_serialize.cpp
/// @brief Various common types serialization implementation
/// @author uentity
/// @version 
/// @date 27.01.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "bs_serialize.h"

#include "table.h"
#include "prop.h"
#include "vartype_table.h"
#include "gis.h"
#include "traj.h"

#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

using namespace blue_sky;
namespace boser = boost::serialization;

/*-----------------------------------------------------------------
 * serialize table_iface
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, table)
	// register conversion to interface
	boser::bs_void_cast_register(
		static_cast< table* >(NULL),
		static_cast< table_iface* >(NULL)
	);
	// serialize data
	ar & t.col_names & t.values;
BLUE_SKY_CLASS_SRZ_FCN_END

// instantiate code using _BYNAME macro
BOOST_SERIALIZATION_ASSUME_ABSTRACT(table_iface)
BLUE_SKY_TYPE_SERIALIZE_DECL(table)
BLUE_SKY_TYPE_SERIALIZE_IMPL(table)

/*-----------------------------------------------------------------
 * serialize prop_iface
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN_T(serialize, prop_storage_, 1)
	ar & t.value;
	ar & t.def_value;
	ar & t.flag;
	ar & t.short_name;
	ar & t.description;
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN_T(serialize, prop_impl, 1)
	ar & t.data;
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, prop)
	boser::bs_void_cast_register(
		static_cast< prop* >(NULL),
		static_cast< prop_iface* >(NULL)
	);

	ar & t.fp_impl;
	ar & t.i_impl;
	ar & t.s_impl;
	ar & t.b_impl;
BLUE_SKY_CLASS_SRZ_FCN_END

// instantiate code using _BYNAME macro
BOOST_SERIALIZATION_ASSUME_ABSTRACT(prop_iface)
BLUE_SKY_TYPE_SERIALIZE_DECL(prop)
BLUE_SKY_TYPE_SERIALIZE_IMPL(prop)

/*-----------------------------------------------------------------
 * serialize vartype_table_iface
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN_T(save, vartype_table, 1)
	// save columns number
	const t_long num_cols = t.get_n_cols();
	ar << num_cols;
	// save data
	for(t_long i = 0; i < num_cols; ++i) {
		ar << (const std::wstring&)t.get_col_name(i);
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
	std::wstring col_name;
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

BOOST_SERIALIZATION_ASSUME_ABSTRACT(vartype_table_iface)
// instantiate for t_float
BLUE_SKY_TYPE_SERIALIZE_DECL_T(vartype_table, 1)
BLUE_SKY_TYPE_SERIALIZE_EXPORT_T(vartype_table, t_float)

/*-----------------------------------------------------------------
 * serialize gis_iface
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, gis)
	boser::bs_void_cast_register(
		static_cast< gis* >(NULL),
		static_cast< gis_iface* >(NULL)
	);

	ar & t.sp_table & t.sp_prop;
BLUE_SKY_CLASS_SRZ_FCN_END

// instantiate code using _BYNAME macro
BOOST_SERIALIZATION_ASSUME_ABSTRACT(gis_iface)
BLUE_SKY_TYPE_SERIALIZE_DECL(gis)
BLUE_SKY_TYPE_SERIALIZE_IMPL(gis)

/*-----------------------------------------------------------------
 * serialize traj_iface
 *----------------------------------------------------------------*/
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, traj)
	boser::bs_void_cast_register(
		static_cast< traj* >(NULL),
		static_cast< traj_iface* >(NULL)
	);

	ar & t.sp_table & t.sp_prop;
BLUE_SKY_CLASS_SRZ_FCN_END

// instantiate code using _BYNAME macro
BOOST_SERIALIZATION_ASSUME_ABSTRACT(traj_iface)
BLUE_SKY_TYPE_SERIALIZE_DECL(traj)
BLUE_SKY_TYPE_SERIALIZE_IMPL(traj)

