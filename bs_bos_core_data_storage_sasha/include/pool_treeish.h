/// @file pool_treeish.h
/// @brief Treeish implementation of BOS pool with different arrays
///
/// Pool is represented as data storage node (bs_node).
/// @author uentity
/// @date 2009-08-11


#ifndef BS_ARRAY_POOL_TREEISH_H
#define BS_ARRAY_POOL_TREEISH_H

#include "bs_kernel.h"
#include "bs_tree.h"
#include "bs_prop_base.h"
#include "throw_exception.h"

#include "throw_exception.h"

//#define D_ELEM_TYPE float
//#define I_ELEM_TYPE unsigned char

namespace blue_sky {
// tag int and float array types
typedef unsigned char             uchar_t;
typedef float                     float_t;

typedef bs_array< float_t >  bs_array_fp;
typedef bs_array< uchar_t >  bs_array_i;

typedef smart_ptr< bs_array_fp > sp_array_fp;
typedef smart_ptr< bs_array_i > sp_array_i;

namespace pool {

//typedef ulong size_t;

enum dimens_idx  //! indexes for dimension parameters
{
	ARRAY_POOL_NX_A,
	ARRAY_POOL_NX_B,
	ARRAY_POOL_NY_A,
	ARRAY_POOL_NY_B,
	ARRAY_POOL_NZ_A,
	ARRAY_POOL_NZ_B,
	ARRAY_POOL_TOTAL
};

template <typename size_t>
struct size_triple {
	size_t nx, ny, nz;

	size_triple() : nx(0), ny(0), nz(0) {}

	size_triple(size_t nx_, size_t ny_, size_t nz_)
		: nx(nx_), ny(ny_), nz(nz_)
	{}

	size_triple(const size_t* pbuf)
		: nx(pbuf[0]), ny(pbuf[1]), nz(pbuf[2])
	{}
};

template <typename size_t>
struct array_size {
	//T def_value;
	size_t a, b, c;
	size_t bias_a, bias_b, bias_c;

	array_size()
		: a(0), b(0), c(0), bias_a(0), bias_b(0), bias_c(0)
	{}

	array_size(size_t a_, size_t b_, size_t c_, size_t bias_a_, size_t bias_b_, size_t bias_c_)
		: a(a_), b(b_), c(c_), bias_a(bias_a_), bias_b(bias_b_), bias_c(bias_c_)
	{}

	template< class int_t >
	array_size(const int_t* pbuf)
		: a(pbuf[ARRAY_POOL_NX_A]), b(pbuf[ARRAY_POOL_NY_A]), c(pbuf[ARRAY_POOL_NZ_A]),
		bias_a(pbuf[ARRAY_POOL_NX_B]), bias_b(pbuf[ARRAY_POOL_NY_B]), bias_c(pbuf[ARRAY_POOL_NZ_B])
	{}

	//bool operator==(const array_size& lhs) {
	//	return (a == lhs.a
	//			&& b == lhs.b
	//			&& c == lhs.c
	//			&& bias_a == lhs.bias_a
	//			&& bias_b == lhs.bias_b
	//			&& bias_c == lhs.bias_c);
	//}
};

/*
class array_info_base {
public:
	typedef smart_ptr< array_info_base, false > sp_arrinfo_base;

	array_size size;

	// ctors
	array_info_base() : size() {}
	array_info_base(const array_size& sz) : size(sz) {}
	// we need to make copies
	virtual sp_arrinfo_base clone() const = 0;

	virtual size_t length(const sp_obj& array) const = 0;
	virtual bool resize(const sp_obj& array, size_t new_size) const = 0;
	// fill array with default values
	virtual bool reset(const sp_obj& array) const = 0;
	// obtain raw pointer to array's buffer
	virtual void* carray(const sp_obj& array) const = 0;
	// fab method for creating array's instances
	virtual sp_obj create_array(size_t nlen) const = 0;
	// obvious
	virtual type_descriptor array_td() const = 0;

	// comparison
	//virtual bool equal(const sp_arrinfo_base& lhs) const {
	//	return size == lhs->size;
	//}

	//bool operator==(const sp_arrinfo_base& lhs) const {
	//	return equal(lhs);
	//}
	//bool operator!=(const sp_arrinfo_base& lhs) const {
	//	return !equal(lhs);
	//}
};
typedef array_info_base::sp_arrinfo_base sp_arrinfo_base;

template< class T >
struct array_info : public array_info_base {
	typedef bs_array< T > array_t;
	typedef smart_ptr< array_t > sp_array;
	typedef lsmart_ptr< sp_array > lsp_array;

	T def_value;

	// ctors
	array_info() : def_value() {}

	array_info(const array_size& sz) 
		: array_info_base(sz), def_value()
	{}

	array_info(const array_size& sz, const T& default_value)
		: array_info_base(sz), def_value(default_value)
	{}

	// override virtual base functions
	sp_arrinfo_base clone() const {
		return new array_info< T >(*this);
	}

	size_t length(const sp_obj& array) const {
		sp_array p_array(array, bs_dynamic_cast());
		if(p_array)
			return static_cast< size_t >(p_array->size());
		else 
			return size_t(-1);
	}

	bool resize(const sp_obj& array, size_t new_size) const {
		sp_array p_array(array, bs_dynamic_cast());
		if(p_array) {
			p_array->resize(new_size);
			return true;
		}
		else return false;
	}

	bool reset(const sp_obj& array) const {
		lsp_array p_array(sp_array(array, bs_dynamic_cast()));
		if(p_array) {
			std::fill(p_array->begin(), p_array->end(), def_value);
			return true;
		}
		else return false;
	}

	void* carray(const sp_obj& array) const {
		// very unsafe - ignores locking strategy
		sp_array p_array(array, bs_dynamic_cast());
		if(p_array)
			return static_cast< void* >(&(*p_array)[0]);
		else
			return NULL;
	}

	sp_obj create_array(size_t nlen) const {
		lsp_array p_arr(sp_array(BS_KERNEL.create_object(array_t::bs_type()), bs_static_cast()));
		p_arr->resize(nlen);
		std::fill(p_arr->begin(), p_arr->end(), def_value);
		return p_arr;
	}

	type_descriptor array_td() const {
		return array_t::bs_type();
	}
};

typedef array_info< uchar_t > array_info_i;
typedef array_info< float_t > array_info_fp;

}
*/

/// @brief Class that represents subnode containing arrays of the same type
template <typename index_t, typename item_t>
class BS_API_PLUGIN bs_pool_subnode : public bs_node {
public:
	typedef smart_ptr< bs_pool_subnode, true > sp_pool_subnode;
  
	typedef index_t size_t;
	typedef pool::size_triple<size_t>     size_triple;
	typedef smart_ptr<size_triple, true>  sp_size_triple;
	typedef pool::array_size array_size;

	size_t get_dimens(const std::string& array_name, size_triple& sz) const;
	size_t get_nx(const std::string& array_name) const;
	size_t get_ny(const std::string& array_name) const;
	size_t get_nz(const std::string& array_name) const;

	size_t get_nlen(const std::string& array_name) const;

	// check if given type can be stored in this pool
	bool accepts(const type_descriptor& td) const;

	bool add_array(const std::string& array_name, const sp_obj& array, const array_info_base& info) const;
	bool rem_array(const std::string& array_name) const;
	sp_obj create_array(const std::string& array_name, const array_info_base& info) const;
	// set_array can be used to add, create or resize array with given name
	sp_obj set_array(const std::string& array_name, const array_info_base& info, const sp_obj& array = NULL) const;

	// search functions
	sp_arrinfo_base find_array(const std::string& array_name) const;
	sp_obj find_array_obj(const std::string& array_name) const;

	bool contain(const std::string& array_name) const;
	bool is_empty(const std::string& array_name) const;

	// get pointer to array's buffer
	void* carray(const std::string& array_name) const;

	// check that array has size > 0
	// throws in case of null model_size
	sp_obj get_non_empty(const std::string& array_name) const;

	// change model size
	void init(const sp_size_triple& model_size) const;
	// syntax sugar for creating subnode
	static sp_pool_subnode create(const type_descriptor& td, const sp_size_triple& model_size);

private:
	class pool_subnode_impl;
	smart_ptr< pool_subnode_impl, false > pimpl_;

	size_t calc_nlen(const array_size& sz) const;

	BLUE_SKY_TYPE_DECL(bs_pool_subnode);
};
typedef bs_pool_subnode::sp_pool_subnode sp_pool_subnode;

/// @brief Class representing arrays pool as node with bs_pool_subnodes
class BS_API_PLUGIN bs_pool_node : public bs_node {
	friend class bs_pool_subnode;
public:
	typedef pool::size_t size_t;
	typedef pool::size_triple size_triple;
	typedef pool::sp_size_triple sp_size_triple;
	typedef pool::array_size array_size;
	typedef pool::array_info_base array_info_base;
	typedef array_info_base::sp_arrinfo_base sp_arrinfo_base;
	typedef pool::array_info_fp array_info_fp;
	typedef pool::array_info_i array_info_i;
	//using pool::array_info;

	void init(const size_triple& model_size) const;

	bool register_pool(const type_descriptor& td, const std::string& subnode_name) const;

	sp_pool_subnode find_pool(const type_descriptor& td) const;
	std::string find_pool_name(const type_descriptor& td) const;

	template< class T >
	sp_pool_subnode find_pool() const {
		return find_subnode(bs_array< T >::bs_type());
	}

	// syntax sugar for accessing always existing int & real pools
	sp_pool_subnode d_pool() const;
	sp_pool_subnode i_pool() const;

	bool contain(const type_descriptor& td) const;

	bool add_array(const std::string& array_name, const sp_obj& array, const array_info_base& info) const;

	void push_array(const std::string& subnode_name, const std::string& array_name,
			const sp_obj& array, const array_info_base& info) const;

	sp_obj create_array(const std::string& array_name, const array_info_base& info) const;
	// added for compatibility with old bs_array_map syntax
	template< class T >
	smart_ptr< bs_array< T > > create_array(const std::string& array_name, const array_size& sz, const T& def_value) {
		typedef bs_array< T > array_t;
		typedef smart_ptr< array_t > sp_array;
		return
			sp_array(create_array(array_name, pool::array_info< T >(sz, def_value)), bs_static_cast());
	}

	// create syntax sugar for int & real arrays
	sp_array_i create_array_i(const std::string& array_name, const array_info_i& info) const {
		return sp_array_i(create_array(array_name, info), bs_static_cast());
	}

	sp_array_fp create_array_fp(const std::string& array_name, const array_info_fp& info) const {
		return sp_array_fp(create_array(array_name, info), bs_static_cast());
	}

	sp_obj get_non_empty(const std::string& array_name, const type_descriptor& td) const;

	// WARNING: ignores locking
	template< class T >
	bs_array< T >& get_non_empty(const std::string& array_name) const {
		typedef bs_array< T > array_t;
		typedef smart_ptr< array_t > sp_array;
		return *const_cast< array_t* >(sp_array(get_non_empty(array_name, array_t::bs_type()), bs_static_cast()).get());
	}
	// get_non_empty syntax sugar for int & real arrays
	bs_array_fp& get_non_empty_fp(const std::string& array_name) const {
		return get_non_empty< float_t >(array_name);
	}

	bs_array_i& get_non_empty_i(const std::string& array_name) const {
		return get_non_empty< uchar_t >(array_name);
	}

	// get pointer to array's buffer
	// WARNING: ignores locking
	void* carray(const std::string& array_name, const type_descriptor& td) const;

	template< class T >
	T* carray(const std::string& array_name) const {
		typedef bs_array< T > array_t;
		typedef smart_ptr< array_t > sp_array;
		return (T*)carray(array_name, array_t::bs_type());
	}

	uchar_t* carray_i(const std::string& array_name) const {
		return carray< uchar_t >(array_name);
	}
	float_t* carray_fp(const std::string& array_name) const {
		return carray< float_t >(array_name);
	}

private:
	class pool_node_impl;
	smart_ptr< pool_node_impl, false > pimpl_;

	sp_size_triple model_size() const;
	size_t calc_nlen(const array_size& sz) const;

	BLUE_SKY_TYPE_DECL(bs_pool_node);
};

bool register_bs_pool_node(const plugin_descriptor& pd);

} // blue_sky namespace

#endif // file guard

