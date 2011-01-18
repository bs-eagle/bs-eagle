/// @file py_pool_treeish.h
/// @brief Python export of bs_pool_node functionality
/// @author uentity
/// @date 2009-08-13

#ifndef PY_POOL_TREEISH_H
#define PY_POOL_TREEISH_H

#include "pool_treeish.h"
#include "py_bs_tree.h"

namespace blue_sky {
namespace python {

class py_pool_subnode : public py_bs_node {
public:
	typedef pool::size_t size_t;
	typedef pool::size_triple size_triple;
	typedef pool::sp_size_triple sp_size_triple;
	typedef pool::array_size array_size;
	typedef pool::array_info_base array_info_base;
	typedef array_info_base::sp_arrinfo_base sp_arrinfo_base;
	typedef pool::array_info_fp array_info_fp;
	typedef pool::array_info_i array_info_i;
	
	// ctors
	py_pool_subnode(const type_descriptor& td, const sp_size_triple& model_size);
	py_pool_subnode(const py_pool_subnode&);

	void init(const bs_pool_node::sp_size_triple& model_size) const;
	bool accepts(const type_descriptor& td) const;

	size_t get_dimens(const std::string& array_name, size_triple& sz) const;
	size_t get_nx(const std::string& array_name) const;
	size_t get_ny(const std::string& array_name) const;
	size_t get_nz(const std::string& array_name) const;

	size_t get_nlen(const std::string& array_name) const;

	// search functions
	bool contain(const std::string& array_name) const;
	bool is_empty(const std::string& array_name) const;	
	sp_arrinfo_base find_array(const std::string& array_name) const;
	sp_obj find_array_obj(const std::string& array_name) const;

	bool add_array(const std::string& array_name, const sp_obj& array, const array_info_base& info) const;
	//bool add_array_i(const std::string& array_name, sp_i_array array, const array_size& sz) const;
	bool rem_array(const std::string& array_name);

	sp_obj create_array(const std::string& array_name, const array_info_base& info) const;
	// set_array can be used to add, create or resize array with given name
	sp_obj set_array(const std::string& array_name, const array_info_base& info, const sp_obj& array = NULL);

	// get pointer to array's buffer
	void* carray(const std::string& array_name) const;

	// check that array has size > 0
	// throws in case of null model_size
	sp_obj get_non_empty(const std::string& array_name) const;
};

class py_pool_node : public py_bs_node {
public:
	typedef bs_pool_subnode::sp_pool_subnode sp_pool_subnode;

	typedef pool::size_t size_t;
	typedef pool::size_triple size_triple;
	typedef pool::sp_size_triple sp_size_triple;
	typedef pool::array_size array_size;
	typedef pool::array_info_base array_info_base;
	typedef array_info_base::sp_arrinfo_base sp_arrinfo_base;
	typedef pool::array_info_fp array_info_fp;
	typedef pool::array_info_i array_info_i;

	// ctors
	py_pool_node();
	py_pool_node(const py_pool_node&);
	
	void init(const size_triple& model_size) const;

	bool register_pool(const type_descriptor& td, const std::string& subnode_name) const;
	std::string find_pool_name(const type_descriptor& td) const;
	sp_pool_subnode find_pool(const type_descriptor& td) const;

	sp_pool_subnode d_pool() const;
	sp_pool_subnode i_pool() const;

	bool contain(const type_descriptor& td) const;

	bool add_array(const std::string& array_name, const sp_obj& array, const array_info_base& info) const;

	void push_array(const std::string& subnode_name, const std::string& array_name,
			const sp_obj& array, const array_info_base& info) const;

	sp_obj create_array(const std::string& array_name, const array_info_base& info) const;

	bool create_array_fp(const std::string& array_name, const pool::array_info< float_t >& info) const;
	bool create_array_i(const std::string& array_name, const pool::array_info< uchar_t >& info) const;

	sp_obj get_non_empty(const std::string& array_name, const type_descriptor& td) const;
	// get_non_empty syntax sugar for int & real arrays
	bs_array_fp& get_non_empty_fp(const std::string& array_name) const;
	bs_array_i& get_non_empty_i(const std::string& array_name) const;

	// get pointer to array's buffer
	void* carray(const std::string& array_name, const type_descriptor& td) const;
	uchar_t* carray_i(const std::string& array_name) const;
	float_t* carray_fp(const std::string& array_name) const;
};

void py_export_bs_pool_node();

}	// namespace python
}	// namespace blue_sky

#endif	// file guard

