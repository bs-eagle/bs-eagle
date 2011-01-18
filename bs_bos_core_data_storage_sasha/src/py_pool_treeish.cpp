/// @file py_pool_treeish.cpp
/// @brief Python exports of bs_pool_node functionality (implementation)
/// @author uentity
/// @date 2009-08-13

#include "py_pool_treeish.h"
#include "py_bs_tree.h"

#include <boost/python/list.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/to_python_converter.hpp>
//#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/pure_virtual.hpp>

#define SUBNODE static_cast< const bs_pool_subnode* >(sp.get())
#define SUBNODE_L static_cast< bs_pool_subnode* >(sp.lock())

#define POOL static_cast< const bs_pool_node* >(sp.get())
#define POOL_L static_cast< bs_pool_node* >(sp.lock())

using namespace blue_sky;
using namespace blue_sky::pool;
using namespace blue_sky::python;
//using namespace boost::python;
namespace bp = boost::python;

/*-----------------------------------------------------------------------------
 *  py_pool_subnode
 *-----------------------------------------------------------------------------*/
py_pool_subnode::py_pool_subnode(const type_descriptor& td, const sp_size_triple& model_size)
	: py_bs_node(sp_obj(bs_pool_subnode::create(td, model_size)))
{}

py_pool_subnode::py_pool_subnode(const py_pool_subnode& p)
	: py_bs_node(BS_KERNEL.create_object_copy(p.sp))
{}

void py_pool_subnode::init(const bs_pool_node::sp_size_triple& sz) const {
	return SUBNODE->init(sz);
}

pool::size_t py_pool_subnode::get_dimens(const std::string& array_name, size_triple& sz) const {
	return SUBNODE->get_dimens(array_name, sz);
}
pool::size_t py_pool_subnode::get_nx(const std::string& array_name) const {
	return SUBNODE->get_nx(array_name);
}

pool::size_t py_pool_subnode::get_ny(const std::string& array_name) const {
	return SUBNODE->get_ny(array_name);
}

pool::size_t py_pool_subnode::get_nz(const std::string& array_name) const {
	return SUBNODE->get_nz(array_name);
}

pool::size_t py_pool_subnode::get_nlen(const std::string& array_name) const {
	return SUBNODE->get_nlen(array_name);
}

bs_pool_node::sp_arrinfo_base py_pool_subnode::find_array(const std::string& array_name) const {
	return SUBNODE->find_array(array_name);
}

bool py_pool_subnode::add_array(const std::string& array_name, const sp_obj& array, const array_info_base& info) const {
	return SUBNODE->add_array(array_name, array, info);
}

bool py_pool_subnode::rem_array(const std::string& array_name) {
	return SUBNODE->rem_array(array_name);
}

// search functions
bool py_pool_subnode::contain(const std::string& array_name) const {
	return SUBNODE->contain(array_name);
}

bool py_pool_subnode::is_empty(const std::string& array_name) const {
	return SUBNODE->is_empty(array_name);
}

sp_obj py_pool_subnode::create_array(const std::string& array_name, const array_info_base& info) const {
	return SUBNODE->create_array(array_name, info);
}
sp_obj py_pool_subnode::set_array(const std::string& array_name, const array_info_base& info, const sp_obj& array) {
	return SUBNODE->set_array(array_name, info, array);
}

void* py_pool_subnode::carray(const std::string& array_name) const {
	return SUBNODE->carray(array_name);
}

sp_obj py_pool_subnode::get_non_empty(const std::string& array_name) const {
	return SUBNODE->get_non_empty(array_name);
}

/*-----------------------------------------------------------------------------
 *  py_pool_node
 *-----------------------------------------------------------------------------*/
py_pool_node::py_pool_node()
	: py_bs_node(BS_KERNEL.create_object(bs_pool_node::bs_type()))
{}

py_pool_node::py_pool_node(const py_pool_node& p)
	: py_bs_node(BS_KERNEL.create_object_copy(p.sp))
{}

void py_pool_node::init(const size_triple& model_size) const {
	return POOL->init(model_size);
}

bool py_pool_node::register_pool(const type_descriptor& td, const std::string& subnode_name) const {
	return POOL->register_pool(td, subnode_name);
}

sp_pool_subnode py_pool_node::find_pool(const type_descriptor& td) const {
	return POOL->find_pool(td);
}

sp_pool_subnode py_pool_node::d_pool() const {
	return POOL->d_pool();
}

sp_pool_subnode py_pool_node::i_pool() const {
	return POOL->i_pool();
}

std::string py_pool_node::find_pool_name(const type_descriptor& td) const {
	return POOL->find_pool_name(td);
}

bool py_pool_node::contain(const type_descriptor& td) const {
	return POOL->contain(td);
}

bool py_pool_node::add_array(const std::string& array_name, const sp_obj& array, const array_info_base& info) const {
	return POOL->add_array(array_name, array, info);
}

void py_pool_node::push_array(const std::string& subnode_name, const std::string& array_name, 
		const sp_obj& array, const array_info_base& info) const {
	return POOL->push_array(subnode_name, array_name, array, info);
}

sp_obj py_pool_node::create_array(const std::string& array_name, const array_info_base& info) const {
	return POOL->create_array(array_name, info);
}

bool py_pool_node::create_array_fp(const std::string& array_name, const array_info_fp& info)  const {
	return POOL->create_array(array_name, info);
}

bool py_pool_node::create_array_i(const std::string& array_name, const array_info_i& info) const {
	return POOL->create_array(array_name, info);
}

sp_obj py_pool_node::get_non_empty(const std::string& array_name, const type_descriptor& td) const {
	return POOL->get_non_empty(array_name, td);
}

bs_array_fp& py_pool_node::get_non_empty_fp(const std::string& array_name) const {
	return POOL->get_non_empty_fp(array_name);
}

bs_array_i& py_pool_node::get_non_empty_i(const std::string& array_name) const {
	return POOL->get_non_empty_i(array_name);
}

void* py_pool_node::carray(const std::string& array_name, const type_descriptor& td) const {
	return POOL->carray(array_name, td);
}

uchar_t* py_pool_node::carray_i(const std::string& array_name) const {
	return POOL->carray_i(array_name);
}

blue_sky::float_t* py_pool_node::carray_fp(const std::string& array_name) const {
	return POOL->carray_fp(array_name);
}

namespace blue_sky { namespace python {
using namespace blue_sky::pool;
typedef blue_sky::bs_pool_subnode::sp_size_triple sp_size_triple;

namespace {
using namespace bp;

// register conversion from Python sequence to type described by conv_traits
// based on code by Roman Yakovenko
// http://www.language-binding.net/pyplusplus/troubleshooting_guide/automatic_conversion/automatic_conversion.html/
template< class conv_traits >
class from_py_seq_converter {
public:
	static void construct(PyObject* py_obj, converter::rvalue_from_python_stage1_data* data) {
		typedef converter::rvalue_from_python_storage< typename conv_traits::type > storage_t;
		storage_t* the_storage = reinterpret_cast< storage_t* >(data);
		void* memory_chunk = the_storage->storage.bytes;

		// convert PyObject -> python::object
		python::object py_sequence(handle<>(borrowed(py_obj)));

		// create object using placement new
		conv_traits::create_type(memory_chunk, py_sequence);
		data->convertible = memory_chunk;
	}

	static void* convertible(PyObject* py_obj) {
		if( !PySequence_Check( py_obj ) ){
			return 0;
		}

		if( !PyObject_HasAttrString( py_obj, "__len__" ) ){
			return 0;
		}

		python::object py_sequence( handle<>( borrowed( py_obj ) ) );

		if( len( py_sequence ) != conv_traits::ctor_args_num() ){
			return 0;
		}

		if(conv_traits::check_ctor_args(py_sequence)){
			return py_obj;
		}
		else {
			return 0;
		}
	}

	static void register2py() {
		converter::registry::push_back( &convertible
				, &construct
				, type_id< typename conv_traits::type >() );
	}
};

struct size_triple_conv_traits {
	typedef size_triple type;

	static int ctor_args_num() { return 3; }

	static void create_type(void* memory_chunk, const bp::object& py_sequence) {
		size_triple* sz = new(memory_chunk) size_triple(extract< size_t >(py_sequence[0]), extract< size_t >(py_sequence[1]),
				extract< size_t >(py_sequence[2]));
		(void)sz;
	}

	static bool check_ctor_args(const bp::object& py_sequence) {
		bool args_ok = true;
		for(ulong i = 0; i < 3; ++i)
			args_ok &= extract< size_t >(py_sequence[i]).check();
		return args_ok;
	}
};

struct array_size_conv_traits {
	typedef array_size type;

	static int ctor_args_num() { return 6; }

	static void create_type(void* memory_chunk, const bp::object& py_sequence) {
		array_size* sz = new(memory_chunk) array_size(extract< size_t >(py_sequence[0]), extract< size_t >(py_sequence[1]),
				extract< size_t >(py_sequence[2]), extract< size_t >(py_sequence[3]), extract< size_t >(py_sequence[4]), 
				extract< size_t >(py_sequence[5]));
		(void)sz;
	}

	static bool check_ctor_args(const bp::object& py_sequence) {
		bool args_ok = true;
		for(ulong i = 0; i < 6; ++i)
			args_ok &= extract< size_t >(py_sequence[i]).check();
		return args_ok;
	}
};

template< class T >
struct array_info_conv_traits {
	typedef blue_sky::pool::array_info< T > type;

	static int ctor_args_num() { return 2; }

	static void create_type(void* memory_chunk, const bp::object& py_obj) {
		type* info = new(memory_chunk) type(extract< array_size >(py_obj[0]),
				extract< const T& >(py_obj[1]));
		(void)info;
	}

	static bool check_ctor_args(const bp::object& py_obj) {
		bool args_ok = true;
		args_ok &= extract< array_size >(py_obj[0]).check();
		args_ok &= extract< const T& >(py_obj[1]).check();
		return args_ok;
	}
};

struct size_triple2py_tuple {
	static PyObject* convert(const size_triple& sz) {
		bp::tuple t = make_tuple(sz.nx, sz.ny, sz.nz);
		return incref(t.ptr());

		//bp::list values;
		////add all items to "values" list
		//values.append(sz.nx); values.append(sz.ny), values.append(sz.nz);
		////create Python tuple from the list
		//return incref( python::tuple( values ).ptr() );
	}
};

struct array_size2py_tuple {
	static PyObject* convert(const array_size& sz) {
		bp::list values;
		//add all items to "values" list
		values.append(sz.a); values.append(sz.b), values.append(sz.c);
		values.append(sz.bias_a); values.append(sz.bias_b), values.append(sz.bias_c);
		//create Python tuple from the list
		return incref( python::tuple( values ).ptr() );
	}
};

class array_info_base_py : public array_info_base, public bp::wrapper< array_info_base > {
public:
	// ctors
	array_info_base_py() {}
	array_info_base_py(const array_size& sz) : array_info_base(sz) {}
	
	// we need to make copies
	sp_arrinfo_base clone() const {
		return this->get_override("clone")();
	}

	size_t length(const sp_obj& array) const {
		return this->get_override("length")(array);
	}

	bool resize(const sp_obj& array, size_t new_size) const {
		return this->get_override("resize")(array, new_size);
	}

	bool reset(const sp_obj& array) const {
		return this->get_override("reset")(array);
	}

	void* carray(const sp_obj& array) const {
		return this->get_override("carray")(array);
	}

	sp_obj create_array(size_t nlen) const {
		return this->get_override("create_array")(nlen);
	}

	type_descriptor array_td() const {
		return this->get_override("array_td")();
	}
};

}	// hidden namespace

void py_export_bs_pool_node() {
	using namespace bp;
	// export size_triple
	bp::class_< size_triple >("size_triple")
		.def(bp::init< size_t, size_t, size_t >())
		.def_readwrite("nx", &size_triple::nx)
		.def_readwrite("ny", &size_triple::ny)
		.def_readwrite("nz", &size_triple::nz)
		;
	// register conversion size_triple -> Python tuple
	//bp::to_python_converter< size_triple, size_triple2py_tuple >();

	// export array_size
	bp::class_< array_size >("array_size")
		.def(bp::init< size_t, size_t, size_t, size_t, size_t, size_t >())
		.def_readwrite("a", &array_size::a)
		.def_readwrite("b", &array_size::b)
		.def_readwrite("c", &array_size::c)
		.def_readwrite("bias_a", &array_size::bias_a)
		.def_readwrite("bias_b", &array_size::bias_b)
		.def_readwrite("bias_c", &array_size::bias_c)
		;
	// register conversion size_triple -> Python tuple
	//bp::to_python_converter< bs_pool_node::array_size, array_size2py_tuple >();

	// export array_info
	bp::class_< array_info_base_py, boost::noncopyable >("array_info_base")
		.def(bp::init< const array_size& >())
		.def_readwrite("size", &array_info_base_py::size)
		.def("clone", pure_virtual(&array_info_base_py::clone))
		.def("length", pure_virtual(&array_info_base_py::length))
		.def("resize", pure_virtual(&array_info_base_py::resize))
		.def("reset", pure_virtual(&array_info_base_py::reset))
		//.def("carray", pure_virtual(&array_info_base_py::carray))
		.def("create_array", pure_virtual(&array_info_base_py::create_array))
		.def("array_td", pure_virtual(&array_info_base_py::array_td))
		;

	bp::class_< array_info< float_t >, bases< array_info_base > >("array_info_fp")
		.def(bp::init< const array_size&, const float_t& >())
		.def_readwrite("def_value", &array_info< float_t >::def_value)
		.def("clone", &array_info< float_t >::clone)
		.def("length", &array_info< float_t >::length)
		.def("resize", &array_info< float_t >::resize)
		.def("reset", &array_info< float_t >::reset)
		//.def("carray", &array_info< float_t >::carray)
		.def("create_array", &array_info< float_t >::create_array)
		.def("array_td", &array_info< float_t >::array_td)
		;

	bp::class_< array_info< uchar_t >, bases< array_info_base > >("array_info_i")
		.def(bp::init< const array_size&, const uchar_t& >())
		.def_readwrite("def_value", &array_info< uchar_t >::def_value)
		.def("clone", &array_info< uchar_t >::clone)
		.def("length", &array_info< uchar_t >::length)
		.def("resize", &array_info< uchar_t >::resize)
		.def("reset", &array_info< uchar_t >::reset)
		//.def("carray", &array_info< uchar_t >::carray)
		.def("create_array", &array_info< uchar_t >::create_array)
		.def("array_td", &array_info< uchar_t >::array_td)
		;

	// register from Python converters
	//from_py_converters::register2py();
	from_py_seq_converter< size_triple_conv_traits >::register2py();
	from_py_seq_converter< array_size_conv_traits >::register2py();
	from_py_seq_converter< array_info_conv_traits< float_t > >::register2py();
	from_py_seq_converter< array_info_conv_traits< uchar_t > >::register2py();

	// export pool_subnode class
	bp::class_< py_pool_subnode, bp::bases< py_bs_node > >("pool_subnode",
		init< const type_descriptor&, const sp_size_triple& >())
		.def("init", &py_pool_subnode::init)
		.def("get_dimens", &py_pool_subnode::get_dimens)
		.def("get_nx", &py_pool_subnode::get_nx)
		.def("get_ny", &py_pool_subnode::get_ny)
		.def("get_nz", &py_pool_subnode::get_nz)
		.def("get_nlen", &py_pool_subnode::get_nlen)
		.def("find_array", &py_pool_subnode::find_array, bp::return_value_policy< bp::return_by_value >())
		.def("add_array", &py_pool_subnode::add_array)
		.def("rem_array", &py_pool_subnode::rem_array)
		.def("create_array", &py_pool_subnode::create_array)
		.def("set_array", &py_pool_subnode::set_array)
		.def("contain", &py_pool_subnode::contain)
		.def("is_empty", &py_pool_subnode::is_empty)
		//.def("carray", &py_pool_subnode::carray)
		.def("get_non_empty", &py_pool_subnode::get_non_empty)
		;

	// export pool_node class
	bp::class_< py_pool_node, bp::bases< py_bs_node > >("pool_node")
		.def("init", &py_pool_node::init)
		.def("register_pool", &py_pool_node::register_pool)
		.def("find_pool", &py_pool_node::find_pool)
		.def("d_pool", &py_pool_node::d_pool)
		.def("i_pool", &py_pool_node::i_pool)
		.def("find_pool_name", &py_pool_node::find_pool_name)
		.def("contain", &py_pool_node::contain)
		.def("add_array", &py_pool_node::add_array)
		.def("push_array", &py_pool_node::push_array)
		.def("create_array", &py_pool_node::create_array)
		.def("create_array", &py_pool_node::create_array_i)
		.def("create_array", &py_pool_node::create_array_fp)
		.def("get_non_empty", &py_pool_node::get_non_empty)
		.def("get_non_empty", &py_pool_node::get_non_empty_fp, bp::return_value_policy< bp::return_by_value >())
		.def("get_non_empty", &py_pool_node::get_non_empty_i, bp::return_value_policy< bp::return_by_value >())
		//.def("carray", &py_pool_node::carray)
		//.def("carray", &py_pool_node::carray_fp)
		//.def("carray", &py_pool_node::carray_i)
		;
}

} }	// namespace blue_sky { python {

