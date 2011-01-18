/// @file pool_treeish.cpp
/// @brief Treeish implementation of BOS pool with different array
/// @author uentity
/// @date 2009-08-12

#include "bs_bos_core_data_storage_stdafx.h"

#include "pool_treeish.h"
#include <map>
// debug
//#include <iostream>


using namespace blue_sky;
using namespace std;
using namespace blue_sky::pool;

// instantiate used str_val_table specialization
namespace blue_sky {
BS_TYPE_IMPL_T_MEM(str_val_table, bs_pool_node::sp_arrinfo_base);
}

namespace {
// types limiter for pool subnode
struct pool_subn_typelim : public bs_node::restrict_types {
	pool_subn_typelim(const type_descriptor& td)
		: td_(td), descr_(string("This node can only contain ") + td_.stype_)
	{}

	const char* sort_name() const {
		return descr_.c_str();
	}

	types_v accept_types() const {
		types_v t;
		t.push_back(td_);
		return t;
	}

	bool accepts(const sp_link& l) const {
		if(l && upcastable_eq()(l->data()->bs_resolve_type(), td_))
			return true;
		return false;
	}

private:
	type_descriptor td_;
	const string descr_;
};
}	//hidden namespace

/*-----------------------------------------------------------------------------
 *  pool_node_impl definition
 *-----------------------------------------------------------------------------*/
class bs_pool_subnode::pool_subnode_impl {
public:
	// methods
	size_t calc_nlen(const array_size& sz) const {
		return (model_size_->nx * sz.a + sz.bias_a)
			* (model_size_->ny * sz.b + sz.bias_b)
			* (model_size_->nz * sz.c + sz.bias_c);
	}

	void init(const sp_size_triple& model_size) {
		model_size_ = model_size;
	}

	const sp_arrinfo_base find(const string& array_name) const {
		array_info_dict::const_iterator p_ai = arrays_info_->find(array_name);
		if(p_ai != arrays_info_->end())
			return p_ai->second->clone();
		else return NULL;
	}

	size_t get_dimens(const std::string& array_name, size_triple& dimens) const {
		sp_arrinfo_base ai = find(array_name);
		if(ai) {
			const array_size& sz = ai->size;
			dimens = size_triple(
					model_size_->nx * sz.a + sz.bias_a,
					model_size_->ny * sz.b + sz.bias_b,
					model_size_->nz * sz.c + sz.bias_c
					);
			return dimens.nx * dimens.ny * dimens.nz;
		}
		else {
			dimens = size_triple();
			return -1;
		}
	}

	size_t get_nx(const string& array_name) const {
		sp_arrinfo_base ai = find(array_name);
		if(ai)
			return model_size_->nx * ai->size.a + ai->size.bias_a;
		else
			return -1;
	}

	size_t get_ny(const string& array_name) const {
		sp_arrinfo_base ai = find(array_name);
		if(ai)
			return model_size_->ny * ai->size.b + ai->size.bias_b;
		else
			return -1;
	}

	size_t get_nz(const string& array_name) const {
		sp_arrinfo_base ai = find(array_name);
		if(ai)
			return model_size_->nz * ai->size.c + ai->size.bias_c;
		else
			return -1;
	}

	size_t get_nlen(const string& array_name) const {
		sp_arrinfo_base ai = find(array_name);
		if(ai)
			return calc_nlen(ai->size);
		else
			return -1;
	}

	bool add_array(const std::string& array_name, const sp_obj& array, const array_info_base& info) {
		info.resize(array, calc_nlen(info.size));
		bs_node::insert_ret_t res = self_->bs_node::insert(array, array_name);
		//cout << "insert_ret_t = " << int(res.second) << endl;
		if(res.second == bs_node::ins_ok) {
			(*arrays_info_)[array_name] = info.clone();
			return true;
		}
		else return false;
	}

	sp_obj create_array(const std::string& array_name, const array_info_base& info) {
		if(!self_->accepts(info.array_td()))
			return NULL;
		// create array & add it
		sp_obj p_arr = info.create_array(calc_nlen(info.size));
		if(add_array(array_name, p_arr, info))
			return p_arr;
		else
			return NULL;
	}

	sp_obj set_array(const std::string& array_name, const array_info_base& info, const sp_obj& array) {
		// check if array already exists
		sp_obj p_arr = self_->find_array_obj(array_name);
		if(p_arr) {
			(*arrays_info_)[array_name] = info.clone();
			if(array)
				self_->bs_node::erase(array_name);
			else {
				info.resize(p_arr, calc_nlen(info.size));
				return p_arr;
			}
		}
		if(array) {
			add_array(array_name, array, info);
			return array;
		}
		else
			return create_array(array_name, info);
	}

	sp_obj get_non_empty(const std::string& array_name) {
		sp_arrinfo_base info = find(array_name);
		if(!info) {
			// no size for array - we can only throw
			bs_throw_exception("No array_info found for array (" + array_name + ")");
		}

		size_t nlen = calc_nlen(info->size);
		sp_obj p_array = self_->find_array_obj(array_name);

		if(p_array) info->resize(p_array, nlen);
		if(nlen == 0) {
			bs_throw_exception("Array (" + array_name + ") has null size");
		}

		if(!p_array) {
			// create new array with obtained size
			if(!(p_array = create_array(array_name, *info))) {
				bs_throw_exception("Array (" + array_name + ") has null size or cannot be created in this pool");
			}
		}
		return p_array;
	}

	bool rem_array(const std::string& array_name) const {
		if(self_->bs_node::erase(array_name)) {
			arrays_info_->erase(array_name);
			return true;
		}
		return false;
	}

	// get pointer to array's buffer
	void* carray(const std::string& array_name) const {
		sp_obj p_arr = self_->find_array_obj(array_name);
		sp_arrinfo_base ai = find(array_name);
		if(p_arr && ai)
			return ai->carray(p_arr);
		return NULL;
	}

	bool contain(const string& array_name) const {
		return arrays_info_->find(array_name) != arrays_info_->end();
	}

	bool is_empty(const string& array_name) const {
		size_t nlen = get_nlen(array_name);
		return (nlen > 0 && nlen != size_t(-1));
	}

	//ctors
	pool_subnode_impl()
		: model_size_(new size_triple),
		arrays_info_(BS_KERNEL.create_object(array_info_dict::bs_type())), self_(NULL)
	{}
	
	pool_subnode_impl(const type_descriptor& td, const sp_size_triple& model_size)
		: elem_type_(td), model_size_(model_size),
		arrays_info_(BS_KERNEL.create_object(array_info_dict::bs_type())), self_(NULL)
	{
		if(!model_size_)
			model_size_ = new size_triple;
	}

	// copy ctor
	pool_subnode_impl(const pool_subnode_impl& p)
		: elem_type_(p.elem_type_), model_size_(new size_triple),
		arrays_info_(BS_KERNEL.create_object_copy(p.arrays_info_)), self_(NULL)
	{
		if(p.model_size_)
			*model_size_ = *p.model_size_;
	}

	// variables
	typedef str_val_table< sp_arrinfo_base > array_info_dict;
	typedef smart_ptr< array_info_dict > sp_array_info_dict;
	//typedef map< type_descriptor, string > subn_register_t;

	type_descriptor elem_type_;
	sp_size_triple model_size_;
	sp_array_info_dict arrays_info_;
	//subn_register_t subn_register_;
	
	// pointer to parent class
	bs_pool_subnode* self_;
};

/*-----------------------------------------------------------------------------
 *  pool_node_impl definition
 *-----------------------------------------------------------------------------*/
class bs_pool_node::pool_node_impl {
public:
	void init(const size_triple& model_size) {
		*model_size_ = model_size;
	}

	// methods
	size_t calc_nlen(const array_size& sz) const {
		return (model_size_->nx * sz.a + sz.bias_a)
			* (model_size_->ny * sz.b + sz.bias_b)
			* (model_size_->nz * sz.c + sz.bias_c);
	}

	bool register_pool(const type_descriptor& td, const string& subn_name) {
		pair< subn_register_t::iterator, bool > res = subn_register_.insert(
				subn_register_t::value_type(td, subn_name)
				);
		if(res.second) {
			// add subnode
			self_->insert(
					bs_pool_subnode::create(td, model_size_),
						subn_name,
						false
					);
		}
		return res.second;
	}

	std::string find_pool_name(const type_descriptor& td) const {
		// search for subnode name in registry
		subn_register_t::const_iterator p_name = subn_register_.find(td);
		if(p_name != subn_register_.end())
			return p_name->second;
		else 
			return "";
	}

	sp_pool_subnode find_pool(const type_descriptor& td) const {
		// search for corresponding subnode in registry
		subn_register_t::const_iterator p_name = subn_register_.find(td);
		if(p_name == subn_register_.end())
			return NULL;
		//cout << "Subnode name '" << p_name->second << "' found for array type " << p_name->first.stype_ << endl;

		// insert array to subnode found
		bs_node::n_iterator p_subn = self_->bs_node::find(p_name->second);
		if(p_subn == self_->end()) {
			//cout << "add_array: array subnode wasn't found" << endl;
			return NULL;
		}
		
		return sp_pool_subnode(p_subn->node(), bs_dynamic_cast());
	}

	bool add_array(const std::string& array_name, const sp_obj& array, const array_info_base& info) const {
		sp_pool_subnode pp_subn = find_pool(array->bs_resolve_type());
		if(pp_subn)
			return pp_subn->add_array(array_name, array, info);
		else
			return false;
	}

	void push_array(const string& subn_name, const string& array_name, const sp_obj& array, const array_info_base& info) {
		register_pool(array->bs_resolve_type(), subn_name);
		add_array(array_name, array, info);
	}

	bool contain(const type_descriptor& td) const {
		return subn_register_.find(td) != subn_register_.end();
	}

	// ctor
	pool_node_impl()
		: model_size_(new size_triple), subn_register_(), self_(NULL)
	{}

	// copy ctor - make a copy of model_size_;
	pool_node_impl(const pool_node_impl& src)
		: model_size_(new size_triple), subn_register_(src.subn_register_), self_(NULL)
	{
		if(src.model_size_) {
			*model_size_ = *src.model_size_;
		}
	}

	// variables
	typedef map< type_descriptor, string > subn_register_t;

	sp_size_triple model_size_;
	subn_register_t subn_register_;

	// pointer to parent class
	bs_pool_node* self_;
};

/*-----------------------------------------------------------------------------
 *  bs_pool_subnode implementation - forwarding to pimpl
 *-----------------------------------------------------------------------------*/
// ctor
bs_pool_subnode::bs_pool_subnode(bs_type_ctor_param param) 
	: pimpl_(new pool_subnode_impl, mutex(), bs_static_cast())
{
	pimpl_.lock()->self_ = this;
}

//bs_pool_subnode::bs_pool_subnode(bs_type_ctor_param param) 
//	: pimpl_(NULL, mutex(), bs_static_cast())
//{
//	pimpl_.lock()->self_ = this;
//	// parse params
//	smart_ptr< str_data_table > sp_vt(param, bs_dynamic_cast());
//	if(!sp_vt)
//		throw bs_exception("bs_pool_subnode: invalid parameters passed to constructor");
//	else {
//		type_descriptor elem_td = sp_vt->extract_value< s_traits_ptr >("elem_td", type_descriptor());
//		sp_pool_node pool = sp_vt->extract_value< s_traits_ptr >("pool_node", NULL);
//		if(elem_td.is_nil() || pool == NULL) {
//			throw bs_exception("bs_pool_subnode: invalid parameters passed to constructor");
//		}
//		pimpl_ = new node_impl(elem_td, pool->model_size());
//	}
//}

// copy ctor
bs_pool_subnode::bs_pool_subnode(const bs_pool_subnode& p)
	: bs_refcounter(), bs_node(p), pimpl_(new pool_subnode_impl(*p.pimpl_), mutex(), bs_static_cast())
{
	pimpl_.lock()->self_ = this;
}

void bs_pool_subnode::init(const sp_size_triple& model_size) const {
	pimpl_.lock()->model_size_ = model_size;
}

// sweet create function
sp_pool_subnode bs_pool_subnode::create(const type_descriptor& td, const sp_size_triple& model_size) {
	sp_pool_subnode subn(BS_KERNEL.create_object(bs_type()), bs_static_cast());
	if(subn) {
		// restrict holded types
		subn->set_sort(new pool_subn_typelim(td));
		subn->init(model_size);
	}
	return subn;
}

pool::size_t bs_pool_subnode::calc_nlen(const array_size& sz) const {
	return pimpl_->calc_nlen(sz);
}

pool::size_t bs_pool_subnode::get_dimens(const std::string& array_name, size_triple& dimens) const {
	return pimpl_->get_dimens(array_name, dimens);
}

pool::size_t bs_pool_subnode::get_nx(const string& array_name) const {
	return pimpl_->get_nx(array_name);
}

pool::size_t bs_pool_subnode::get_ny(const string& array_name) const {
	return pimpl_->get_ny(array_name);
}

pool::size_t bs_pool_subnode::get_nz(const string& array_name) const {
	return pimpl_->get_nz(array_name);
}

pool::size_t bs_pool_subnode::get_nlen(const string& array_name) const {
	return pimpl_->get_nlen(array_name);
}

bool bs_pool_subnode::add_array(const std::string& array_name, const sp_obj& array, const array_info_base& info) const {
	return pimpl_->add_array(array_name, array, info);
}

sp_obj bs_pool_subnode::create_array(const std::string& array_name, const array_info_base& info) const {
	return pimpl_->create_array(array_name, info);
}

sp_obj bs_pool_subnode::set_array(const std::string& array_name, const array_info_base& info, const sp_obj& array) const {
	return pimpl_->set_array(array_name, info, array);
}

bool bs_pool_subnode::rem_array(const std::string& array_name) const {
	return pimpl_->rem_array(array_name);
}

bool bs_pool_subnode::contain(const std::string& array_name) const {
	return pimpl_->contain(array_name);
}

bool bs_pool_subnode::is_empty(const std::string& array_name) const {
	return pimpl_->is_empty(array_name);
}

bool bs_pool_subnode::accepts(const type_descriptor& td) const {
	bs_node::sort_traits::types_v sup_t = bs_node::get_sort()->accept_types();
	if(std::find_if(sup_t.begin(), sup_t.end(), std::bind2nd(upcastable_eq(), td)) != sup_t.end())
		return true;
	return false;
}

sp_arrinfo_base bs_pool_subnode::find_array(const string& array_name) const {
	return pimpl_->find(array_name);
}

sp_obj bs_pool_subnode::find_array_obj(const string& array_name) const {
	bs_node::n_iterator p_array = bs_node::find(array_name);
	if(p_array != bs_node::end())
		return p_array->data();
	return NULL;
}

sp_obj bs_pool_subnode::get_non_empty(const string& array_name) const {
	return pimpl_->get_non_empty(array_name);
}

void* bs_pool_subnode::carray(const std::string& array_name) const {
	return pimpl_->carray(array_name);
}
/*-----------------------------------------------------------------------------
 *  bs_pool_node implementation - forwards to pimpl
 *-----------------------------------------------------------------------------*/
// ctor
bs_pool_node::bs_pool_node(bs_type_ctor_param) 
	: pimpl_(new pool_node_impl, mutex(), bs_static_cast())
{
	pimpl_.lock()->self_ = this;
	// append d_pool & i_pool subnodes
	pimpl_.lock()->register_pool(bs_array_i::bs_type(), "i_pool");
	pimpl_.lock()->register_pool(bs_array_fp::bs_type(), "d_pool");
}

// copy ctor
bs_pool_node::bs_pool_node(const bs_pool_node& src)
	: bs_refcounter(), bs_node(src), pimpl_(new pool_node_impl(*src.pimpl_), mutex(), bs_static_cast())
{
	pimpl_.lock()->self_ = this;
}

void bs_pool_node::init(const size_triple& model_size) const {
	pimpl_.lock()->init(model_size);
}

pool::size_t bs_pool_node::calc_nlen(const array_size& sz) const {
	return pimpl_->calc_nlen(sz);
}

bool bs_pool_node::register_pool(const type_descriptor& td, const string& subn_name) const {
	return pimpl_.lock()->register_pool(td, subn_name);
}

sp_pool_subnode bs_pool_node::find_pool(const type_descriptor& td) const {
	return pimpl_->find_pool(td);
}

string bs_pool_node::find_pool_name(const type_descriptor& td) const {
	return pimpl_->find_pool_name(td);
}

bool bs_pool_node::contain(const type_descriptor& td) const {
	return pimpl_->contain(td);
}

bool bs_pool_node::add_array(const std::string& array_name, const sp_obj& array, const array_info_base& info) const {
	return pimpl_->add_array(array_name, array, info);
}

sp_obj bs_pool_node::create_array(const std::string& array_name, const array_info_base& info) const {
	sp_pool_subnode p_pool = find_pool(info.array_td());
	if(!p_pool)
		return NULL;
	else
		return p_pool->create_array(array_name, info);
}

void bs_pool_node::push_array(const std::string& subnode_name, const std::string& array_name, 
		const sp_obj& array, const array_info_base& info) const 
{
	return pimpl_.lock()->push_array(subnode_name, array_name, array, info);
}

sp_pool_subnode bs_pool_node::d_pool() const {
	return find_pool(bs_array_fp::bs_type());
}

sp_pool_subnode bs_pool_node::i_pool() const {
	return find_pool(bs_array_i::bs_type());
}

sp_obj bs_pool_node::get_non_empty(const std::string& array_name, const type_descriptor& td) const {
	sp_pool_subnode p_pool = find_pool(td);
	if(!p_pool)
		return NULL;
	else
		return p_pool->get_non_empty(array_name);
}

void* bs_pool_node::carray(const std::string& array_name, const type_descriptor& td) const {
	sp_pool_subnode p_pool = find_pool(td);
	if(!p_pool)
		return NULL;
	else
		return p_pool->carray(array_name);
}


namespace blue_sky {
BLUE_SKY_TYPE_STD_CREATE(bs_pool_subnode);
BLUE_SKY_TYPE_STD_COPY(bs_pool_subnode);

BLUE_SKY_TYPE_STD_CREATE(bs_pool_node);
BLUE_SKY_TYPE_STD_COPY(bs_pool_node);

BLUE_SKY_TYPE_IMPL(bs_pool_subnode, objbase, "bs_pool_subnode", "Array pool subnode containing arrays of single type", "");
BLUE_SKY_TYPE_IMPL(bs_pool_node, objbase, "bs_pool_node", "Array pool represented as storage node", "");

bool register_bs_pool_node(const plugin_descriptor& pd) {
	return BS_KERNEL.register_type(pd, bs_pool_subnode::bs_type());
	return BS_KERNEL.register_type(pd, bs_pool_node::bs_type());
}

}	// namespace blue_sky

