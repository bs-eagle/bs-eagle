/// @file tops_iterator_traits.h
/// @brief Triats for tops_iterator
/// @author uentity
/// @version 1.0
/// @date 21.11.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef TOPS_ITERATOR_TRAITS_M6AE7JT3
#define TOPS_ITERATOR_TRAITS_M6AE7JT3

#include "conf.h"
#include "rs_smesh_iface.h"
#include <boost/pool/pool_alloc.hpp>

namespace blue_sky { namespace wpi {

template< class iterator_type >
struct carray_ti_traits : public iterator_type {
	typedef iterator_type iterator_t;
	typedef typename iterator_t::value_type value_type;
	typedef typename iterator_t::reference reference;
	//typedef const value_type& const_reference;

	enum { n_cell_pts = 24 };
	typedef void ctor_param_t;

	carray_ti_traits(ctor_param_t* = NULL) {}
	//carray_ti_traits(ctor_param_t* = NULL) : just_born(true) {}

	reference ss(ulong offs) {
		return data_[offs];
	}
	//const_reference ss(ulong idx) const {
	//	return data_[idx];
	//}

	bool switch_buf(ulong cell_id) const {
		return false;
	}

	void assign(const carray_ti_traits& rhs) {
		std::copy(&rhs.data_[0], &rhs.data_[n_cell_pts], &data_[0]);
	}

	void operator=(const carray_ti_traits& rhs) {
		assign(rhs);
	}

	value_type data_[n_cell_pts];
	// flag to help bufpool determine if buffer was just created
	//bool just_born;
};

// bufpool strategy collects pool of unique cells buffers
// and even more reduces memory consuption
template< class iterator_type >
struct bufpool_ti_traits : public iterator_type {
	typedef iterator_type iterator_t;
	typedef typename iterator_t::value_type value_type;
	typedef typename iterator_t::reference reference;
	typedef typename iterator_t::pointer pointer;
	//typedef const value_type& const_reference;

	// use simple traits as buffer holder
	typedef carray_ti_traits< iterator_t > cell_buffer;
	enum { n_cell_pts = cell_buffer::n_cell_pts };

	// helper structure that actually contain cell data buffer
	//typedef value_type cell_buffer[n_cell_pts];

	typedef boost::fast_pool_allocator<
		std::pair< ulong, cell_buffer >,
		boost::default_user_allocator_new_delete,
		boost::details::pool::null_mutex
		> cell_buf_allocator;
	typedef std::map< ulong, cell_buffer, std::less< ulong >, cell_buf_allocator > cell_buf_storage;
	typedef typename cell_buf_storage::value_type cbs_value_type;
	typedef typename cell_buf_storage::iterator cbs_iterator;
	typedef std::pair< cbs_iterator, bool > cbs_ins_res;

	typedef cell_buf_storage ctor_param_t;

	bufpool_ti_traits(ctor_param_t* pstore = NULL) : store_(pstore) {}

	reference ss(ulong offs) {
		return *(data_ + offs);
	}

	// return true, if given cell is FOUND in cache
	// true prevents recalc
	bool switch_buf(ulong cell_id) {
		static cell_buffer t;
		if(!store_) return false;

		//cell_buffer& b = (*store_)[cell_id];
		//data_ = &b.data_[0];
		//const bool res = !b.just_born;
		//b.just_born = false;
		//return res;
		////if(b.just_born) {
		////	b.just_born = false;
		////	return false;
		////}
		////return true;

		// check if given cell is already cached
		//cbs_iterator p_cb = store_->find(cell_id);
		//if(p_cb != store_->end()) {
		//	data_ = &p_cb->second.data_[0];
		//	return true;
		//}
		//else {
		//	data_ = &(*store_)[cell_id].data_[0];
		//	return false;
		//}

		cbs_ins_res r = store_->insert(cbs_value_type(cell_id, t));
		data_ = &r.first->second.data_[0];
		return !r.second;
	}

	void assign(const bufpool_ti_traits& rhs) {
		store_ = rhs.store_;
		data_ = rhs.data_;
	}

	cell_buf_storage* store_;
	pointer data_;
};

}} /* blue_sky::wpi */

#endif /* end of include guard: TOPS_ITERATOR_TRAITS_M6AE7JT3 */

