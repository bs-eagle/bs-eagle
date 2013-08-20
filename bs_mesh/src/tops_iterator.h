/// @file tops_iterator.h
/// @brief Simulate iterator over tops array with on-line calculation of needed cell vertices
/// @author uentity
/// @version 1.0
/// @date 02.11.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef TOPS_ITERATOR_Q81ZI8GM
#define TOPS_ITERATOR_Q81ZI8GM

#include "tops_iterator_traits.h"

namespace blue_sky { namespace wpi {

// iterator to random access tops array calculated on the fly
template< template< class, uint > class ti_strategy = carray_ti_traits, uint D = 3 >
class tops_iterator :
	public ti_strategy< std::iterator< std::random_access_iterator_tag, t_float >, D >
{
public:
	typedef ti_strategy< std::iterator< std::random_access_iterator_tag, t_float >, D > strat_t;
	typedef typename strat_t::ctor_param_t strat_ctor_param_t;
	typedef typename strat_t::iterator_t iterator_t;

	typedef typename iterator_t::value_type value_type;
	typedef typename iterator_t::pointer pointer;
	typedef typename iterator_t::reference reference;
	typedef typename iterator_t::difference_type difference_type;

	typedef const value_type& const_reference;
	typedef const value_type* const_pointer;

	typedef t_ulong ulong;
	typedef t_uint uint;

	typedef smart_ptr< rs_smesh_iface, true > sp_smesh;

	enum { n_cell_pts = strat_t::n_cell_pts };

	tops_iterator() : cid_(0), offs_(0) {}

	tops_iterator(strat_ctor_param_t strat_param, const ulong pos = 0)
		: strat_t(strat_param)
	{
		// force cell switching
		cid_ = pos / n_cell_pts + 1;
		switch_pos(pos);
	}
	// std copy ctor is fine

	// common iterator operations
	reference operator*() {
		switch_pos(pos_);
		return strat_t::ss(offs_);
		//return data_[offs_];
	}
	const_reference operator*() const {
		return *const_cast< tops_iterator& >(*this);
	}

	pointer operator->() {
		return &(**this);
	}
	const_pointer operator->() const {
		return &(**this);
	}

	tops_iterator& operator++() {
		++pos_;
		//if(++offs_ == n_cell_pts)
		//	switch_cell(cid_ + 1);
		return *this;
	}
	tops_iterator operator++(int) {
		tops_iterator t = *this;
		++(*this);
		return t;
	}

	tops_iterator& operator--() {
		--pos_;
		//if(offs_ == 0) {
		//	switch_cell(cid_ - 1);
		//	offs_ = n_cell_pts - 1;
		//}
		//else
		//	--offs_;
		return *this;
	}
	tops_iterator operator--(int) {
		tops_iterator t = *this;
		--(*this);
		return t;
	}

	bool operator==(const tops_iterator& rhs) const {
		return pos_ == rhs.pos_;
		//return cid_ == rhs.cid_ && offs_ == rhs.offs_;
	}
	bool operator!=(const tops_iterator& rhs) const {
		return !(*this == rhs);
	}

	void operator=(const tops_iterator& rhs) {
		//mesh_ = rhs.mesh_;
		cid_ = rhs.cid_;
		offs_ = rhs.offs_;
		pos_ = rhs.pos_;
		strat_t::assign(rhs);
		//std::copy(&rhs.data_[0], &rhs.data_[n_cell_pts], &data_[0]);
	}

	// random-access operations
	// Element access operator can destroy iterator position!
	//reference operator[](const difference_type n) {
	//	switch_pos(pos_);
	//	if(ulong(offs_ + n) < n_cell_pts)
	//		return strat_t::ss(offs_ + n);
	//		//return data_[offs_ + n];
	//	else {
	//		switch_pos(pos_ + n);
	//		return strat_t::ss(offs_);
	//		//return data_[offs_];
	//	}
	//}

	value_type operator[](const difference_type n) const {
		return *(*this + n);
	}

	tops_iterator& operator+=(const difference_type n) {
		pos_ += n;
		//if(fit2data(n))
		//	offs_ += n;
		//else
		//	switch_pos(cid_ * n_cell_pts + offs_ + n);
		return *this;
	}
	tops_iterator& operator-=(const difference_type n) {
		pos_ -= n;
		return *this;
		//return this->operator+=(-n);
	}

	tops_iterator operator+(const difference_type n) const {
		tops_iterator t(*this);
		t += n;
		return t;
	}
	tops_iterator operator-(const difference_type n) const {
		tops_iterator t(*this);
		t -= n;
		return t;
	}

	difference_type operator-(const tops_iterator& rhs) const {
		return difference_type(pos_ - rhs.pos_);
		//return (cid_ - rhs.cid_)* n_cell_pts + offs_ - rhs.offs_;
	}

	bool operator<(const tops_iterator& rhs) const {
		return pos_ < rhs.pos_;
		//return cid_ * n_cell_pts + offs_ < rhs.cid_ * n_cell_pts + rhs.offs_;
	}
	bool operator>(const tops_iterator& rhs) const {
		return pos_ > rhs.pos_;
		//return cid_ * n_cell_pts + offs_ > rhs.cid_ * n_cell_pts + rhs.offs_;
	}
	bool operator<=(const tops_iterator& rhs) const {
		return !(*this > rhs);
	}
	bool operator>=(const tops_iterator& rhs) const {
		return !(*this < rhs);
	}

	// if tops_iterator has an underlying backend storage
	// then return index in this storage current iterator points to
	// this can be useful, for ex. for structured grid backend representation
	ulong backend_index() {
		switch_pos(pos_);
		return strat_t::backend_index(offs_);
	}

private:
	// id of currently calculated cell
	ulong cid_, offs_, pos_;

	void switch_cell(const ulong cell_id) {
		strat_t::switch_cell(cell_id);
		//ti_switch_cell(*this, cell_id);
	}

	inline void switch_pos(const ulong pos) {
		// check if new position fits in current data buffer
		offs_ = pos - cid_ * n_cell_pts;
		if(offs_ < n_cell_pts)
			return;

		// if we're here then we need to actually switch to new cell
		cid_ = pos / n_cell_pts;
		switch_cell(cid_);
		offs_ = pos - cid_ * n_cell_pts;
		pos_ = pos;
	}
};

}} // blue_sky::wpi

#endif /* end of include guard: TOPS_ITERATOR_Q81ZI8GM */

