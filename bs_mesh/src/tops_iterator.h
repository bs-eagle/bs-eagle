/// @file tops_iterator.h
/// @brief Simulate iterator over tops array with on-line calculation of needed cell vertices
/// @author uentity
/// @version 1.0
/// @date 02.11.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef TOPS_ITERATOR_Q81ZI8GM
#define TOPS_ITERATOR_Q81ZI8GM

#include "mesh_grdecl.h"

// iterator to random access tops array calculated on the fly
class tops_iterator : public std::iterator< std::random_access_iterator_tag, t_float > {
public:
	typedef std::iterator< std::random_access_iterator_tag, t_float > base_t;
	typedef typename base_t::value_type value_type;
	typedef typename base_t::pointer pointer;
	typedef typename base_t::reference reference;
	typedef typename base_t::difference_type difference_type;

	typedef const value_type& const_reference;
	typedef const value_type* const_pointer;

	typedef typename mesh_grdecl::element_t element_t;
	typedef typename element_t::fpoint3d_t fpoint3d_t;

    typedef t_ulong ulong;
    typedef t_uint uint;

	enum { n_cell_pts = 24 };

	tops_iterator() : mesh_(NULL), cid_(0), offs_(0) {}
	tops_iterator(const mesh_grdecl& mesh, const ulong pos = 0)
		: mesh_(&mesh)
	{
		switch_pos(pos);
	}
	// std copy ctor is fine

	// common iterator operations
	reference operator*() {
		return data_[offs_];
	}
	const_reference operator*() const {
		return data_[offs_];
	}

	pointer operator->() {
		return &data_[offs_];
	}
	const_pointer operator->() const {
		return &data_[offs_];
	}

	tops_iterator& operator++() {
		if(++offs_ == n_cell_pts)
			switch_cell(cid_ + 1);
		return *this;
	}
	tops_iterator operator++(int) {
		tops_iterator t = *this;
		++(*this);
		return t;
	}

	tops_iterator& operator--() {
		if(offs_ == 0) {
			switch_cell(cid_ - 1);
			offs_ = n_cell_pts - 1;
		}
		else
			--offs_;
		return *this;
	}
	tops_iterator operator--(int) {
		tops_iterator t = *this;
		--(*this);
		return t;
	}

	bool operator==(const tops_iterator& rhs) const {
		return cid_ == rhs.cid_ && offs_ == rhs.offs_;
	}
	bool operator!=(const tops_iterator& rhs) const {
		return !(*this == rhs);
	}

	void operator=(const tops_iterator& rhs) {
		mesh_ = rhs.mesh_;
		cid_ = rhs.cid_;
		offs_ = rhs.offs_;
		std::copy(&rhs.data_[0], &rhs.data_[n_cell_pts], &data_[0]);
	}

	// random-access operations
	// Element access operator destroys iterator position!
	reference operator[](const ulong n) {
		if(offs_ + n < n_cell_pts)
			return data_[offs_ + n];
		else {
			switch_pos(cid_ * n_cell_pts + offs_ + n);
			return data_[offs_];
		}
	}

	tops_iterator& operator+=(const ulong n) {
		switch_pos(cid_ * n_cell_pts + offs_ + n);
		return *this;
	}
	tops_iterator& operator-=(const ulong n) {
		switch_pos(cid_ * n_cell_pts + offs_ - n);
		return *this;
	}

	tops_iterator operator+(const ulong n) {
		tops_iterator t(*this);
		t += n;
		return t;
	}
	tops_iterator operator-(const ulong n) {
		tops_iterator t(*this);
		t -= n;
		return t;
	}

	ulong operator-(const tops_iterator& rhs) const {
		return (cid_ - rhs.cid_)* n_cell_pts + offs_ - rhs.offs_;
	}

	bool operator<(const tops_iterator& rhs) const {
		return cid_ * n_cell_pts + offs_ < rhs.cid_ * n_cell_pts + rhs.offs_;
	}
	bool operator>(const tops_iterator& rhs) const {
		return cid_ * n_cell_pts + offs_ > rhs.cid_ * n_cell_pts + rhs.offs_;
	}
	bool operator<=(const tops_iterator& rhs) const {
		return !(*this > rhs);
	}
	bool operator>=(const tops_iterator& rhs) const {
		return !(*this < rhs);
	}

private:
	// ref to mesh object
	const mesh_grdecl* mesh_;
	// id of currently calculated cell
	ulong cid_, offs_;
	// cache of calculated cell corner coord
	value_type data_[n_cell_pts];

	void switch_cell(const ulong cell_id) {
		const ulong plane_sz = mesh_->nx * mesh_->ny;
		const ulong z = cell_id / plane_sz;
		const ulong y = (cell_id - z * plane_sz) / mesh_->nx;
		// cid_ == -1 marks end() interator
		// offset is reset to 0
		if(z >= ulong(mesh_->nz))
			cid_ = ulong(-1);
		else
			cid_ = cell_id;
		offs_ = 0;

		fpoint3d_t corner;
		element_t e;
		if(cell_id < ulong(mesh_->n_elements)) {
			mesh_->calc_element(cell_id - mesh_->nx * (z * mesh_->ny + y), y, z, e);
			for(uint c = 0; c < 8; ++c) {
				corner = e.get_corners()[c];
				data_[c * 3] = corner.x;
				data_[c * 3 + 1] = corner.y;
				data_[c * 3 + 2] = corner.z;
			}
		}
	}

	void switch_pos(const ulong pos) {
		const ulong cell_id = pos / n_cell_pts;
		switch_cell(cell_id);
		offs_ = pos - cell_id * n_cell_pts;
	}

	void advance(const ulong delta) {
	}
};

#endif /* end of include guard: TOPS_ITERATOR_Q81ZI8GM */

