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

/*-----------------------------------------------------------------
 * Strategies that describe tops_iterator local buffer managing policy
 *----------------------------------------------------------------*/

template< class iterator_type >
struct carray_ti_buf_traits : public iterator_type {
	typedef iterator_type iterator_t;
	typedef typename iterator_t::value_type value_type;
	typedef typename iterator_t::reference reference;

	enum { n_cell_pts = 24 };
	typedef void ctor_param_t;

	carray_ti_buf_traits(ctor_param_t* = NULL) {}
	//carray_ti_buf_traits(ctor_param_t* = NULL) : just_born(true) {}

	reference ss(const ulong offs) {
		return data_[offs];
	}

	bool switch_buf(const ulong) const {
		return false;
	}

	void assign(const carray_ti_buf_traits& rhs) {
		std::copy(&rhs.data_[0], &rhs.data_[n_cell_pts], &data_[0]);
	}

	value_type data_[n_cell_pts];
	// flag to help bufpool determine if buffer was just created
	//bool just_born;
};

// bufpool strategy collects pool of unique cells buffers
// and even more reduces memory consuption
template< class iterator_type >
struct bufpool_ti_buf_traits : public iterator_type {
	typedef iterator_type iterator_t;
	typedef typename iterator_t::value_type value_type;
	typedef typename iterator_t::reference reference;
	typedef typename iterator_t::pointer pointer;

	// use simple traits as buffer holder
	typedef carray_ti_buf_traits< iterator_t > cell_buffer;
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

	bufpool_ti_buf_traits(ctor_param_t* pstore = NULL) 
		: store_(pstore)
	{}

	reference ss(const ulong offs) {
		return *(data_ + offs);
	}

	// return true, if given cell is FOUND in cache
	// true prevents recalc
	bool switch_buf(const ulong cell_id) {
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

	void assign(const bufpool_ti_buf_traits& rhs) {
		store_ = rhs.store_;
		data_ = rhs.data_;
	}

	cell_buf_storage* store_;
	pointer data_;
};

/*-----------------------------------------------------------------
 * switch_cell calculating cell vertices through mesh object
 *----------------------------------------------------------------*/
template< class buf_traits >
struct mesh_ti_sc_traits : public buf_traits {
	typedef buf_traits buf_traits_t;
	typedef typename buf_traits::ctor_param_t buf_ctor_param_t;
	typedef typename buf_traits::pointer pointer;

	// if buf_ctor_param_t is void, we can omit it
	template< class Sc_param, class Buf_param, class = void >
	struct deduce_ctor_param {
		typedef std::pair< Sc_param, Buf_param > type;

		static Sc_param& sc_param(type& t) {
			return t.first;
		}
		static Buf_param& buf_param(type& t) {
			return t.second;
		}
	};
	template< class Sc_param, class unused >
	struct deduce_ctor_param< Sc_param, void*, unused > {
		typedef Sc_param type;

		static Sc_param& sc_param(type& t) {
			return t;
		}
		static void* buf_param(type& t) {
			return NULL;
		}
	};
	typedef deduce_ctor_param< rs_smesh_iface*, buf_ctor_param_t* > ctor_param_deducer;
	// actually deduce traits ctor param
	typedef typename ctor_param_deducer::type ctor_param_t;

	// empty ctor
	mesh_ti_sc_traits() : mesh_(NULL) {}

	// std stor
	mesh_ti_sc_traits(ctor_param_t p)
		: buf_traits(ctor_param_deducer::buf_param(p)),
		  mesh_(ctor_param_deducer::sc_param(p))
	{}

	void switch_cell(const ulong cell_id) {
		typedef grd_ecl::fpoint3d fpoint3d_t;

		// obtain mesh dimensions
		rs_smesh_iface::index_point3d_t dims = mesh_->get_dimens();
		const ulong sz = dims[0] * dims[1] * dims[2];

		// are we inside mesh bounds?
		if(cell_id >= sz)
			return;

		if(!buf_traits::switch_buf(cell_id)) {
			const ulong plane_sz = dims[0] * dims[1];
			const ulong z = cell_id / plane_sz;
			const ulong y = (cell_id - z * plane_sz) / dims[0];

			const grd_ecl::fpoint3d_vector corners = mesh_->calc_element(
				cell_id - z * plane_sz - y * dims[0], y, z
			);
			// copy corners to plain data array
			pointer pdata = &buf_traits::ss(0);
			for(uint c = 0; c < 8; ++c) {
				const fpoint3d_t& corner = corners[c];
				*pdata++ = corner.x;
				*pdata++ = corner.y;
				*pdata++ = corner.z;
			}
		}
	}

	void assign(const mesh_ti_sc_traits& rhs) {
		buf_traits::assign(rhs);
		mesh_ = rhs.mesh_;
	}

	ulong backend_index(const ulong offs) const {
		// meaningless operation, beacause no backend buffer exists
		// cells are retrieved 'online' from mesh
		return offs;
	}

protected:
	rs_smesh_iface* mesh_;
};

// shortcuts for use in clients
template< class iterator_type >
struct carray_ti_traits : public mesh_ti_sc_traits< carray_ti_buf_traits< iterator_type > > {
	typedef carray_ti_buf_traits< iterator_type > buf_traits_t;
	typedef mesh_ti_sc_traits< buf_traits_t > sc_traits_t;
	typedef typename sc_traits_t::ctor_param_t ctor_param_t;

	// empty ctor
	carray_ti_traits() {}

	carray_ti_traits(ctor_param_t p)
		: sc_traits_t(p)
	{}
};

template< class iterator_type >
struct bufpool_ti_traits : public mesh_ti_sc_traits< bufpool_ti_buf_traits< iterator_type > > {
	typedef bufpool_ti_buf_traits< iterator_type > buf_traits_t;
	typedef mesh_ti_sc_traits< buf_traits_t > sc_traits_t;
	typedef typename sc_traits_t::ctor_param_t ctor_param_t;

	// empty ctor
	bufpool_ti_traits() {}

	bufpool_ti_traits(ctor_param_t p)
		: sc_traits_t(p)
	{}
};

/*-----------------------------------------------------------------
 * Strategy that calculate cell vertices using structured grid representation
 *----------------------------------------------------------------*/
// we SHOULD store final buffer, because later in wpi::cell direct
// reinterpret_cast is used to access cell's vertex coords
// so storing only offset indices in sgrid is not enough and lead to ERRORS
// thats why we inherit from carray_ti_buf_traits
template< class iterator_type >
struct sgrid_ti_traits : public carray_ti_buf_traits< iterator_type > {
	typedef carray_ti_buf_traits< iterator_type > base_t;
	typedef iterator_type iterator_t;
	typedef typename iterator_t::value_type value_type;
	typedef typename iterator_t::reference reference;
	typedef typename v_float::iterator vf_iterator;

	enum { n_cell_pts = base_t::n_cell_pts };
	typedef ulong cellv_index[n_cell_pts];

	// traits should be initialized with iterator to beginning of sgrid
	// + dimeshions of original mesh
	struct sgrid_handle {
		vf_iterator sgrid;
		ulong nx, ny, size;
	};
	typedef const sgrid_handle& ctor_param_t;

	using base_t::data_;

	// empty ctor
	sgrid_ti_traits() : h_(dumb_sh()), cell_offs_(0) {}

	// std ctor
	sgrid_ti_traits(const sgrid_handle& h)
		: h_(h), cell_offs_(0)
	{
		const ulong line =   (h_.nx + 1) * 3;
		const ulong plane =  (h_.nx + 1) * (h_.ny + 1) * 3;
		const ulong lplane = (h_.nx + 1) * (h_.ny + 2) * 3;
		const cellv_index t = {
			// vert 0, 1
			0, 1, 2, 3, 4, 5,
			// vert 2, 3
			line, line + 1, line + 2, line + 3, line + 4, line + 5,
			// vert 4, 5
			plane, plane + 1, plane + 2, plane + 3, plane + 4, plane + 5,
			// vert 6, 7
			lplane, lplane + 1, lplane + 2, lplane + 3, lplane + 4, lplane + 5
		};
		std::copy(&t[0], &t[n_cell_pts], &vidx_[0]);
	}

	void switch_cell(const ulong cell_id) {
		// are we inside mesh bounds?
		if(cell_id >= h_.size)
			return;

		// calc offset of cell beginning in structured grid
		const ulong plane_sz = h_.nx * h_.ny;
		const ulong z = cell_id / plane_sz;
		const ulong y = (cell_id - plane_sz * z) / h_.nx;
		const ulong x = cell_id - plane_sz * z - h_.nx * y;

		cell_offs_ = (z * (h_.nx + 1) * (h_.ny + 1) + y * (h_.nx + 1) + x) * 3;
		// copy cell data from sgrid to local buffer
		for(ulong i = 0; i < n_cell_pts; ++i)
			data_[i] = *(h_.sgrid + cell_offs_ + vidx_[i]);
	}

	// subscript comes from base buffer class
	//reference ss(const ulong offs) {
	//	return *(start_ + cell_offs_ + vidx_[offs]);
	//}

	// assign only iterators belonging to the same mesh!
	// or, more precisely, to mesh with identical dimensions
	void assign(const sgrid_ti_traits& rhs) {
		//nx_ = rhs.nx_; ny_ = rhs.ny_; sz_ = rhs.sz_;
		//start_ = rhs.start_;
		cell_offs_ = rhs.cell_offs_;
		std::copy(&rhs.vidx_[0], &rhs.vidx_[n_cell_pts], &vidx_[0]);
		// assign data buffer
		base_t::assign(rhs);
	}

	ulong backend_index(const ulong offs) const {
		// return index of current iterator inside structured grid buffer
		return cell_offs_ + vidx_[offs];
	}

	// don't copy handle in every tops_iterator -- store refernce to it
	// because sgrid_handle for all iterators lives in trimesh
	// and all iterators are invalid when trimesh is destroyed
	const sgrid_handle& h_;
	ulong cell_offs_;
	cellv_index vidx_;

private:
	// just for empty ctor
	static const sgrid_handle& dumb_sh() {
		static const sgrid_handle h = {
			vf_iterator(), 0, 0, 0
		};
		return h;
	}
};

}} /* blue_sky::wpi */

#endif /* end of include guard: TOPS_ITERATOR_TRAITS_M6AE7JT3 */

