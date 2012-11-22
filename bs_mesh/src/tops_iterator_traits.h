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

	carray_ti_buf_traits(ctor_param_t*)
	{}
	//carray_ti_buf_traits(ctor_param_t* = NULL) : just_born(true) {}

	reference ss(ulong offs) {
		return data_[offs];
	}

	bool switch_buf(ulong cell_id) const {
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

protected:
	rs_smesh_iface* mesh_;
};

// shortcuts for use in clients
template< class iterator_type >
struct carray_ti_traits : public mesh_ti_sc_traits< carray_ti_buf_traits< iterator_type > > {
	typedef carray_ti_buf_traits< iterator_type > buf_traits_t;
	typedef mesh_ti_sc_traits< buf_traits_t > sc_traits_t;
	typedef typename sc_traits_t::ctor_param_t ctor_param_t;

	carray_ti_traits(ctor_param_t p)
		: sc_traits_t(p)
	{}
};

template< class iterator_type >
struct bufpool_ti_traits : public mesh_ti_sc_traits< bufpool_ti_buf_traits< iterator_type > > {
	typedef bufpool_ti_buf_traits< iterator_type > buf_traits_t;
	typedef mesh_ti_sc_traits< buf_traits_t > sc_traits_t;
	typedef typename sc_traits_t::ctor_param_t ctor_param_t;

	bufpool_ti_traits(ctor_param_t p)
		: sc_traits_t(p)
	{}
};

/*-----------------------------------------------------------------
 * Strategy that calculate cell vertices using structured grid representation
 *----------------------------------------------------------------*/
template< class buf_traits >
struct sgrid_ti_sc_traits : public buf_traits {
	typedef buf_traits buf_traits_t;
	typedef typename buf_traits::ctor_param_t buf_ctor_param_t;
	typedef typename buf_traits::pointer pointer;
	typedef typename v_float::const_iterator cvf_iterator;

	// traits should be initialized with iterator to beginning of sgrid
	// + dimeshions of original mesh
	struct sgrid_handle {
		cvf_iterator sgrid;
		ulong nx, ny, size;
	};

	typedef typename mesh_ti_sc_traits< buf_traits >::template
		deduce_ctor_param< sgrid_handle, buf_ctor_param_t* > ctor_param_deducer;
	// actually deduce ctor parameter
	typedef typename ctor_param_deducer::type ctor_param_t;

	sgrid_ti_sc_traits(ctor_param_t p)
		: buf_traits(ctor_param_deducer::buf_param(p)),
		  start_(ctor_param_deducer::sc_param(p).sgrid),
		  nx_(ctor_param_deducer::sc_param(p).nx),
		  ny_(ctor_param_deducer::sc_param(p).ny),
		  sz_(ctor_param_deducer::sc_param(p).size)
	{}

	template< class inp_iterator, class outp_iterator >
	static void copy_points(const inp_iterator& src, outp_iterator& dst, ulong n = 1) {
		dst = std::copy(src, src + n*3, dst);
	}

	void switch_cell(const ulong cell_id) {
		// are we inside mesh bounds?
		if(cell_id >= sz_)
			return;

		if(buf_traits::switch_buf(cell_id)) return;

		// calc offset of cell beginning in structured grid
		const ulong plane_sz = nx_ * ny_;
		const ulong z = cell_id / plane_sz;
		const ulong y = (cell_id - plane_sz * z) / nx_;
		const ulong x = cell_id - plane_sz * z - nx_ * y;

		const cvf_iterator start = start_ +
			(z * (nx_ + 1) * (ny_ + 1) + y * (nx_ + 1) + x) * 3;

		// fill local tops_iterator buffer
		pointer pdata = &buf_traits::ss(0);
		// vert 0, 1;
		pdata = std::copy(start, start + 6, pdata);
		//*pdata++ = *start; *pdata++ = *(start + 1);
		// vert 2, 3
		pdata = std::copy(start + (nx_ + 1)*3, start + (nx_ + 3)*3, pdata);
		//*pdata++ = *(start + (nx_ + 1)); *pdata++ = *(start + (nx_ + 2));
		// vert 4, 5
		pdata = std::copy(start + (nx_ + 1)*(ny_ + 1)*3, start + ((nx_ + 1)*(ny_ + 1) + 2)*3, pdata);
		//*pdata++ = *(start + (nx_ + 1)*(ny_ + 1));
		// vert 5
		//*pdata++ = *(start + ((nx_ + 1)*(ny_ + 1) + 1));
		// vert 6
		std::copy(start + (nx_ + 1)*(ny_ + 2)*3, start + ((nx_ + 1)*(ny_ + 2) + 2)*3, pdata);
		//*pdata++ = *(start + (nx_ + 1)*(ny_ + 2));
		// vert 7
		//*pdata++ = *(start + ((nx_ + 1)*(ny_ + 2) + 1));
	}

protected:
	// pointer to the beginning of structured grid array
	cvf_iterator start_;
	ulong nx_, ny_, sz_;
};

// shortcuts for use in clients
template< class iterator_type >
struct carray_sgrid_ti_traits : public sgrid_ti_sc_traits< carray_ti_buf_traits< iterator_type > > {
	typedef carray_ti_buf_traits< iterator_type > buf_traits_t;
	typedef sgrid_ti_sc_traits< buf_traits_t > sc_traits_t;
	typedef typename sc_traits_t::ctor_param_t ctor_param_t;

	carray_sgrid_ti_traits(ctor_param_t p)
		: sc_traits_t(p)
	{}
};

template< class iterator_type >
struct bufpool_sgrid_ti_traits : public sgrid_ti_sc_traits< bufpool_ti_buf_traits< iterator_type > > {
	typedef bufpool_ti_buf_traits< iterator_type > buf_traits_t;
	typedef sgrid_ti_sc_traits< buf_traits_t > sc_traits_t;
	typedef typename sc_traits_t::ctor_param_t ctor_param_t;

	bufpool_sgrid_ti_traits(ctor_param_t p)
		: sc_traits_t(p)
	{}
};

}} /* blue_sky::wpi */

#endif /* end of include guard: TOPS_ITERATOR_TRAITS_M6AE7JT3 */

