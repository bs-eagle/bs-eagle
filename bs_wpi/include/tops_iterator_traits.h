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
template< class iterator_type, uint cell_buf_sz_ >
struct carray_ti_buf_traits : public iterator_type {
	typedef iterator_type iterator_t;
	typedef typename iterator_t::value_type value_type;
	typedef typename iterator_t::reference reference;

	enum { cell_buf_sz = cell_buf_sz_ };
	typedef void ctor_param_t;

	carray_ti_buf_traits(ctor_param_t* = NULL) {
		std::fill(&data_[0], &data_[cell_buf_sz], 0);
	}
	//carray_ti_buf_traits(ctor_param_t* = NULL) : just_born(true) {}

	reference ss(const ulong offs) {
		return data_[offs];
	}

	bool switch_buf(const ulong) const {
		return false;
	}

	void assign(const carray_ti_buf_traits& rhs) {
		std::copy(&rhs.data_[0], &rhs.data_[cell_buf_sz], &data_[0]);
	}

	value_type data_[cell_buf_sz];
	// flag to help bufpool determine if buffer was just created
	//bool just_born;
};

// bufpool strategy collects pool of unique cells buffers
// and even more reduces memory consuption
template< class iterator_type, uint cell_buf_sz_ >
struct bufpool_ti_buf_traits : public iterator_type {
	typedef iterator_type iterator_t;
	typedef typename iterator_t::value_type value_type;
	typedef typename iterator_t::reference reference;
	typedef typename iterator_t::pointer pointer;

	// use simple traits as buffer holder
	typedef carray_ti_buf_traits< iterator_t, cell_buf_sz_ > cell_buffer;
	enum { cell_buf_sz = cell_buffer::cell_buf_sz };

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
template< class buf_traits, uint D >
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
			for(uint c = 0; c < (1 << D); ++c) {
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
template< class iterator_type, uint D >
struct carray_ti_traits :
	public mesh_ti_sc_traits< carray_ti_buf_traits< iterator_type, (1 << D) * 3 >, D >
{
	typedef carray_ti_buf_traits< iterator_type, (1 << D) * 3 > buf_traits_t;
	typedef mesh_ti_sc_traits< buf_traits_t, D > sc_traits_t;
	typedef typename sc_traits_t::ctor_param_t ctor_param_t;

	// n_cell_pts should come somewhere, for ex from strategy
	// n_cell_pts = 2^D * N_coord_per_point, defines buffer size fo one cell
	// for current 2D & 3D strategies assume that N_coord_per_point = 3, so
	// n_cell_pts = 24 for 3D, = 12 for 2D
	// but we can have for ex. 2 coords per point in 2D, in that case n_cell_pts = 8
	//
	// update: we SHOULD have n_cell_pts = 24, because client code relies on fact that
	// tops_iterator steps over 3D tops array, where each cell occupy 24 elements
	// this is used in index calculations, etc
	enum { n_cell_pts = 24 };

	// empty ctor
	carray_ti_traits() {}

	carray_ti_traits(ctor_param_t p)
		: sc_traits_t(p)
	{}
};

template< class iterator_type, uint D >
struct bufpool_ti_traits :
	public mesh_ti_sc_traits< bufpool_ti_buf_traits< iterator_type, (1 << D) * 3 >, D >
{
	typedef bufpool_ti_buf_traits< iterator_type, (1 << D) * 3 > buf_traits_t;
	typedef mesh_ti_sc_traits< buf_traits_t, D > sc_traits_t;
	typedef typename sc_traits_t::ctor_param_t ctor_param_t;

	enum { n_cell_pts = 24 };

	// empty ctor
	bufpool_ti_traits() {}

	bufpool_ti_traits(ctor_param_t p)
		: sc_traits_t(p)
	{}
};

/*-----------------------------------------------------------------
 * Handles describing underlying structured grid for corresponding traits
 *----------------------------------------------------------------*/
// standard handle - describes structured grid represented as continuous array
// handle should be initialized with iterator to beginning of sgrid
// + dimenshions of original mesh
template< uint D >
struct sgrid_array_handle {
	enum { cell_buf_sz = (1 << D) * 3 };
	typedef ulong cellv_index[cell_buf_sz];
	typedef typename v_float::iterator vf_iterator;

private :
	// helper for vidx initialization
	template< uint D_, class = void >
	struct init_vidx {
		typedef ulong vertex_pos_i[D_];

		init_vidx(const ulong nx, const ulong ny, cellv_index& res) {
			const ulong line =   (nx + 1) * 3;
			const ulong plane =  (nx + 1) * (ny + 1) * 3;
			const ulong lplane = (nx + 1) * (ny + 2) * 3;
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
			std::copy(&t[0], &t[cell_buf_sz], &res[0]);
		}

		// convert mesh cell id to offset of cell's coord in srgid buffer
		static ulong cell_id2offs(const ulong nx, const ulong ny, const ulong cell_id) {
			// calc offset of cell beginning in structured grid
			const ulong plane_sz = nx * ny;
			const ulong z = cell_id / plane_sz;
			const ulong y = (cell_id - plane_sz * z) / nx;
			const ulong x = cell_id - plane_sz * z - nx * y;

			return (z * (nx + 1) * (ny + 1) + y * (nx + 1) + x) * 3;
		}

		// convert point offset in sgrid to x, y, z offsets
		static void offs2pos(const ulong nx, const ulong ny, const ulong offs, vertex_pos_i& res) {
			const ulong plane_sz = (nx + 1) * (ny + 1);
			res[2] = offs / plane_sz;
			res[1] = (offs - plane_sz * res[2]) / (nx + 1);
			res[0] = offs - plane_sz * res[2] - (nx + 1) * res[1];
		}
	};

	// specialization for 2D
	template< class unused >
	struct init_vidx< 2, unused > {
		typedef ulong vertex_pos_i[2];

		init_vidx(const ulong nx, const ulong ny, cellv_index& res) {
			const ulong line =   (nx + 1) * 3;
			const cellv_index t = {
				// vert 0, 1
				0, 1, 2, 3, 4, 5,
				// vert 2, 3
				line, line + 1, line + 2, line + 3, line + 4, line + 5
			};
			std::copy(&t[0], &t[cell_buf_sz], &res[0]);
		}

		// convert mesh cell id to offset of cell's coord in srgid buffer
		static ulong cell_id2offs(const ulong nx, const ulong ny, const ulong cell_id) {
			// calc offset of cell beginning in structured grid
			const ulong y = cell_id / nx;
			return (y * (nx + 1) + cell_id - nx * y) * 3;
		}

		// convert point offset in sgrid to x, y offsets
		static void offs2pos(const ulong nx, const ulong /* ny */, const ulong offs, vertex_pos_i& res) {
			res[1] = offs / (nx + 1);
			res[0] = offs - (nx + 1)*res[1];
		}
	};

public :
	typedef init_vidx< D > init_vidx_t;

	vf_iterator sgrid;
	ulong nx, ny, size;
	cellv_index vidx;

	// empty ctor
	sgrid_array_handle()
		: sgrid(vf_iterator()), nx(0), ny(0), size(0)
	{}

	// normal ctor
	sgrid_array_handle(const vf_iterator& sgrid_, const ulong nx_, const ulong ny_, const ulong sz)
		: sgrid(sgrid_), nx(nx_), ny(ny_), size(sz) {
		// calc v_idx
		init_vidx_t(nx, ny, vidx);
	}

	void init(const vf_iterator& sgrid_, const ulong nx_, const ulong ny_, const ulong sz) {
		sgrid = sgrid_;
		nx = nx_; ny = ny_; size = sz;
		init_vidx_t(nx, ny, vidx);
	}

	// convert mesh cell id to offset of cell's coord in srgid buffer
	ulong cell_id2offs(const ulong cell_id) const {
		return init_vidx_t::cell_id2offs(nx, ny, cell_id);
	}
}; // eof sgrid_handle

// sgrid abstraction handle representing rectangular grid
// it actually differs from sgrid_array_handle only in type of sgrid member
// which is not a simple iterator now
template< uint D >
struct rgrid_handle {
	typedef sgrid_array_handle< D > sgrid_handle;
	typedef typename sgrid_handle::cellv_index cellv_index;
	typedef typename sgrid_handle::init_vidx_t init_vidx_t;
	typedef typename init_vidx_t::vertex_pos_i vertex_pos_i;

	// modified version of dim_subscript from coord_zcorn_tools
	// doesn't store accumulated sum, thus allow random access
	struct dim_subscript {
		dim_subscript(const spv_float& dim = NULL, const t_float offset = 0) {
			init(dim, offset);
		}

		void init(const spv_float& dim, const t_float offset) {
			dim_ = dim; offset_ = offset;
			if(!dim_)
				ss_fcn_ = &dim_subscript::ss_null_dim;
			else if(dim_->size() == 1)
				ss_fcn_ = &dim_subscript::ss_const_dim;
			else
				ss_fcn_ = &dim_subscript::ss_array_dim;
		}

		// choke for empty rgrid ctor
		t_float ss_null_dim(const ulong) const {
			return 0.0;
		}

		t_float ss_const_dim(const ulong idx) const {
			return offset_ + dim_->ss(0) * idx;
		}

		t_float ss_array_dim(const ulong idx) const {
			t_float res = offset_;
			for(ulong i = 0; i < idx; ++i)
				res += dim_->ss(i);
			return res;
		}

		t_float operator[](const ulong idx) const {
			return (this->*ss_fcn_)(idx);
		}

	private:
		spv_float dim_;
		t_float (dim_subscript::*ss_fcn_)(const ulong) const;
		t_float offset_;
	};

	// replace iterator over sgrid with this struct
	// support minimal operators needed
	struct r2s_adapter {
		r2s_adapter(const rgrid_handle& h, const ulong offs = 0)
			: h_(h), offs_(offs)
		{}

		// returns copy of r2s_adapter
		r2s_adapter operator+(const t_long offs) const {
			return r2s_adapter(h_, offs_ + offs);
		}

		r2s_adapter& operator++() {
			offs_++;
			return *this;
		}

		t_float operator*() const {
			// calc point offset
			// assume that we have N_coords_per_point = 3
			const ulong pt_offs = offs_ / 3;
			const ulong c_offs = offs_ - pt_offs * 3;
			// convert offset to x, y, z coords
			vertex_pos_i pos;
			init_vidx_t::offs2pos(h_.dim[0], h_.dim[1], pt_offs, pos);

			if(c_offs < D)
				return h_.d_ss[c_offs][pos[c_offs]];
			else
				return 0;
		}

		const rgrid_handle& h_;
		// offset is raw, not point offset!
		ulong offs_;
	};

	// fake sgrid
	r2s_adapter sgrid;
	// dimensions
	ulong dim[D];
	ulong size;
	// deltas subscript along each dimension
	dim_subscript d_ss[D];
	// like sgrid_array_handle
	cellv_index vidx;

	// empty ctor
	rgrid_handle() : sgrid(*this) {
		std::fill(&dim[0], &dim[D], 0);
		size = 0;
	}

	// normal ctor via init
	rgrid_handle(const ulong(& dim_)[D], const spv_float(& deltas_)[D],
		const t_float(& min_pt_)[D]) : sgrid(*this)
	{
		init();
	}

	void init(const ulong(& dim_)[D], const spv_float(& deltas_)[D],
		const t_float(& min_pt_)[D])
	{
		size = 1;
		for(uint i = 0; i < D; ++i) {
			dim[i] = dim_[i];
			size *= dim[i];
			d_ss[i].init(deltas_[i], min_pt_[i]);
		}

		// init vidx
		init_vidx_t(dim[0], dim[1], vidx);
	}

	// convert mesh cell id to offset of cell's coord in srgid buffer
	ulong cell_id2offs(const ulong cell_id) const {
		return init_vidx_t::cell_id2offs(dim[0], dim[1], cell_id);
	}
};

/*-----------------------------------------------------------------
 * Strategy that calculate cell vertices using structured grid representation
 *----------------------------------------------------------------*/
// we SHOULD store final buffer, because later in wpi::cell direct
// reinterpret_cast is used to access cell's vertex coords
// so storing only offset indices in sgrid is not enough and lead to ERRORS
// thats why we inherit from carray_ti_buf_traits
//
// NOTE: sgrid trait actually depends on dimensionality
// so we need to pass here dimensions num
// strictly speaking we also need N_coord_per_point, but assume it = 3
template< class iterator_type, uint D, template< uint > class sgrid_handle_t = sgrid_array_handle >
struct sgrid_ti_traits_impl : public carray_ti_buf_traits< iterator_type, (1 << D) * 3 > {
	typedef carray_ti_buf_traits< iterator_type, (1 << D) * 3 > base_t;
	typedef iterator_type iterator_t;
	typedef typename iterator_t::value_type value_type;
	typedef typename iterator_t::reference reference;
	typedef typename v_float::iterator vf_iterator;

	enum { n_cell_pts = 24 };
	enum { cell_buf_sz = base_t::cell_buf_sz };

	typedef sgrid_handle_t< D > sgrid_handle;
	typedef const sgrid_handle& ctor_param_t;
	using base_t::data_;

	// empty ctor
	sgrid_ti_traits_impl() : h_(NULL), cell_offs_(-1) {}

	// std ctor
	sgrid_ti_traits_impl(const sgrid_handle& h)
		: h_(&h), cell_offs_(h_->size)
	{}

	void switch_cell(const ulong cell_id) {
		// are we inside mesh bounds?
		if(cell_id >= h_->size)
			return;

		const ulong new_offs = h_->cell_id2offs(cell_id);
		if(new_offs == cell_offs_)
			return;
		cell_offs_ = new_offs;
		// copy cell data from sgrid to local buffer
		for(ulong i = 0; i < cell_buf_sz; ++i)
			data_[i] = *(h_->sgrid + (cell_offs_ + h_->vidx[i]));
	}

	// subscript comes from base buffer class

	// assign only iterators belonging to the same mesh!
	// or, more precisely, to mesh with identical dimensions
	void assign(const sgrid_ti_traits_impl& rhs) {
		h_ = rhs.h_;
		cell_offs_ = rhs.cell_offs_;
		// assign data buffer
		base_t::assign(rhs);
	}

	ulong backend_index(const ulong offs) const {
		// return index of current iterator inside structured grid buffer
		return cell_offs_ + h_->vidx[offs];
	}

	// don't copy handle in every tops_iterator -- store refernce to it
	// because sgrid_handle for all iterators lives in trimesh
	// and all iterators are invalid when trimesh is destroyed
	//
	// update: storing references doen't work with google btree associative
	// containers, because they construct elements using empty ctor + assignment operator
	// that leads to mesh handle uninitialized in most iterators
	// so we actually have to store a COPY of handle in every sgrid_traits
	//
	// update1: store pointer to sgrid_handle
	const sgrid_handle* h_;
	ulong cell_offs_;
};

// classic sgrid
template< class iterator_type, uint D >
struct sgrid_ti_traits : public sgrid_ti_traits_impl< iterator_type, D > {
	typedef sgrid_ti_traits_impl< iterator_type, D > impl_t;
	typedef typename impl_t::sgrid_handle sgrid_handle;

	sgrid_ti_traits() : impl_t() {}

	sgrid_ti_traits(const sgrid_handle& h) : impl_t(h) {}
};

// rectangular grid
template< class iterator_type, uint D >
struct rgrid_ti_traits : public sgrid_ti_traits_impl< iterator_type, D, rgrid_handle > {
	typedef sgrid_ti_traits_impl< iterator_type, D, rgrid_handle > impl_t;
	typedef typename impl_t::sgrid_handle sgrid_handle;

	rgrid_ti_traits() : impl_t() {}

	rgrid_ti_traits(const sgrid_handle& h) : impl_t(h) {}
};

}} /* blue_sky::wpi */

#endif /* end of include guard: TOPS_ITERATOR_TRAITS_M6AE7JT3 */

