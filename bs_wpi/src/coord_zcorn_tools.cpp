/// @file coord_zcorn_tools.cpp
/// @brief COORD & ZCORN generation and refinement tools
/// @author uentity
/// @date 2011-05-06
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.
#include "coord_zcorn_tools.h"
//#include "mesh_grdecl.h"
#include "tops_iterator.h"
#include "i_cant_link_2_mesh.h"

#define BOUND_MERGE_THRESHOLD 0.8
#define DEFAULT_SMOOTH_RATIO 0.1

/*-----------------------------------------------------------------
 * coord_zcorn_tools implementation
 *----------------------------------------------------------------*/
namespace blue_sky { namespace coord_zcorn_tools {

// other typedefs
typedef bs_array< fp_stor_t, bs_vector_shared > fp_storvec_t;
typedef smart_ptr< fp_storvec_t > spfp_storvec_t;
typedef v_long int_arr_t;
typedef spv_long spi_arr_t;
typedef std::set< fp_t > fp_set;

typedef fp_storarr_t::iterator a_iterator;
typedef fp_storarr_t::const_iterator ca_iterator;
typedef fp_storvec_t::iterator v_iterator;
typedef fp_storvec_t::const_iterator cv_iterator;
typedef fp_set::iterator fps_iterator;

/*-----------------------------------------------------------------
 * helper iterators
 *----------------------------------------------------------------*/
// iterator that jumps over given offset instead of fixed +1
template< class iterator_t, int step_size = 1 >
class slice_iterator : public std::iterator<
						typename std::iterator_traits< iterator_t >::iterator_category,
						typename std::iterator_traits< iterator_t >::value_type,
						typename std::iterator_traits< iterator_t >::difference_type,
						typename std::iterator_traits< iterator_t >::pointer,
						typename std::iterator_traits< iterator_t >::reference
						>
{
	typedef std::iterator_traits< iterator_t > traits_t;

public:
	typedef typename traits_t::value_type value_type;
	typedef typename traits_t::pointer pointer;
	typedef typename traits_t::reference reference;
	typedef typename traits_t::difference_type difference_type;

	// set step in compile-time
	slice_iterator() : p_(), step_(step_size) {}
	slice_iterator(const iterator_t& i) : p_(i), step_(step_size) {}
	slice_iterator(const slice_iterator& i) : p_(i.p_), step_(i.step_) {}

	// set step in runtime
	slice_iterator(difference_type step) : p_(), step_(step) {}
	slice_iterator(const iterator_t& i, difference_type step) : p_(i), step_(step) {}

	reference operator[](difference_type i) const {
		return p_[i * step_];
	}

	void operator+=(difference_type n) {
		p_ += n * step_;
	}
	void operator-=(difference_type n) {
		p_ -= n * step_;
	}

	friend slice_iterator operator+(const slice_iterator& lhs, difference_type n) {
		return slice_iterator(lhs.p_ + (n * lhs.step_), lhs.step_);
	}
	friend slice_iterator operator+(difference_type n, const slice_iterator& lhs) {
		return slice_iterator(lhs.p_ + (n * lhs.step_), lhs.step_);
	}
	friend slice_iterator operator-(const slice_iterator& lhs, difference_type n) {
		return slice_iterator(lhs.p_ - (n * lhs.step_), lhs.step_);
	}
	friend difference_type operator-(const slice_iterator& lhs, const slice_iterator& rhs) {
		return (lhs.p_ - rhs.p_) / lhs.step_;
	}

	slice_iterator& operator++() {
		p_ += step_;
		return *this;
	}
	slice_iterator operator++(int) {
		slice_iterator tmp = *this;
		p_ += step_;
		return tmp;
	}

	slice_iterator& operator--() {
		p_ -= step_;
		return *this;
	}
	slice_iterator operator--(int) {
		slice_iterator tmp = *this;
		p_ -= step_;
		return tmp;
	}

	pointer operator->() {
		return p_.operator->();
	}
	reference operator*() {
		return *p_;
	}

	slice_iterator& operator=(const slice_iterator& lhs) {
		p_ = lhs.p_;
		return *this;
	}

	slice_iterator& operator=(const iterator_t& lhs) {
		p_ = lhs;
		return *this;
	}

	bool operator<(const slice_iterator& rhs) {
		return p_ < rhs.p_;
	}
	bool operator>(const slice_iterator& rhs) {
		return p_ > rhs.p_;
	}
	bool operator<=(const slice_iterator& rhs) {
		return p_ <= rhs.p_;
	}
	bool operator>=(const slice_iterator& rhs) {
		return p_ >= rhs.p_;
	}
	bool operator==(const slice_iterator& rhs) {
		return p_ == rhs.p_;
	}
	bool operator!=(const slice_iterator& rhs) {
		return p_ != rhs.p_;
	}

	iterator_t& backend() {
		return p_;
	}
	const iterator_t& backend() const {
		return p_;
	}

	difference_type step() const {
		return step_;
	}

private:
	iterator_t p_;
	const difference_type step_;
};
typedef slice_iterator< a_iterator, 6 > dim_iterator;
typedef slice_iterator< ca_iterator, 6 > cdim_iterator;

// iterator that calc sumulative sum when doing *p
template< class iterator_t >
class cumsum_iterator : public std::iterator< std::bidirectional_iterator_tag, typename std::iterator_traits< iterator_t >::value_type > {
	typedef std::iterator_traits< cumsum_iterator > traits_t;

public:
	typedef typename traits_t::value_type value_type;
	typedef typename traits_t::pointer pointer;
	typedef typename traits_t::reference reference;
	typedef typename traits_t::difference_type difference_type;

	// set step in compile-time
	cumsum_iterator(value_type offset = 0) : p_(), sum_(offset) {}
	cumsum_iterator(const iterator_t& i, value_type offset = 0) : p_(i), sum_(offset) {}
	cumsum_iterator(const cumsum_iterator& i) : p_(i.p_), sum_(i.sum_) {}

	cumsum_iterator& operator++() {
		sum_ += *p_++;
		return *this;
	}
	cumsum_iterator& operator++(int) {
		cumsum_iterator tmp = *this;
		sum_ += *p_++;
		return tmp;
	}

	cumsum_iterator& operator--() {
		sum_ -= *p_--;
		return *this;
	}
	cumsum_iterator& operator--(int) {
		cumsum_iterator tmp = *this;
		sum_ -= *p_--;
		return tmp;
	}

	pointer operator->() const {
		return p_.operator->();
	}
	value_type operator*() const {
		return sum_;
	}

	cumsum_iterator& operator=(const cumsum_iterator& lhs) {
		p_ = lhs.p_;
		sum_ = lhs.sum_;
		return *this;
	}

	friend difference_type operator-(const cumsum_iterator& lhs, const cumsum_iterator& rhs) {
		return lhs.p_ - rhs.p_;
	}

	bool operator<(const cumsum_iterator& rhs) {
		return p_ < rhs.p_;
	}
	bool operator>(const cumsum_iterator& rhs) {
		return p_ > rhs.p_;
	}
	bool operator<=(const cumsum_iterator& rhs) {
		return p_ <= rhs.p_;
	}
	bool operator>=(const cumsum_iterator& rhs) {
		return p_ >= rhs.p_;
	}
	bool operator==(const cumsum_iterator& rhs) {
		return p_ == rhs.p_;
	}
	bool operator!=(const cumsum_iterator& rhs) {
		return p_ != rhs.p_;
	}

	iterator_t& backend() {
		return p_;
	}
	const iterator_t& backend() const {
		return p_;
	}

private:
	iterator_t p_;
	value_type sum_;
};

/*-----------------------------------------------------------------
 * helper structure to get dx[i] regardless of whether it given by one number or by array of
 * numbers
 *----------------------------------------------------------------*/
struct dim_subscript {
	dim_subscript(const fp_storarr_t& dim, fp_t offset = .0)
		: dim_(dim), offset_(offset), sum_(offset)
	{
		if(dim_.size() == 1)
			ss_fcn_ = &dim_subscript::ss_const_dim;
		else
			ss_fcn_ = &dim_subscript::ss_array_dim;
	}

	fp_stor_t ss_const_dim(const t_ulong idx) {
		return static_cast< fp_stor_t >(fp_t(offset_ + dim_[0] * idx));
	}

	fp_stor_t ss_array_dim(const t_ulong idx) {
		if(idx > 0) sum_ += dim_[idx - 1];
		return static_cast< fp_stor_t>(sum_);
	}

	void reset() { sum_ = offset_; }

	fp_t operator[](t_ulong idx) {
		return (this->*ss_fcn_)(idx);
	}

private:
	const fp_storarr_t& dim_;
	fp_stor_t (dim_subscript::*ss_fcn_)(const t_ulong);
	const fp_t offset_;
	fp_t sum_;
};

/*-----------------------------------------------------------------
 * misc functions to help generating COORD & ZCORN
 *----------------------------------------------------------------*/
template< class array_t >
void resize_zcorn(array_t& zcorn, int_t nx, int_t ny, int_t new_nx, int_t new_ny) {
	using namespace std;
	typedef typename array_t::iterator v_iterator;
	typedef typename array_t::value_type value_t;
	typedef typename array_t::size_type size_t;

	const int_t nz = (int_t)(zcorn.size() >> 3) / (nx  * ny);
	//const int_t delta = 4 * (new_nx * new_ny - nx * ny);

	// cache z-values for each plane
	const int_t plane_sz = nx * ny * 4;
	vector< value_t > z_cache(nz * 2);
	slice_iterator< v_iterator > pz(zcorn.begin(), plane_sz);
	slice_iterator< v_iterator > pz_end(zcorn.begin(), plane_sz);
	copy(pz, pz + nz *2, z_cache.begin());

	// refill zcorn with cached values & new plane size
	const int_t new_plane_sz = (new_nx * new_ny * 4);
	zcorn.resize(new_nx * new_ny * nz * 8);
	v_iterator p_newz = zcorn.begin();
	for(typename vector< value_t >::const_iterator p_cache = z_cache.begin(), end = z_cache.end(); p_cache != end; ++p_cache) {
		fill_n(p_newz, new_plane_sz, *p_cache);
  p_newz += new_plane_sz;
	}
}

template< class ret_array_t, class array_t >
void coord2deltas(const array_t& src, ret_array_t& res) {
	typedef typename array_t::const_iterator carr_iterator;

	if(src.size() < 2) return;
	res.resize(src.size() - 1);
	typename ret_array_t::iterator p_res = res.begin();
	carr_iterator a = src.begin(), b = a, end = src.end();
	for(++b; b != end; ++b)
		*p_res++ = *b - *a++;
}

spfp_storarr_t gen_coord(int_t nx, int_t ny, spfp_storarr_t dx, spfp_storarr_t dy, fp_t x0, fp_t y0) {
	using namespace std;
	typedef fp_storarr_t::value_type value_t;

	// DEBUG
	BSOUT << "gen_coord: init stage" << bs_end;
	// create subscripters
	dim_subscript dxs(*dx, x0);
	dim_subscript dys(*dy, y0);

	// if dimension offset is given as array, then size should be taken from array size
	if(dx->size() > 1) nx = (int_t) dx->size();
	if(dy->size() > 1) ny = (int_t) dy->size();

	// create arrays
	spfp_storarr_t coord = BS_KERNEL.create_object(fp_storarr_t::bs_type());
// FIXME: raise exception
	if(!coord) return NULL;

	// DEBUG
	//BSOUT << "gen_coord: creation starts..." << bs_end;
	// fill coord
	// coord is simple grid
	coord->resize((nx + 1)*(ny + 1)*6, value_t(0));
	fp_storarr_t::iterator pcd = coord->begin();
	for(int_t iy = 0; iy <= ny; ++iy) {
		fp_t cur_y = dys[iy];
		for(int_t ix = 0; ix <= nx; ++ix) {
			pcd[0] = pcd[3] = dxs[ix];
			pcd[1] = pcd[4] = cur_y;
			pcd[5] = 1; // pcd[2] = 0 from init
			pcd += 6;
		}
		dxs.reset();
	}
	// DEBUG
	//BSOUT << "gen_coord: creation finished" << bs_end;

	return coord;
}

spfp_storarr_t gen_coord2(spfp_storarr_t x, spfp_storarr_t y) {
	using namespace std;
	typedef fp_storarr_t::value_type value_t;
  int_t nx, ny;
  value_t *ys, *xs;

	// DEBUG
	//BSOUT << "gen_coord: init stage" << bs_end;

	// if dimension offset is given as array, then size should be taken from array size
	if(x->size() > 1) nx = (int_t) x->size() - 1;
	if(y->size() > 1) ny = (int_t) y->size() - 1;

  xs = &(*x)[0];
  ys = &(*y)[0];


	// create arrays
	spfp_storarr_t coord = BS_KERNEL.create_object(fp_storarr_t::bs_type());
// FIXME: raise exception
	if(!coord) return NULL;

	// DEBUG
	//BSOUT << "gen_coord: creation starts..." << bs_end;
	// fill coord
	// coord is simple grid
	coord->resize((nx + 1)*(ny + 1)*6, value_t(0));
	fp_storarr_t::iterator pcd = coord->begin();
	for(int_t iy = 0; iy <= ny; ++iy) {
		fp_t cur_y = ys[iy];
		for(int_t ix = 0; ix <= nx; ++ix) {
			pcd[0] = pcd[3] = xs[ix];
			pcd[1] = pcd[4] = cur_y;
			pcd[5] = 1; // pcd[2] = 0 from init
			pcd += 6;
		}
	}
	// DEBUG
	//BSOUT << "gen_coord: creation finished" << bs_end;

	return coord;
}

coord_zcorn_pair gen_coord_zcorn(int_t nx, int_t ny, int_t nz, spv_float dx, spv_float dy, spv_float dz, fp_stor_t x0, fp_stor_t y0, fp_stor_t z0) {
	using namespace std;
	typedef coord_zcorn_pair ret_t;
	typedef v_float::value_type value_t;
    spv_float null_arr = 0;

	// DEBUG
	//BSOUT << "gen_coord_zcorn: init stage" << bs_end;
	// create subscripter
	if(!dx || !dy || !dz) return ret_t(null_arr, null_arr);
	if(!dx->size() || !dy->size() || !dz->size()) return ret_t(null_arr, null_arr);

	// if dimension offset is given as array, then size should be taken from array size
	if(dx->size() > 1) nx = (int_t) dx->size();
	if(dy->size() > 1) ny = (int_t) dy->size();
	if(dz->size() > 1) nz = (int_t) dz->size();

	// create zcorn array
	spv_float zcorn = BS_KERNEL.create_object(v_float::bs_type());
	// FIXME: raise exception
	if(!zcorn) return ret_t(null_arr, null_arr);

	// fill zcorn
	// very simple case
	dim_subscript dzs(*dz, z0);
	zcorn->resize(nx * ny * nz * 8);

	// DEBUG
	//BSOUT << "gen_coord_zcorn: ZCORN creating starts..." << bs_end;
	v_float::iterator pcd = zcorn->begin();
	const t_long plane_size = nx * ny * 4;
	t_float z_cache = dzs[0];
	for(t_long iz = 1; iz <= nz; ++iz) {
		fill_n(pcd, plane_size, z_cache);
		pcd += plane_size;
		z_cache = dzs[iz];
		fill_n(pcd, plane_size, z_cache);
		pcd += plane_size;
	}
	// DEBUG
	//BSOUT << "gen_coord_zcorn: ZCORN creating finished" << bs_end;
	//BSOUT << "gen_coord_zcorn: COORD creating starts..." << bs_end;

	return ret_t(gen_coord(nx, ny, dx, dy, x0, y0), zcorn);
}

struct extract_plane {
	uint_t nx, ny;

	extract_plane(uint_t nx_, uint_t ny_) : nx(nx_), ny(ny_) {};

	template< class input_iter, class outp_iter >
	void copy_row(input_iter& start, outp_iter& tgt) {
		typedef slice_iterator< input_iter > input_slice;
		typedef slice_iterator< outp_iter > outp_slice;

		// some consts
		const uint_t row_sz = nx * 24;
		const uint_t srow_sz = (nx + 1)*3;

		input_slice pin(start, 24);
		outp_slice pout(tgt, 3);
		// add top left corner of every cell
		// x
		copy(pin, pin + nx, pout);
		// y
		pin = start + 1; pout.backend() = tgt + 1;
		copy(pin, pin + nx, pout);
		// z
		pin = start + 2; pout = tgt + 2;
		copy(pin, pin + nx, pout);

		// last point is top right corner
		tgt = std::copy(start + row_sz - 21, start + row_sz - 18, tgt + srow_sz - 3);
		//tgt += srow_sz;
		start += row_sz;
	}

	template< class input_iter, class outp_iter >
	void go(input_iter& start, outp_iter& tgt) {
		// copy first ny rows from plane
		for(t_ulong y = 0; y < ny; ++y)
			copy_row(start, tgt);

		// last row is point 2 on cell cube
		// don't mess original start
		input_iter start_ = start - nx*24 + 2*3;
		copy_row(start_, tgt);
	}
};

// convert tops array to structured grid representation
template< class src_iterator_t >
spv_float tops2struct_grid_impl(const uint_t nx, const uint_t ny, const uint_t nz, src_iterator_t tops) {
	typedef v_float::iterator v_iterator;
	//typedef v_float::const_iterator cv_iterator;
	using namespace std;

	// sanity check
	if(!nx || !ny)
		return spv_float();

	const uint_t plane_sz = nx * ny * 24;
	//const uint_t nz = uint_t(tops->size() / plane_sz);
	const uint_t splane_sz = (nx + 1)*(ny + 1)*3;

	// resulting array
	spv_float res = BS_KERNEL.create_object(v_float::bs_type());
	res->resize(splane_sz * (nz + 1));
	v_iterator tgt = res->begin();
	src_iterator_t start = tops;

	// loop over z layers
	extract_plane ex(nx, ny);
	for(t_ulong i = 0; i < nz; ++i)
		ex.go(start, tgt);
	// last plane is bootom of the whole cube
	start -= plane_sz - 4*3;
	ex.go(start, tgt);

	return res;
}

// specialization for plain tops array
spv_float tops2struct_grid(const uint_t nx, const uint_t ny, spv_float tops) {
	return tops2struct_grid_impl(
		nx, ny, uint_t(tops->size() / (nx * ny * 24)), tops->begin()
	);
}

// direct conversion from COORD & ZCORN using tops_iterator
spv_float tops2struct_grid(uint_t nx, uint_t ny, spv_float coord, spv_float zcorn) {
	typedef smart_ptr< rs_smesh_iface, true > sp_smesh;
	typedef wpi::tops_iterator< wpi::carray_ti_traits, 3 > iterator_t;

	// build mesh_grdecl around given mesh
	sp_himesh handy = BS_KERNEL.create_object("handy_mesh_iface");
	sp_smesh mesh = handy->make_mesh_grdecl(nx, ny, coord, zcorn);
	assert(mesh);

	return tops2struct_grid_impl(
		nx, ny, uint_t(mesh->get_n_elements() / (nx * ny)), iterator_t(mesh)
	);
}

spv_float tops2struct_grid(smart_ptr< rs_smesh_iface > mesh) {
	if(!mesh) return NULL;
	typedef wpi::tops_iterator< wpi::carray_ti_traits, 3 > iterator_t;

	// obtain mesh dimensions
	rs_smesh_iface::index_point3d_t dims = mesh->get_dimens();
	// go
	return tops2struct_grid_impl(
		uint_t(dims[0]), uint_t(dims[1]), uint_t(dims[2]), iterator_t(mesh.get())
	);
}

/*-----------------------------------------------------------------
 * helper structure related to first refine_mesh algorithm
 *----------------------------------------------------------------*/
struct proc_ray {
	template< int_t direction, class = void >
	struct dir_ray {
		enum { dir = direction };
		bool is_bound;

		dir_ray(bool is_bound_ = false) : is_bound(is_bound_) {}

		template< class ray_t >
		static typename ray_t::iterator end(ray_t& ray) {
			return ray.end();
		}

		template< class ray_t >
		static typename ray_t::iterator last(ray_t& ray) {
			return --ray.end();
		}

		template< class ray_t >
		static typename ray_t::iterator closest_bound(ray_t& ray, fp_t v) {
			return std::upper_bound(ray.begin(), ray.end(), v);
		}

		static fp_t min(fp_t f, fp_t s) {
			return std::min(f, s);
		}

		static fp_t max(fp_t f, fp_t s) {
			return std::max(f, s);
		}
	};

	template< class unused >
	struct dir_ray< -1, unused > {
		enum { dir = -1 };
		bool is_bound;

		dir_ray(bool is_bound_ = false) : is_bound(is_bound_) {}

		template< class ray_t >
		static typename ray_t::iterator end(ray_t& ray) {
			return --ray.begin();
		}

		template< class ray_t >
		static typename ray_t::iterator last(ray_t& ray) {
			return ray.begin();
		}

		template< class ray_t >
		static typename ray_t::iterator closest_bound(ray_t& ray, fp_t v) {
			typename ray_t::iterator t = std::upper_bound(ray.begin(), ray.end(), v);
			if(t == ray.begin()) return t;
			return --t;
		}

		static fp_t min(fp_t f, fp_t s) {
			return std::max(f, s);
		}

		static fp_t max(fp_t f, fp_t s) {
			return std::min(f, s);
		}
	};

	typedef const fp_t& predicate_fcn_t(const fp_t&, const fp_t&);
	template< class ray_t >
	static fp_t find_cell(ray_t& coord, predicate_fcn_t p) {
		typedef typename ray_t::iterator ray_iterator;
		if(coord.size() < 2) return 0;
		// find max cell size
		fp_t cell_sz = 0;
		ray_iterator a = coord.begin(), b = a, end = coord.end();
		for(++b; b != end; ++b) {
			cell_sz = p(*b - *a++, cell_sz);
			//if(tmp > max_cell) max_cell = tmp;
		}
		return cell_sz;
	}

	template< class ray_t >
	static fp_t find_max_cell(ray_t& coord) {
		return find_cell< ray_t > (coord, std::max< fp_t >);
	}

	template< class ray_t >
	static fp_t find_min_cell(ray_t& coord) {
		return find_cell< ray_t >(coord, std::min< fp_t >);
	}

	template< class ray_t, class dir_ray_t >
	static void go(ray_t& ray, fp_t start_point, fp_stor_t d, fp_stor_t a, dir_ray_t dr) {
		using namespace std;
		const int_t dir = static_cast< int_t >(dir_ray_t::dir);
		// find where new point fall to
		const fp_t max_sz = find_max_cell(ray) * 0.5;
		const fp_t max_front = *dr.last(ray);
		fp_t cell_sz = d;
		fp_t wave_front = dr.min(start_point + dir * 0.5 * d, max_front);

		// make refined ray
		// add refinement grid
		fp_set ref_ray;
		while(cell_sz <= max_sz && abs(wave_front - max_front) > 0) {
			ref_ray.insert(wave_front);
			cell_sz *= a;
			wave_front = dr.min(wave_front + dir * cell_sz, max_front);
		}

		// merge refinement with original grid
		copy(ray.begin(), ray.end(), insert_iterator< fp_set >(ref_ray, ref_ray.begin()));
		// copy results back
		ray.clear();
		copy(ref_ray.begin(), ref_ray.end(), insert_iterator< ray_t >(ray, ray.begin()));
	}

	template< class ray_t >
	static int_t kill_tight_cells(ray_t& ray, fp_t min_cell_sz) {
		typedef typename ray_t::iterator ray_iterator;

		if(ray.size() < 2) return 0;
		int_t merge_cnt = 0;
		ray_iterator a = ray.begin(), b = a, end = ray.end();
		for(++b; b != end; ++b) {
			if(*b - *a < min_cell_sz) {
				ray.erase(a++);
				++merge_cnt;
			}
			else ++a;
		}
		return merge_cnt;
	}

	template< class ray_t >
	static int_t band_filter(ray_t& ray, fp_t smooth_ratio) {
		typedef typename ray_t::size_type size_t;
		typedef std::set< size_t, std::greater< size_t > > idx_set;
		typedef typename idx_set::const_iterator idx_iterator;

		if(ray.size() < 2) return 0;

		// find cells to remove
		idx_set dying;
		//std::vector< size_t > dying;
		//
		// left to right walk
		int_t n = 0;
		for(size_t i = 0; i < ray.size() - 1; ++i) {
			// dont check dead bands
			if(dying.find(i) != dying.end()) continue;
			// main condition
			if(ray[i] * smooth_ratio > ray[i + 1]) {
				dying.insert(i + 1);
				//dying.push_back(i + 1);
				// check if cell after bound is less than before bound
				if(i + 2 < ray.size() && ray[i + 2] < ray[i])
					ray[i + 2] += ray[i + 1];
				else
					ray[i] += ray[i + 1];
				++n;
			}
		}

		// right to left walk
		for(size_t i = ray.size() - 1; i >= 1; --i) {
			// dont check dead bands
			if(dying.find(i) != dying.end()) continue;
			// main condition
			if(ray[i] * smooth_ratio > ray[i - 1]) {
				dying.insert(i - 1);
				//dying.push_back(i + 1);
				// check if cell after bound is less than before bound
				if(i - 2 < ray.size() && ray[i - 2] < ray[i])
					ray[i - 2] += ray[i - 1];
				else
					ray[i] += ray[i - 1];
				++n;
			}
		}

		// remove dead cells
		for(idx_iterator p_id = dying.begin(), end = dying.end(); p_id != end; ++p_id) {
			ray.erase(ray.begin() + *p_id);
		}
		//for(size_t i = 0; i < dying.size(); ++i)
		//	ray.erase(ray.begin() + i);
		return n;
	}
};

template< class ray_t >
void refine_mesh_impl(ray_t& coord, fp_stor_t point, fp_stor_t d, fp_stor_t a) {
	using namespace std;
	typedef typename fp_storvec_t::iterator v_iterator;
	typedef typename fp_storvec_t::const_iterator cv_iterator;

	// sanity check
	if(!coord.size()) return;

	// process coord in both directions
	proc_ray::go(coord, point, d, a, typename proc_ray::template dir_ray< 1 >());
	proc_ray::go(coord, point, d, a, typename proc_ray::template dir_ray< -1 >());
}

/*-----------------------------------------------------------------
 * convert points given in (i, j) format to absolute fp coordinates
 *----------------------------------------------------------------*/
spfp_storarr_t point_index2coord(int_t nx, int_t ny, spfp_storarr_t coord, spi_arr_t points_pos,
		bool relative = false)
{
	using namespace std;
	typedef slice_iterator< v_iterator, 6 > dim_iterator;
	typedef int_arr_t::value_type pos_t;
	const int_t ydim_step = 6 * (nx + 1);

	// (left, top)
	dim_iterator px(coord->begin());
	dim_iterator py(coord->begin() + 1, ydim_step);
	const fp_stor_t x0 = *min_element(px, px + (nx + 1));
	const fp_stor_t y0 = *min_element(py, py + (ny + 1));

	// create resulting array
	spfp_storarr_t res = BS_KERNEL.create_object(fp_storarr_t::bs_type());
	if(!res) return res;
	const int_t points_num = points_pos->size() >> 1;
	res->resize(points_num * 2);

	int_arr_t::const_iterator p = points_pos->begin();
	fp_storarr_t::iterator r = res->begin();
	pos_t i, j;
	for(int_t t = 0; t < points_num; ++t) {
		// points = {(i, j)}
		i = *p++; j = *p++;

		// find x coordinate
		dim_iterator px = coord->begin();
		advance(px, min(nx, i));
		*r = (*(px + 1) + *px) * 0.5;
		if(relative) *r -= x0;
		++r;

		// find y coordinate
		dim_iterator py(coord->begin() + 1, ydim_step);
		advance(py, min(ny, j));
		*r = (*(py + 1) + *py) * 0.5;
		if(relative) *r -= y0;
		++r;
	}

	return res;
}

spfp_storarr_t point_index2coord(int_t nx, int_t ny, spfp_storarr_t coord, spi_arr_t points_pos, spfp_storarr_t points_param) {
	using namespace std;
	typedef slice_iterator< v_iterator, 6 > dim_iterator;
	typedef int_arr_t::value_type pos_t;
	const int_t ydim_step = 6 * (nx + 1);

	// create resulting array
	spfp_storarr_t res = BS_KERNEL.create_object(fp_storarr_t::bs_type());
	if(!res) return res;
	const int_t points_num = points_pos->size() >> 1;
	res->resize(points_num * 6);

	int_arr_t::const_iterator p = points_pos->begin();
	fp_storarr_t::const_iterator pp = points_param->begin();
	fp_storarr_t::iterator r = res->begin();
	pos_t i, j;
	for(int_t t = 0; t < points_num; ++t) {
		// points = {(i, j)}
		i = *p++; j = *p++;

		// find x coordinate
		dim_iterator px = coord->begin();
		advance(px, min(nx, i));
		*r++ = (*(px + 1) + *px) * 0.5;

		// find y coordinate
		dim_iterator py(coord->begin() + 1, ydim_step);
		advance(py, min(ny, j));
		*r++ = (*(py + 1) + *py) * 0.5;

		// fill other points params (dx, dy, ax, ay)
		copy(pp, pp + 4, r);
	}

	return res;
}

/*-----------------------------------------------------------------
 * new mesh generation algo based on well points
 *----------------------------------------------------------------*/
template< class ray_t, class dir_ray_t >
void make_wave(ray_t& ray, fp_t start_point, fp_stor_t d, fp_stor_t a,
		fp_stor_t max_sz, fp_stor_t field_len, fp_stor_t closest, dir_ray_t dr)
{
	using namespace std;
	const int_t dir = static_cast< int_t >(dir_ray_t::dir);
	// DEBUG
	//BSOUT << dir << ": center " << start_point << "; closest " << closest << bs_end;

	// try to calc, how many wave fronts can we insert
	const uint_t N_max = uint_t(floor(std::log(max_sz / d) / std::log(a) + 1));
	fp_t S = (dir * (closest - start_point) + d) * 0.5;
	uint_t N = uint_t(floor(std::log(S * (a - 1)/d + 1) / std::log(a)));
	N = min< uint_t >(N, N_max);
	// insert one bound anyway
	//N = max< uint_t>(N, 1);
	// calc tail to half of distance to nearest bound
	fp_t tail = S - d * (std::pow(a, double(N)) - 1)/(a - 1);

	// refined ray stored here
	fp_set ref_ray;
	// farthest bound
	const fp_t max_front = dr.max(0, field_len);
	fp_t cell_sz = d;
	fp_t wave_front = dr.min(start_point - dir * 0.5 * d, max_front);

	if(tail < d * 0.5 && N > 0) {
		--N;
		if(S >= d && dir > 0 && !dr.is_bound)
			ref_ray.insert(wave_front + S);
	}

	for(uint_t i = 0; i < N; ++i) {
		wave_front += dir * cell_sz;
		//if(i == N -1 && tail < d * 0.5)
		//	wave_front += dir * tail;
		wave_front = dr.min(wave_front, max_front);
		if(abs(wave_front - max_front) < 1e-10)
			break;
		ref_ray.insert(wave_front);
		cell_sz = min< fp_t >(cell_sz * a, max_sz);
	}

	// merge refinement with original grid - suppose that ray is a set
	copy(ref_ray.begin(), ref_ray.end(), insert_iterator< ray_t >(ray, ray.begin()));
	//copy(ray.begin(), ray.end(), insert_iterator< fp_set >(ref_ray, ref_ray.begin()));
	//// copy results back
	//ray.clear();
	//copy(ref_ray.begin(), ref_ray.end(), insert_iterator< ray_t >(ray, ray.begin()));
}

template< class delta_t >
void fill_gaps(delta_t& d, fp_stor_t cell_sz, fp_stor_t min_sz,
		fp_t max_sz_tol = 0.3, bool strict_max_sz = false)
{
	using namespace std;
	typedef typename delta_t::iterator d_iterator;
	// fill big gaps with cells of given size

	list< fp_stor_t > refined_d;
	const fp_t m = 1 / cell_sz;
	fp_t cur_sz;
	for(d_iterator pd = d.begin(), end = d.end(); pd != end; ++pd) {
		if(*pd <= cell_sz) {
			refined_d.push_back(*pd);
			continue;
		}
		// we have long gap here
		// how much cells can we insert?
		uint_t N = uint_t(floor(*pd * m));
		fp_t tail = *pd - N * cell_sz;

		// by default insert cells with size = cell_sz
		cur_sz = cell_sz;
		// if we can't make cells larger than cell_sz
		// or have only one cell
		// then increase N and make equal regular cells < cell_sz
		if(N == 1 || strict_max_sz) {
			++N;
			cur_sz = *pd / N;
			tail = 0;
		}
		else if(tail > max_sz_tol * cell_sz && tail >= min_sz) {
			// non-strict mode
			// if tail is relatively big - push it as separate cell
			refined_d.push_back(tail);
			tail = 0;
		}
		// fill gap
		for(uint_t i = 0; i < N; ++i) {
			if((i == 0 || i == N - 1) && tail > 0) {
				// spread short tail between first & last cell
				refined_d.push_back(cur_sz + tail * 0.5);
			}
			else
				refined_d.push_back(cur_sz);
		}
	}

	// copy refined delta back
	d.resize(refined_d.size());
	copy(refined_d.begin(), refined_d.end(), d.begin());
}

template< class delta_t, class hit_idx_t >
void find_hit_idx(
	const delta_t& delta_x, const delta_t& delta_y, hit_idx_t& hit_idx, spfp_storarr_t points_pos,
	fp_stor_t x0 = 0, fp_stor_t y0 = 0)
{
	using namespace std;
	typedef typename hit_idx_t::iterator hit_iterator;
	typedef typename delta_t::const_iterator delta_iterator;
	typedef cumsum_iterator< delta_iterator > cs_iterator;
	typedef typename cs_iterator::difference_type diff_t;

	uint_t cnt = points_pos->size() >> 1;
	hit_idx.resize(cnt * 2);
	hit_iterator p_hit = hit_idx.begin();
	a_iterator pp = points_pos->begin();
	for(uint_t i = 0; i < cnt; ++i) {
		cs_iterator p_id = lower_bound(
			cs_iterator(delta_x.begin(), x0),
			cs_iterator(delta_x.end(), x0), *pp++);
		*p_hit++ = max< diff_t >(p_id - delta_x.begin() - 1, 0);
		p_id = lower_bound(
			cs_iterator(delta_y.begin(), y0),
			cs_iterator(delta_y.end(), y0), *pp++);
		*p_hit++ = max< diff_t >(p_id - delta_y.begin() - 1, 0);
	}
}

// here we assume that coord are nondescending in X and Y directions
template< class coord_t, class hit_idx_t >
void find_hit_idx(uint_t nx, uint_t ny, const coord_t& coord,
	hit_idx_t& hit_idx, spfp_storarr_t points_pos)
{
	using namespace std;
	typedef typename hit_idx_t::iterator hit_iterator;
	typedef typename cdim_iterator::difference_type diff_t;
	const int_t ydim_step = 6 * (nx + 1);

	uint_t cnt = points_pos->size() >> 1;
	hit_idx.resize(cnt * 2);
	hit_iterator p_hit = hit_idx.begin();
	ca_iterator pp = points_pos->begin();
	for(uint_t i = 0; i < cnt; ++i) {
		fp_t crd_x = *pp++, crd_y = *pp++;
		// find x id
		cdim_iterator p_xid = lower_bound(
			cdim_iterator(coord.begin()),
			cdim_iterator(coord.begin()) + (nx + 1), crd_x);
		*p_hit++ = max< diff_t >((p_xid - coord.begin()) - 1, 0);

		// find y id
		// lower bound don't work here because step is set in runtime
		// so use stupid search
		cdim_iterator p_yid(coord.begin() + 1, ydim_step);
		for(uint_t i = 0; i <= ny; ++i) {
			if(*p_yid >= crd_y) break;
			++p_yid;
		}
		*p_hit++ = max< diff_t >((p_yid - (coord.begin() + 1)) - 1, 0);
	}
}

// x & y coord vectors define initial mesh,
// wich will be refined by wave algorithm
coord_zcorn_pair wave_mesh_deltas_s1_impl(
	fp_stor_t max_dx, fp_stor_t max_dy,
	spfp_storarr_t points_pos, spfp_storarr_t points_param,
	const fp_set& x, const fp_set& y)
{
	using namespace std;
	typedef fp_storvec_t::iterator v_iterator;
	typedef std::map< fp_stor_t, fp_storarr_t::const_iterator > pmap_t;

	// helper to produce waves in + and - directions from give point
	struct make_2side_wave {
		static void go(const pmap_t& pmap, fp_set& mesh, fp_stor_t len, fp_stor_t max_d, uint cid) {
			static char dir[] = {'x', 'y'};
			// store processed points here
			fp_set p_ready;

			// process points in X direction
			int_t cnt = 0;
			pmap_t::const_iterator lower = pmap.begin();
			pmap_t::const_iterator upper = pmap.begin();
			++upper;
			for(pmap_t::const_iterator p = pmap.begin(), end = pmap.end(); p != end; ++p) {
				fp_storarr_t::const_iterator p_param = p->second;
				fp_stor_t d = p_param[cid];
				fp_stor_t a = p_param[cid + 2];

				// point coord value
				fp_stor_t val = p->first;
				// DEBUG
				BSOUT << "point[" << ++cnt << "] at (" << dir[cid] << " = " << val << "), d" << dir[cid] << " = " << d
					<< ", a" << dir[cid] << " = " << a << bs_end;
				// process only new points
				if(d != 0 && p_ready.find(val) == p_ready.end()) {
					make_wave(mesh, val, d, a, max_d, len,
							upper == end ? len + (len - val) : upper->first,
							proc_ray::dir_ray< 1 >(upper == end));
					make_wave(mesh, val, d, a, max_d, len,
							p == pmap.begin() ? -val : lower->first,
							proc_ray::dir_ray< -1 >(p == pmap.begin()));

					p_ready.insert(val);
				}
				if(p != pmap.begin())
					++lower;
#ifdef _WIN32
				if(upper != end)
#endif
					++upper;
			}
		}
	};

	// DEBUG
	BSOUT << "wave_mesh_deltas: init stage" << bs_end;
	// sanity check
	if(!points_pos || x.size() < 2 || y.size() < 2) return coord_zcorn_pair();

	fp_storarr_t::const_iterator pp = points_pos->begin(); //, p_end = points_pos->end();
	const ulong p_num = points_pos->size() >> 1;
	//p_end -= (p_end - pp) % 2;

	// prepare x & y coord vectors to be without offset
	const fp_t x_offs = *x.begin();
	const fp_t y_offs = *y.begin();
	fp_set x0;
	// x0[i] = x[i] - x_offs
	transform(x.begin(), x.end(), insert_iterator< fp_set >(x0, x0.begin()),
		bind2nd(std::minus< fp_t >(), x_offs));
	// y0[i] = y[i] - y_offs
	fp_set y0;
	transform(y.begin(), y.end(), insert_iterator< fp_set >(y0, y0.begin()),
		bind2nd(std::minus< fp_t >(), y_offs));
	const fp_t len_x = *(--x0.end());
	const fp_t len_y = *(--y0.end());

	// DEBUG
	BSOUT << "wave_mesh_deltas: points processing starts..." << bs_end;
	BSOUT << "len_x = " << len_x << ", len_y = " << len_y << bs_end;

	// params array: {(dx, dy, ax, ay)}
	// if params specified only once - then params is equal for all points
	fp_stor_t dx, dy, ax, ay;
	bool const_params = false;
	fp_storarr_t::const_iterator p_param = points_param->begin();
	if(points_param->size() == 4) {
		const_params = true;
		dx = *(p_param++); dy = *(p_param++);
		ax = *(p_param++); ay = *(p_param++);
		p_param = points_param->begin();
	}
	else if(points_param->size() < p_num * 4) {
		BSERR << "wave_mesh_deltas_s1: wrong size of points params array!" << bs_end;
		return coord_zcorn_pair();
	}

	// points array: {(x, y}}
	// first pass - build set of increasing px_coord & py_coord
	// for each point save pointer to point parameters
	//fp_set px_coord, py_coord;
	pmap_t px_coord, py_coord;
	for(ulong i = 0; i < p_num; ++i) {
		px_coord[*pp++] = p_param;
		py_coord[*pp++] = p_param;
		if(!const_params)
			p_param += 4;
		//px_coord.insert(*pp++);
		//py_coord.insert(*pp++);
	}

	// process points in X direction
	make_2side_wave::go(px_coord, x0, len_x, max_dx, 0);

	// process points in Y direction
	make_2side_wave::go(py_coord, y0, len_y, max_dy, 1);

	// make deltas from coordinates
	//vector< fp_stor_t > delta_x, delta_y;
	//coord2deltas(x0, delta_x);
	//coord2deltas(y0, delta_y);

	// copy delta_x & delta_y to bs_arrays
	//nx = (int_t)  delta_x.size();
	//ny = (int_t)  delta_y.size();
	spfp_storarr_t adx = BS_KERNEL.create_object(fp_storarr_t::bs_type()),
	               ady = BS_KERNEL.create_object(fp_storarr_t::bs_type());
	coord2deltas(x0, *adx);
	coord2deltas(y0, *ady);
	//adx->resize(delta_x.size()); ady->resize(delta_y.size());
	//copy(delta_x.begin(), delta_x.end(), adx->begin());
	//copy(delta_y.begin(), delta_y.end(), ady->begin());

	return coord_zcorn_pair(adx, ady);
}

// std wave mesh implementation
coord_zcorn_pair wave_mesh_deltas_s1(
	fp_stor_t max_dx, fp_stor_t max_dy,
	fp_stor_t len_x, fp_stor_t len_y, spfp_storarr_t points_pos, spfp_storarr_t points_param)
{
	// make intial mesh with bounds
	fp_set x, y;
	x.insert(0); x.insert(len_x);
	y.insert(0); y.insert(len_y);
	// call implementation
	return wave_mesh_deltas_s1_impl(max_dx, max_dy, points_pos, points_param, x, y);
}

// mesh refine with wave algorithm
// returns refined coord & zcorn
coord_zcorn_pair refine_wave_mesh(
	int_t& nx, int_t& ny,
	spfp_storarr_t coord, spfp_storarr_t zcorn,
	fp_stor_t max_dx, fp_stor_t max_dy,
	spfp_storarr_t points_pos, spfp_storarr_t points_param)
{
	using namespace std;
	// convert COORD -> x & y coord vectors
	fp_set x;
	dim_iterator px(coord->begin());
	copy(px, px + (nx + 1), insert_iterator< fp_set >(x, x.begin()));

	fp_set y;
	dim_iterator py(coord->begin() + 1, 6 * (nx + 1));
	copy(py, py + (ny + 1), insert_iterator< fp_set >(y, y.begin()));

	// remember old dimensions for correct zcorn resize
	int_t old_nx = nx;
	int_t old_ny = ny;

	// (left, top, up)
	const fp_stor_t x0 = *x.begin();
	const fp_stor_t y0 = *y.begin();
	//const fp_stor_t z0 = *min_element(zcorn->begin(), zcorn->end());

	// invoke wave algo on existing mesh
	coord_zcorn_pair deltas = wave_mesh_deltas_s1_impl(max_dx, max_dy, points_pos, points_param, x, y);
	spfp_storarr_t& dx = deltas.first;
	spfp_storarr_t& dy = deltas.second;
	nx = dx->size(); ny = dy->size();

	// generate new coord
	spfp_storarr_t new_coord = gen_coord(nx, ny, dx, dy, x0, y0);

	// rzcorn = copy of zcorn
	spfp_storarr_t new_zcorn = BS_KERNEL.create_object(fp_storarr_t::bs_type());
	new_zcorn->resize(zcorn->size());
	std::copy(new_zcorn->begin(), new_zcorn->end(), zcorn->begin());
	// resize rzcorn
	resize_zcorn(*new_zcorn, old_nx, old_ny, nx, ny);

	// return result
	return coord_zcorn_pair(new_coord, new_zcorn);
}

// mesh refine with wave algorithm
// returns refined coord & zcorn
coord_zcorn_pair refine_wave_mesh_deltas(
	spfp_storarr_t dx, spfp_storarr_t dy,
	fp_stor_t max_dx, fp_stor_t max_dy,
	spfp_storarr_t points_pos, spfp_storarr_t points_param)
{
	using namespace std;
	typedef cumsum_iterator< fp_storarr_t::iterator > cs_iterator;

	// convert deltas -> coord vectors
	fp_set x;
	copy(cs_iterator(dx->begin()), cs_iterator(dx->end()),
		insert_iterator< fp_set >(x, x.begin()));
	// append last bound
	x.insert(*(--x.end()) + dx->ss(dx->size() - 1));

	fp_set y;
	copy(cs_iterator(dy->begin()), cs_iterator(dy->end()),
		insert_iterator< fp_set >(y, y.begin()));
	// append last bound
	y.insert(*(--y.end()) + dy->ss(dy->size() - 1));

	// invoke wave algo on existing mesh
	return wave_mesh_deltas_s1_impl(max_dx, max_dy, points_pos, points_param, x, y);
}

void wave_mesh_deltas_s2(
	fp_stor_t cell_dx, fp_stor_t cell_dy,
	fp_stor_t min_dx, fp_stor_t min_dy,
	spfp_storarr_t dx, spfp_storarr_t dy,
	fp_t max_sz_tol, bool strict_max_sz)
{
	//BSOUT << "fill gaps" << bs_end;
	fill_gaps(*dx, cell_dx, min_dx, max_sz_tol, strict_max_sz);
	fill_gaps(*dy, cell_dy, min_dy, max_sz_tol, strict_max_sz);
}

void wave_mesh_deltas_s2(
	fp_stor_t cell_dx, fp_stor_t cell_dy,
	spfp_storarr_t points_param,
	spfp_storarr_t dx, spfp_storarr_t dy,
	fp_t max_sz_tol, bool strict_max_sz)
{
	using namespace std;

	// find min_dx & min_dy based on params
	typedef slice_iterator< a_iterator, 4 > param_iterator;
	param_iterator wx = points_param->begin();
	param_iterator wy = points_param->begin() + 1;
	uint_t param_len = points_param->size() >> 2;

	wave_mesh_deltas_s2(
		cell_dx, cell_dy,
		*min_element(wx, wx + param_len),
		*min_element(wy, wy + param_len),
		dx, dy,
		max_sz_tol, strict_max_sz
	);
}

coord_zcorn_pair wave_mesh_deltas(
	fp_stor_t max_dx, fp_stor_t max_dy,
	fp_stor_t len_x, fp_stor_t len_y,
	spfp_storarr_t points_pos, spfp_storarr_t points_param,
	spi_arr_t hit_idx)
{
	using namespace std;
	// call stage 1
	coord_zcorn_pair res = wave_mesh_deltas_s1(
		max_dx, max_dy, len_x, len_y,
		points_pos, points_param
	);
	// call stage 2
	wave_mesh_deltas_s2(
		max_dx, max_dy, points_param,
		res.first, res.second
	);
	// find hit index if needed
	if(hit_idx)
		find_hit_idx(*res.first, *res.second, *hit_idx, points_pos);

	// DEBUG
	// check if sum(deltas) = len
	BSOUT << "sum(delta_x) = " << accumulate(
		res.first->begin(), res.first->end(), fp_stor_t(0)) << bs_end;
	BSOUT << "sum(delta_y) = " << accumulate(
		res.second->begin(), res.second->end(), fp_stor_t(0)) << bs_end;

	return res;
}

/*-----------------------------------------------------------------
 * wave_mesh = wave_mesh_deltas + gen_coord_zcorn
 *----------------------------------------------------------------*/
coord_zcorn_pair wave_mesh(
	int_t& nx, int_t& ny, fp_stor_t max_dx, fp_stor_t max_dy,
	fp_stor_t len_x, fp_stor_t len_y, spfp_storarr_t points_pos, spfp_storarr_t points_param,
	int_t nz, spfp_storarr_t dz,
	fp_stor_t x0, fp_stor_t y0, fp_stor_t z0)
{
	// refine coord
	coord_zcorn_pair deltas = wave_mesh_deltas(
		max_dx, max_dy, len_x, len_y,
		points_pos, points_param
	);
	spfp_storarr_t& dx = deltas.first;
	spfp_storarr_t& dy = deltas.second;
	nx = dx->size(); ny = dy->size();

	return gen_coord_zcorn(nx, ny, nz, dx, dy, dz, x0, y0, z0);
}

/*-----------------------------------------------------------------
 * wave_mesh for points in (i,j) format
 *----------------------------------------------------------------*/
coord_zcorn_pair wave_mesh(
	spfp_storarr_t coord,
	int_t& nx, int_t& ny, fp_stor_t max_dx, fp_stor_t max_dy,
	spi_arr_t points_pos, spfp_storarr_t points_param,
	int_t nz, spfp_storarr_t dz,
	fp_stor_t x0, fp_stor_t y0, fp_stor_t z0)
{
	using namespace std;
	// calc field length
	const int_t ydim_step = 6 * (nx + 1);
	fp_stor_t len_x, len_y;
	dim_iterator px(coord->begin());
	dim_iterator py(coord->begin() + 1, ydim_step);
	len_x = *max_element(px, px + (nx + 1))
		- *min_element(px, px + (nx + 1));
	len_y = *max_element(py, py + (ny + 1))
		- *min_element(py, py + (ny + 1));
	// invoke mesh gen
	return wave_mesh(
		nx, ny, max_dx, max_dy, len_x, len_y,
		point_index2coord(nx, ny, coord, points_pos, true),
		points_param, nz, dz, x0, y0, z0
	);
}

/*-----------------------------------------------------------------
 * refine mesh algo based on existing grid
 *----------------------------------------------------------------*/
coord_zcorn_pair refine_mesh_deltas(int_t& nx, int_t& ny, spfp_storarr_t coord,
	spfp_storarr_t points, fp_t cell_merge_thresh, fp_t band_thresh,
	spi_arr_t hit_idx)
{
	using namespace std;

	typedef fp_storvec_t::iterator v_iterator;
	typedef slice_iterator< v_iterator, 6 > dim_iterator;

	// DEBUG
	BSOUT << "refine_mesh: init stage" << bs_end;
	// sanity check
	if(!coord || !points) return coord_zcorn_pair();

	// convert coord & zcorn to shared vectors
	//spfp_storvec_t vcoord = BS_KERNEL.create_object(fp_storvec_t::bs_type());
	//if(vcoord) vcoord->init_inplace(coord->get_container());
	//else return coord_zcorn_pair();

	//spfp_storvec_t vzcorn = BS_KERNEL.create_object(fp_storvec_t::bs_type());
	//if(vzcorn) vzcorn->init_inplace(zcorn->get_container());
	//else return coord_zcorn_pair();

	// build x and y coord maps
	//spfp_storvec_t x = BS_KERNEL.create_object(fp_storvec_t::bs_type());
	fp_set x;
	copy(dim_iterator(coord->begin()), dim_iterator(coord->begin()) + (nx + 1),
			insert_iterator< fp_set >(x, x.begin()));

	//spfp_storvec_t y = BS_KERNEL.create_object(fp_storvec_t::bs_type());
	//y->resize(ny + 1);
	fp_set y;
	const int_t ydim_step = 6 * (nx + 1);
	copy(dim_iterator(coord->begin() + 1, ydim_step), dim_iterator(coord->begin() + 1, ydim_step) + (ny + 1),
			insert_iterator< fp_set >(y, y.begin()));

	fp_storarr_t::const_iterator pp = points->begin(), p_end = points->end();
	// make (p_end - p_begin) % 4 = 0
	p_end -= (p_end - pp) % 6;

	// DEBUG
	BSOUT << "refine_mesh: points processing starts..." << bs_end;
	// process points in turn
	// points array: {(x, y, dx, dy, ax, ay)}
	fp_stor_t x_coord, y_coord, dx, dy, ax, ay;
	fp_stor_t min_dx = 0, min_dy = 0;
	bool first_point = true;
	// store processed points here
	fp_set dx_ready, dy_ready;
	// main cycle
	int_t cnt = 0;
	while(pp != p_end) {
		x_coord = *(pp++); y_coord = *(pp++);
		dx = *(pp++); dy = *(pp++);
		ax = *(pp++); ay = *(pp++);
		if(first_point) {
			min_dx = dx; min_dy = dy;
			first_point = false;
		}
		else {
			min_dx = min(min_dx, dx);
			min_dy = min(min_dx, dy);
		}
		// DEBUG
		BSOUT << "point[" << ++cnt << "] at (" << x_coord << ", " << y_coord << "), dx = " << dx
		<< ", dy = " << dy << ", ax = " << ax << ", ay = " << ay << bs_end;
		// process only new points
		if(dx != 0 && dx_ready.find(x_coord) == dx_ready.end()) {
			refine_mesh_impl(x, x_coord, dx, ax);
			dx_ready.insert(x_coord);
		}
		if(dy != 0 && dy_ready.find(y_coord) == dy_ready.end()) {
			refine_mesh_impl(y, y_coord, dy, ay);
			dy_ready.insert(y_coord);
		}
	}

	// DEBUG
	//BSOUT << "refine_mesh: kill_tight_centers(x)" << bs_end;
	// kill too tight cells in x & y directions
	while(proc_ray::kill_tight_cells(x, cell_merge_thresh * min_dx)) {}
	// DEBUG
	//BSOUT << "refine_mesh: kill_tight_centers(y)" << bs_end;
	while(proc_ray::kill_tight_cells(y, cell_merge_thresh * min_dy)) {}

	// DEBUG
	//BSOUT << "refine_mesh: coord2deltas" << bs_end;
	// make deltas from coordinates
	vector< fp_stor_t > delta_x, delta_y;
	coord2deltas(x, delta_x);
	coord2deltas(y, delta_y);

	// DEBUG
	BSOUT << "refine_mesh: bands filter" << bs_end;
	// filter tight bands from grid
	while(proc_ray::band_filter(delta_x, band_thresh)) {}
	while(proc_ray::band_filter(delta_y, band_thresh)) {}

	// find coord mesh offsets
	fp_stor_t x0 = *min_element(x.begin(), x.end());
	fp_stor_t y0 = *min_element(y.begin(), y.end());
	// find what cells in refined mesh are hit by given points
	if(hit_idx) {
		typedef cumsum_iterator< vector< fp_stor_t >::iterator > cs_iterator;
		typedef cs_iterator::difference_type diff_t;

		hit_idx->resize(cnt * 2);
		int_arr_t::iterator p_hit = hit_idx->begin();
		pp = points->begin();
		for(int_t i = 0; i < cnt; ++i) {
			cs_iterator p_id = lower_bound(cs_iterator(delta_x.begin(), x0), cs_iterator(delta_x.end(), x0), *pp++);
			*p_hit++ = max< diff_t >(p_id - delta_x.begin() - 1, 0);
			p_id = lower_bound(cs_iterator(delta_y.begin(), y0), cs_iterator(delta_y.end(), y0), *pp++);
			*p_hit++ = max< diff_t >(p_id - delta_y.begin() - 1, 0);
			pp += 4;
		}
	}

	// DEBUG
	//BSOUT << "refine_mesh: copy delta_x & delta_y to bs_arrays" << bs_end;
	// copy delta_x & delta_y to bs_arrays
	nx = (int_t)  delta_x.size();
	ny = (int_t)  delta_y.size();
	spfp_storarr_t adx = BS_KERNEL.create_object(fp_storarr_t::bs_type()),
				   ady = BS_KERNEL.create_object(fp_storarr_t::bs_type());
	adx->resize(delta_x.size()); ady->resize(delta_y.size());
	copy(delta_x.begin(), delta_x.end(), adx->begin());
	copy(delta_y.begin(), delta_y.end(), ady->begin());

	return coord_zcorn_pair(adx, ady);
}

coord_zcorn_pair refine_mesh(int_t& nx, int_t& ny, spfp_storarr_t coord, spfp_storarr_t zcorn,
		spfp_storarr_t points, fp_t cell_merge_thresh, fp_t band_thresh,
		spi_arr_t hit_idx)
{
	using namespace std;

	if(!zcorn) return coord_zcorn_pair();

	// remember old dimesnions for correct zcorn resize
	int_t old_nx = nx;
	int_t old_ny = ny;

	// refine coord
	coord_zcorn_pair refine_deltas = refine_mesh_deltas(nx, ny, coord, points, cell_merge_thresh, band_thresh, hit_idx);
	spfp_storarr_t& delta_x = refine_deltas.first;
	spfp_storarr_t& delta_y = refine_deltas.second;

	// DEBUG
	BSOUT << "refine_mesh: update ZCORN" << bs_end;
	// update zcorn
	vector< fp_stor_t > vzcorn(zcorn->begin(), zcorn->end());
	resize_zcorn(vzcorn, old_nx, old_ny, nx, ny);
	// create bs_array from new zcorn
	spfp_storarr_t rzcorn = BS_KERNEL.create_object(fp_storarr_t::bs_type());
	if(!rzcorn) return coord_zcorn_pair();
	rzcorn->resize(vzcorn.size());
	copy(vzcorn.begin(), vzcorn.end(), rzcorn->begin());

	// rebuild grid based on processed x_coord & y_coord
	return coord_zcorn_pair(gen_coord(nx, ny, delta_x, delta_y, coord->ss(0), coord->ss(1)), rzcorn);
}

/*-----------------------------------------------------------------
 * refine_mesh_deltas for points given in (i,j) format
 *----------------------------------------------------------------*/
coord_zcorn_pair refine_mesh_deltas(int_t& nx, int_t& ny, spfp_storarr_t coord,
	spi_arr_t points_pos, spfp_storarr_t points_param, fp_t cell_merge_thresh, fp_t band_thresh,
	spi_arr_t hit_idx)
{
	return refine_mesh_deltas(
		nx, ny, coord,
		point_index2coord(nx, ny, coord, points_pos, points_param),
		cell_merge_thresh, band_thresh, hit_idx
	);
}

/*-----------------------------------------------------------------
 * refine_mesh for points given in (i,j) format
 *----------------------------------------------------------------*/
coord_zcorn_pair refine_mesh(int_t& nx, int_t& ny, spfp_storarr_t coord, spfp_storarr_t zcorn,
		spi_arr_t points_pos, spfp_storarr_t points_param, fp_t cell_merge_thresh, fp_t band_thresh,
		spi_arr_t hit_idx)
{
	return refine_mesh(
		nx, ny, coord, zcorn,
		point_index2coord(nx, ny, coord, points_pos, points_param),
		cell_merge_thresh, band_thresh, hit_idx
	);
}

spi_arr_t find_hit_idx(
	spfp_storarr_t dx, spfp_storarr_t dy, spfp_storarr_t points_pos,
	fp_stor_t x0, fp_stor_t y0)
{
	spi_arr_t hit_idx = BS_KERNEL.create_object(int_arr_t::bs_type());
	find_hit_idx(*dx, *dy, *hit_idx, points_pos, x0, y0);
	return hit_idx;
}

spi_arr_t find_hit_idx(
	uint_t nx, uint_t ny, spfp_storarr_t coord,
	spfp_storarr_t points_pos)
{
	spi_arr_t hit_idx = BS_KERNEL.create_object(int_arr_t::bs_type());
	find_hit_idx(nx, ny, *coord, *hit_idx, points_pos);
	return hit_idx;
}

/*-----------------------------------------------------------------
 * generate structured grid with given parameters
 *----------------------------------------------------------------*/
spv_float gen_sgrid(
	t_ulong nx, t_ulong ny, t_ulong nz,
	spv_float dx, spv_float dy, spv_float dz,
	t_float x0, t_float y0, t_float z0,
	bool zyx_order
) {
	enum { D = 3 };

	dim_subscript ds[] = {
		dim_subscript(*dx, x0),
		dim_subscript(*dy, y0),
		dim_subscript(*dz, z0)
	};
	spv_float deltas[] = { dx, dy, dz };
	uint_t sz[] = { nx + 1, ny + 1, nz + 1 };
	t_ulong sgrid_sz = 1;
	for(uint_t i = 0; i < D; ++i) {
		if(deltas[i]->size() > 1)
			sz[i] = deltas[i]->size();
		sgrid_sz *= sz[i];
	}

	// reorder if specified
	t_uint order[] = { 0, 1, 2 };
	if(zyx_order) {
		std::swap(order[0], order[2]);
		//std::reverse(&ds[0], &ds[3]);
		//std::reverse(&sz[0], &sz[3]);
	}

	// generate mesh
	spv_float sgrid = BS_KERNEL.create_object(v_float::bs_type());
	if(!sgrid) return NULL;
	sgrid->resize(sgrid_sz * D);
	v_float::iterator psg = sgrid->begin();

	t_float vertex[] = { 0, 0, 0 };
	for(t_ulong i = 0; i < sz[order[2]]; ++i) {
		vertex[order[2]] = ds[order[2]][i];
		for(t_ulong j = 0; j < sz[order[1]]; ++j) {
			vertex[order[1]] = ds[order[1]][j];
			for(t_ulong k = 0; k < sz[order[0]]; ++k) {
				vertex[order[0]] = ds[order[0]][k];
				// in ZYX order we'll get [z, y, x] coords order in vertex
				//if(zyx_order)
				//	std::reverse(&vertex[0], &vertex[D]);
				// store vertex value
				psg = std::copy(&vertex[0], &vertex[D], psg);
			}
			ds[order[0]].reset();
		}
		ds[order[1]].reset();
	}
	return sgrid;
}

spv_float gen_sgrid_2d(
	t_ulong nx, t_ulong ny,
	spv_float dx, spv_float dy,
	t_float x0, t_float y0,
	bool yx_order
) {
	enum { D = 2 };

	dim_subscript ds[] = {
		dim_subscript(*dx, x0),
		dim_subscript(*dy, y0)
	};
	spv_float deltas[] = { dx, dy };
	uint_t sz[] = { nx + 1, ny + 1 };
	t_ulong sgrid_sz = 1;
	for(uint_t i = 0; i < D; ++i) {
		if(deltas[i]->size() > 1)
			sz[i] = deltas[i]->size() + 1;
		sgrid_sz *= sz[i];
	}

	// reorder if specified
	t_uint order[] = { 0, 1 };
	if(yx_order) {
		std::swap(order[0], order[1]);
	}

	// generate mesh
	spv_float sgrid = BS_KERNEL.create_object(v_float::bs_type());
	if(!sgrid) return NULL;
	sgrid->resize(sgrid_sz * 3);
	v_float::iterator psg = sgrid->begin();

	t_float vertex[] = { 0, 0, 0 };
	for(t_ulong j = 0; j < sz[order[1]]; ++j) {
		vertex[order[1]] = ds[order[1]][j];
		for(t_ulong k = 0; k < sz[order[0]]; ++k) {
			vertex[order[0]] = ds[order[0]][k];
			// in ZYX order we'll get [z, y, x] coords order in vertex
			//if(zyx_order)
			//	std::reverse(&vertex[0], &vertex[D]);
			// store vertex value
			psg = std::copy(&vertex[0], &vertex[3], psg);
		}
		ds[order[0]].reset();
	}
	return sgrid;
}

}}	// eof blue_sky::coord_zcorn_tools

