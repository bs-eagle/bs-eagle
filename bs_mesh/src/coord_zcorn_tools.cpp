#include "bs_mesh_stdafx.h"
#include "mesh_grdecl.h"
#include <iterator>

#define BOUND_MERGE_THRESHOLD 0.8
#define DEFAULT_SMOOTH_RATIO 0.1

/*-----------------------------------------------------------------
 * helpers for gen_coord_zcorn & refine_mesh
 *----------------------------------------------------------------*/
namespace blue_sky { namespace coord_zcorn_tools {

// shorter aliases
typedef t_long int_t;
typedef t_double fp_t;
typedef t_float fp_stor_t;

// other typedefs
typedef bs_array< fp_stor_t, bs_vector_shared > fp_storvec_t;
typedef v_float fp_storarr_t;
typedef v_long int_arr_t;
typedef smart_ptr< fp_storvec_t > spfp_storvec_t;
typedef spv_float spfp_storarr_t;
typedef spv_long spi_arr_t;
typedef std::pair< spfp_storarr_t, spfp_storarr_t > coord_zcorn_pair;
typedef std::set< fp_t > fp_set;

typedef fp_storvec_t::iterator v_iterator;
typedef fp_storvec_t::const_iterator cv_iterator;
typedef fp_set::iterator fps_iterator;

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

	slice_iterator& operator++() {
		p_ += step_;
		return *this;
	}
	slice_iterator& operator++(int) {
		slice_iterator tmp = *this;
		p_ += step_;
		return tmp;
	}

	slice_iterator& operator--() {
		p_ -= step_;
		return *this;
	}
	slice_iterator& operator--(int) {
		slice_iterator tmp = *this;
		p_ -= step_;
		return tmp;
	}

	pointer operator->() const {
		return p_.operator->();
	}
	reference operator*() const {
		return *p_;
	}

	slice_iterator& operator=(const slice_iterator& lhs) {
		p_ = lhs.p_;
		return *this;
	}

	friend difference_type operator-(const slice_iterator& lhs, const slice_iterator& rhs) {
		return (lhs.p_ - rhs.p_) / lhs.step_;
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

private:
	iterator_t p_;
	const difference_type step_;
};

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

template< class array_t >
static fp_t sum(const array_t& a) {
	fp_t s = 0;
	for(typename array_t::const_iterator i = a.begin(), end = a.end(); i != end; ++i)
		s += *i;
	return s;
}

// helper structure to get dx[i] regardless of whether it given by one number or by array of
// numbers
struct dim_subscript {
	dim_subscript(const fp_storarr_t& dim, fp_t offset = .0)
		: dim_(dim), offset_(offset), sum_(offset)
	{
		if(dim_.size() == 1)
			ss_fcn_ = &dim_subscript::ss_const_dim;
		else
			ss_fcn_ = &dim_subscript::ss_array_dim;
	}

	fp_stor_t ss_const_dim(int_t idx) {
		return static_cast< fp_stor_t >(fp_t(dim_[0] * idx));
	}

	fp_stor_t ss_array_dim(int_t idx) {
		fp_t tmp = sum_;
		if(idx > 0) sum_ += dim_[idx - 1];
		return static_cast< fp_stor_t>(tmp);
	}

	void reset() { sum_ = offset_; }

	fp_t operator[](int_t idx) {
		return (this->*ss_fcn_)(idx);
	}

private:
	const fp_storarr_t& dim_;
	fp_stor_t (dim_subscript::*ss_fcn_)(int_t);
	const fp_t offset_;
	fp_t sum_;
};

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
	BSOUT << "gen_coord: creation starts..." << bs_end;
	// fill coord
	// coord is simple grid
	coord->init((nx + 1)*(ny + 1)*6, value_t(0));
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
	BSOUT << "gen_coord: creation finished" << bs_end;

	return coord;
}

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

void insert_column(int_t nx, int_t ny, fp_storvec_t& coord, fp_storvec_t& zcorn, fp_stor_t where) {
	using namespace std;

	typedef fp_storvec_t::iterator v_iterator;
	typedef slice_iterator< v_iterator, 6 > dim_iterator;

	// reserve mem for insterts
	coord.reserve((nx + 2)*(ny + 1)*6);

	// find a place to insert
	dim_iterator pos = lower_bound(dim_iterator(coord.begin()), dim_iterator(coord.begin()) + (nx + 1), where);
	//if(pos == dim_iterator(coord.begin()) + (nx + 1)) return;
	int_t ins_offset = pos.backend() - coord.begin();

	// process all rows
	fp_stor_t y, z1, z2;
	v_iterator vpos;
	for(int_t i = ny; i >= 0; --i) {
		// save y and z values
		vpos = coord.begin() + i*(nx + 1)*6 + ins_offset;
		if(ins_offset == (nx + 1)*6) {
			// insert at the boundary
			y = *(vpos - 5); z1 = *(vpos - 4); z2 = *(vpos - 1);
		}
		else {
			// insert in the beginning/middle of row
			y = *(vpos + 1); z1 = *(vpos + 2); z2 = *(vpos + 5);
		}
		// insert new vector
		insert_iterator< fp_storvec_t > ipos(coord, vpos);
		*ipos++ = where; *ipos++ = y; *ipos++ = z1;
		*ipos++ = where; *ipos++ = y; *ipos = z2;
	}

	// update zcorn
	resize_zcorn(zcorn, nx, ny, nx + 1, ny);
}

void insert_row(int_t nx, int_t ny, fp_storvec_t& coord, fp_storvec_t& zcorn, fp_stor_t where) {
	using namespace std;
	typedef fp_storvec_t::iterator v_iterator;
	typedef slice_iterator< v_iterator > dim_iterator;
	const int_t ydim_step = 6 * (nx + 1);

	// reserve mem for insterts
	coord.reserve((nx + 1)*(ny + 2)*6);

	// find a place to insert
	const dim_iterator search_end = dim_iterator(coord.begin() + 1, ydim_step) + (ny + 1);
	dim_iterator pos = lower_bound(dim_iterator(coord.begin() + 1, ydim_step), search_end, where);
	//if(pos == search_end) return;
	v_iterator ins_point = pos.backend() - 1;

	// make cache of x values from first row
	spfp_storvec_t cache_x = BS_KERNEL.create_object(fp_storvec_t::bs_type());
	cache_x->resize(nx + 1);
	typedef slice_iterator< v_iterator, 3 > hdim_iterator;
	copy(hdim_iterator(coord.begin()), hdim_iterator(coord.begin() + (nx + 1)*6), cache_x->begin());

	// insert row
	insert_iterator< fp_storvec_t > ipos(coord, ins_point);
	v_iterator p_x = cache_x->begin();
	fp_stor_t z1 = *(coord.begin() + 2), z2 = *(coord.begin() + 5);
	for(int_t i = 0; i <= nx; ++i) {
		*ipos++ = *p_x++; *ipos++ = where; *ipos++ = z1;
		*ipos++ = *p_x++; *ipos++ = where; *ipos++ = z2;
	}

	// update zcorn
	resize_zcorn(zcorn, nx, ny, nx, ny + 1);
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

struct proc_ray {
	template< int_t direction, class = void >
	struct dir_ray {
		enum { dir = direction };

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
	};

	template< class unused >
	struct dir_ray< -1, unused > {
		enum { dir = -1 };

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
			return std::upper_bound(ray.begin(), ray.end(), v)--;
		}

		static fp_t min(fp_t f, fp_t s) {
			return std::max(f, s);
		}
	};

	template< class ray_t, class predicate >
	static fp_t find_cell(ray_t& coord, predicate p) {
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
		return find_cell(coord, std::max< fp_t >);
	}

	template< class ray_t >
	static fp_t find_min_cell(ray_t& coord) {
		return find_cell(coord, std::min< fp_t >);
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

coord_zcorn_pair refine_mesh_deltas(int_t& nx, int_t& ny, spfp_storarr_t coord,
	spfp_storarr_t points, fp_t cell_merge_thresh, fp_t band_thresh,
	spi_arr_t hit_idx = NULL)
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

	// find what cells in refined mesh are hit by given points
	if(hit_idx) {
		typedef cumsum_iterator< vector< fp_stor_t >::iterator > cs_iterator;
		typedef cs_iterator::difference_type diff_t;

		hit_idx->resize(cnt * 2);
		int_arr_t::iterator p_hit = hit_idx->begin();
		pp = points->begin();
		for(int_t i = 0; i < cnt; ++i) {
			cs_iterator p_id = lower_bound(cs_iterator(delta_x.begin(), *x.begin()), cs_iterator(delta_x.end(), *x.begin()), *pp++);
			*p_hit++ = max< diff_t >(p_id - delta_x.begin() - 1, 0);
			p_id = lower_bound(cs_iterator(delta_y.begin(), *y.begin()), cs_iterator(delta_y.end(), *y.begin()), *pp++);
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
		spi_arr_t hit_idx = NULL)
{
	using namespace std;

	if(!zcorn) return coord_zcorn_pair();

	// refine coord
	coord_zcorn_pair refine_deltas = refine_mesh_deltas(nx, ny, coord, points, cell_merge_thresh, band_thresh, hit_idx);
	spfp_storarr_t& delta_x = refine_deltas.first;
	spfp_storarr_t& delta_y = refine_deltas.second;

	// DEBUG
	BSOUT << "refine_mesh: update ZCORN" << bs_end;
	// update zcorn
	vector< fp_stor_t > vzcorn(zcorn->begin(), zcorn->end());
	resize_zcorn(vzcorn, nx, ny, (int_t) delta_x->size(), (int_t) delta_y->size());
	// create bs_array from new zcorn
	spfp_storarr_t rzcorn = BS_KERNEL.create_object(fp_storarr_t::bs_type());
	if(!rzcorn) return coord_zcorn_pair();
	rzcorn->resize(vzcorn.size());
	copy(vzcorn.begin(), vzcorn.end(), rzcorn->begin());

	// rebuild grid based on processed x_coord & y_coord
	return coord_zcorn_pair(gen_coord(nx, ny, delta_x, delta_y, coord->ss(0), coord->ss(1)), rzcorn);
}

}}  // eof namespace blue_sky::coord_zcorn_tools

/*-----------------------------------------------------------------
 * implementation of bs_mesh::gen_coord_zcorn & refine_mesh
 *----------------------------------------------------------------*/
using namespace blue_sky;
namespace czt = blue_sky::coord_zcorn_tools;

std::pair< spv_float, spv_float >
mesh_grdecl::gen_coord_zcorn(t_long nx, t_long ny, t_long nz, spv_float dx, spv_float dy, spv_float dz, t_float x0, t_float y0, t_float z0) {
	using namespace std;
	typedef std::pair< spv_float, spv_float > ret_t;
	typedef v_float::value_type value_t;

	// DEBUG
	BSOUT << "gen_coord_zcorn: init stage" << bs_end;
	// create subscripter
	if(!dx || !dy || !dz) return ret_t(NULL, NULL);
	if(!dx->size() || !dy->size() || !dz->size()) return ret_t(NULL, NULL);

	// if dimension offset is given as array, then size should be taken from array size
	if(dz->size() > 1) nz = (t_long) dz->size();

	// create zcorn array
	spv_float zcorn = BS_KERNEL.create_object(v_float::bs_type());
  // FIXME: raise exception
	if(!zcorn) return ret_t(NULL, NULL);

	// fill zcorn
	// very simple case
	czt::dim_subscript dzs(*dz, z0);
	zcorn->init(nx*ny*nz*8);

	// DEBUG
	BSOUT << "gen_coord_zcorn: ZCORN creating starts..." << bs_end;
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
	BSOUT << "gen_coord_zcorn: ZCORN creating finished" << bs_end;
	BSOUT << "gen_coord_zcorn: COORD creating starts..." << bs_end;

	return ret_t(czt::gen_coord(nx, ny, dx, dy, x0, y0), zcorn);
}

std::pair< spv_float, spv_float >
mesh_grdecl::refine_mesh_deltas(t_long& nx, t_long& ny, spv_float coord,
		spv_float points, spv_long hit_idx, t_double cell_merge_thresh, t_double band_thresh)
{
	return czt::refine_mesh_deltas(nx, ny, coord, points, cell_merge_thresh, band_thresh, hit_idx);
}

std::pair< spv_float, spv_float >
mesh_grdecl::refine_mesh(t_long& nx, t_long& ny, spv_float coord, spv_float zcorn,
		spv_float points, spv_long hit_idx, t_double cell_merge_thresh, t_double band_thresh)
{
	return czt::refine_mesh(nx, ny, coord, zcorn, points, cell_merge_thresh, band_thresh, hit_idx);
}

