/// @file coord_zcorn_tools.cpp
/// @brief COORD & ZCORN generation and refinement tools
/// @author uentity
/// @date 2011-05-06
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "coord_zcorn_tools.h"
#include "mesh_grdecl.h"

#define BOUND_MERGE_THRESHOLD 0.8
#define DEFAULT_SMOOTH_RATIO 0.1

/*-----------------------------------------------------------------
 * coord_zcorn_tools implementation
 *----------------------------------------------------------------*/
namespace blue_sky { namespace coord_zcorn_tools {

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
 * convert points given in (i, j) format to absolute fp coordinates
 *----------------------------------------------------------------*/
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

BS_API_PLUGIN coord_zcorn_pair refine_mesh_deltas_s(
	int_t& nx, int_t& ny, fp_stor_t max_dx, fp_stor_t max_dy,
	fp_stor_t len_x, fp_stor_t len_y, spfp_storarr_t points_pos, spfp_storarr_t points_param)
{
	using namespace std;

	typedef fp_storvec_t::iterator v_iterator;
	typedef slice_iterator< v_iterator, 6 > dim_iterator;

	// DEBUG
	BSOUT << "refine_mesh: init stage" << bs_end;
	// sanity check
	if(!points_pos) return coord_zcorn_pair();

	fp_storarr_t::const_iterator pp = points_pos->begin(), p_end = points_pos->end();
	// make (p_end - p_begin) % 4 = 0
	p_end -= (p_end - pp) % 2;

	// DEBUG
	BSOUT << "refine_mesh: points processing starts..." << bs_end;
	BSOUT << "len_x = " << len_x << ", len_y = " << len_y << bs_end;
	// points array: {(x, y}}
	// first pass - build set of increasing px_coord & py_coord
	fp_set px_coord, py_coord;
	while(pp != p_end) {
		px_coord.insert(*pp++);
		py_coord.insert(*pp++);
	}

	// make intial mesh with bounds
	fp_set x, y;
	x.insert(0); x.insert(len_x);
	y.insert(0); y.insert(len_y);

	// if params specified only once - then params is equal for all points
	fp_stor_t dx, dy, ax, ay;
	fp_storarr_t::const_iterator p_param = points_param->begin();
	bool const_params = false;
	if(points_param->size() == 4) {
		const_params = true;
		dx = *(p_param++); dy = *(p_param++);
		ax = *(p_param++); ay = *(p_param++);
	}

	// points coord
	//fp_stor_t x_coord = 0, y_coord = 0;
	// store processed points here
	fp_set dx_ready, dy_ready;

	// process points in X direction
	int_t cnt = 0;
	fp_set::const_iterator lower = px_coord.begin();
	fp_set::const_iterator upper = px_coord.begin();
	++upper;
	for(fp_set::const_iterator px = px_coord.begin(), end = px_coord.end(); px != end; ++px) {
		if(!const_params) {
			dx = *(p_param++); dy = *(p_param++);
			ax = *(p_param++); ay = *(p_param++);
		}
		// DEBUG
		BSOUT << "point[" << ++cnt << "] at (x = " << *px << "), dx = " << dx
		<< ", ax = " << ax << bs_end;
		// process only new points
		if(dx != 0 && dx_ready.find(*px) == dx_ready.end()) {
			refine_point_s(x, *px, dx, ax, max_dx, len_x,
				upper == end ? len_x + (len_x - *px) : *upper,
				typename proc_ray::dir_ray< 1 >());
			refine_point_s(x, *px, dx, ax, max_dx, len_x,
				px == px_coord.begin() ? -*px : *lower,
				typename proc_ray::dir_ray< -1 >());

			//refine_mesh_impl_s(x, *px, dx, ax, max_dx, len_x);
			dx_ready.insert(*px);
		}
		if(px != px_coord.begin())
			++lower;
		++upper;
	}

	// process points in Y direction
	cnt = 0;
	p_param = points_param->begin();
	lower = py_coord.begin();
	upper = py_coord.begin(); ++upper;
	for(fp_set::const_iterator py = py_coord.begin(), end = py_coord.end(); py != end; ++py) {
		if(!const_params) {
			dx = *(p_param++); dy = *(p_param++);
			ax = *(p_param++); ay = *(p_param++);
		}
		// DEBUG
		BSOUT << "point[" << ++cnt << "] at (y = " << *py << "), dy = " << dy
		<< ", ay = " << ay << bs_end;
		// process only new points
		if(dy != 0 && dy_ready.find(*py) == dy_ready.end()) {
			refine_point_s(y, *py, dy, ay, max_dy, len_y,
				upper == end ? len_y + (len_y - *py) : *upper,
				typename proc_ray::dir_ray< 1 >());
			refine_point_s(y, *py, dy, ay, max_dy, len_y,
				py == py_coord.begin() ? -*py : *lower,
				typename proc_ray::dir_ray< -1 >());

			//refine_mesh_impl_s(y, *py, dy, ay, max_dy, len_y);
			dy_ready.insert(*py);
		}
		if(py != py_coord.begin())
			++lower;
		++upper;
	}

	//proc_ray::kill_tight_cells(x, dx);
	//proc_ray::kill_tight_cells(y, dy);

	// DEBUG
	//BSOUT << "refine_mesh: coord2deltas" << bs_end;
	// make deltas from coordinates
	vector< fp_stor_t > delta_x, delta_y;
	coord2deltas(x, delta_x);
	coord2deltas(y, delta_y);

	// DEBUG
	// check if sum(deltas) = len
	//BSOUT << "sum(delta_x) = " << accumulate(delta_x.begin(), delta_x.end(), fp_stor_t(0)) << bs_end;
	//BSOUT << "sum(delta_y) = " << accumulate(delta_y.begin(), delta_y.end(), fp_stor_t(0)) << bs_end;

	//BSOUT << "fill gaps" << bs_end;
	//fill_gaps(delta_x, max_dx);
	//fill_gaps(delta_y, max_dy);

	// DEBUG
	// check if sum(deltas) = len
	BSOUT << "sum(delta_x) = " << accumulate(delta_x.begin(), delta_x.end(), fp_stor_t(0)) << bs_end;
	BSOUT << "sum(delta_y) = " << accumulate(delta_y.begin(), delta_y.end(), fp_stor_t(0)) << bs_end;


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

}}

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

/*-----------------------------------------------------------------
 * points given in abs coordinates merged with params
 *----------------------------------------------------------------*/
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

/*-----------------------------------------------------------------
 * points given in (i,j) format
 *----------------------------------------------------------------*/
std::pair< spv_float, spv_float >
mesh_grdecl::refine_mesh_deltas(t_long& nx, t_long& ny, spv_float coord,
		spv_long points_pos, spv_float points_param, spv_long hit_idx,
		t_double cell_merge_thresh, t_double band_thresh)
{
	return czt::refine_mesh_deltas(nx, ny, coord, points_pos, points_param,
			cell_merge_thresh, band_thresh, hit_idx);
}

std::pair< spv_float, spv_float >
mesh_grdecl::refine_mesh(t_long& nx, t_long& ny, spv_float coord, spv_float zcorn,
		spv_long points_pos, spv_float points_param, spv_long hit_idx,
		t_double cell_merge_thresh, t_double band_thresh)
{
	return czt::refine_mesh(nx, ny, coord, zcorn, points_pos, points_param,
			cell_merge_thresh, band_thresh, hit_idx);
}

