/// @file surface_reader.cpp
/// @brief Read Petrel surface ASCII file format into BS array
/// @author uentity
/// @version 1.0
/// @date 30.08.2015
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "bs_kernel.h"
#include "conf.h"
#include "pool_iface.h"

#include <fstream>
#include <sstream>
#include <locale>

namespace blue_sky {
typedef smart_ptr< h5_pool_iface > sp_h5_pool;

// hidden details
namespace {

// convert array from fortran to C order
template< class T >
T convert_arrayF2C(const T& src, const ulong nx, const ulong ny, const ulong nz) {
	if(src->size() != nx * ny * nz)
		return src;

	// convert to C order
	T dst = BS_KERNEL.create_object(src->bs_resolve_type());
	dst->resize(nx * ny * nz);
	for (ulong i = 0; i < nx; i++)
		for (ulong j = 0; j < ny; j++)
			for (ulong k = 0; k < nz; k++)
				dst->ss(k + j * nz + i * ny * nz) = src->ss(i + j * nx + k * nx * ny);
	return dst;
}

} // eof hidden namespace

void read_surface(const std::string& fname, ulong nx, ulong ny,
	smart_ptr< h5_pool_iface > pool, const std::string& surf_name = "")
{
	// set C locale for proper numbers reading
	setlocale(LC_NUMERIC, "C");

	// open file
	std::ifstream f(fname.c_str(), std::ios::in);
	if(!f) {
		BSERR << "Error opening file " << fname << bs_end;
		return;
	}

	std::string linebuf;
	std::istringstream line_s;
	std::size_t pos;
	t_float point[3];
	ulong row, col;

	// setup buffer with dimensions passed
	ulong dims[3] = {nx, ny, 1};
	spv_float databuf = BS_KERNEL.create_object(v_float::bs_type());
	databuf->resize(dims[0] * dims[1] * 3);
	std::fill(databuf->begin(), databuf->end(), 0.);

	// actually read file
	ulong n_points = 0;
	while(std::getline(f, linebuf)) {
		line_s.clear();

		// parse some info from comments
		if(linebuf[0] == '#') {
			pos = linebuf.find("Grid_size:");
			if(pos != std::string::npos) {
				line_s.str(linebuf.substr(pos + 11));
				// format is "252 x 151", so read intermediate 'x' char
				if(line_s >> nx >> linebuf[0] >> ny) {
					dims[0] = nx; dims[1] = ny;
					// update data buffer size
					databuf->resize(dims[0] * dims[1] * 3);
					std::fill(databuf->begin(), databuf->end(), 0.);
			}
			}
			continue;
		}

		// if we are here, then start array reading
		line_s.str(linebuf);
		if(line_s >> point[0] >> point[1] >> point[2] >> col >> row) {
			const ulong offs = ((row - 1)*ny + col - 1) * 3;
			if(offs > databuf->size() - 3)
				continue;
			std::copy(&point[0], &point[0] + 3, databuf->begin() + offs);
			// count points actually read from file
			++n_points;
		}
	}

	// restore locale
	setlocale(LC_NUMERIC, "");

	// check if we read nothing
	if(!n_points) {
		BSERR << "No valid data is read from " << fname << bs_end;
		return;
	}
	else {
		BSOUT << "read_surface: succefully read " << n_points << " points from " << fname << bs_end;
	}

	// write array to pool
	std::string surf_name_ = surf_name;
	if(!surf_name.size()) {
		surf_name_ = fname.substr(0, fname.rfind('.'));
		const std::size_t slashpos = std::min(surf_name_.rfind('/'), surf_name_.rfind('\\'));
		if(slashpos != std::string::npos)
			surf_name_ = surf_name_.substr(slashpos + 1);
	}

	//npy_intp h5p_dims[] = { 0, npy_intp(dims[0]), 0, npy_intp(dims[1]), 0, 1 };
	npy_intp h5p_dims[] = { npy_intp(dims[0]), npy_intp(dims[1]) };
	//pool->declare_fp_data(surf_name_, 0, 3, &h5p_dims[0], 1);
	pool->declare_fp_data(surf_name_, 0, 2, &h5p_dims[0], 0);
	pool->set_fp_data(surf_name_, databuf);
}

} /* blue_sky */

