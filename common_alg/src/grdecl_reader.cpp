/// @file grdecl_reader.cpp
/// @brief Read GRDECL Petrel text files into set of BlueSky arrays
/// @author uentity
/// @version 1.0
/// @date 03.04.2013
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "bs_kernel.h"
#include "conf.h"
#include "bos_reader_iface.h"
#include "pool_iface.h"
#include "bs_misc.h"

#include <sstream>
#include <locale>

#define MAX_LINE_LEN 65536

namespace blue_sky {
typedef smart_ptr< bos_reader_iface > sp_bos_reader;
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

template< class T >
inline bool read_array(const sp_bos_reader& br, const char* array_name, const T& dst, const ulong sz) {
	struct array_reader {
		static ulong go(const sp_bos_reader& br, const char* array_name, const spv_float& dst, const ulong sz) {
			return (ulong)br->read_fp_array(array_name, &dst->ss(0), sz);
		}
		static ulong go(const sp_bos_reader& br, const char* array_name, const spv_int& dst, const ulong sz) {
			return (ulong)br->read_int_array(array_name, &dst->ss(0), sz);
		}
	};

	if(!sz) return true;
	dst->resize(sz);
	array_reader AR;
	ulong N;
	if((N = AR.go(br, array_name, dst, sz)) != sz) {
		BSERR << "Error reading " << array_name << " array, array size " << sz <<
			", actually read " << N << " elements" << bs_end;
		return false;
	}
	return true;
}

}

/// @brief Parse GRDECL file and read all contained arrays
///
/// It's assumed that all arrays are of size nx*ny*nz, besides
/// COORD which is (nx + 1)*(ny + 1)*(nz + 1)*2 and
/// ZCORN which is nx*ny*nz*8
/// @param fname - input filename
/// @param nx - number of cells in X dimension
/// @param ny - number of cells in Y dimension
/// @param nz - number of cells in Z dimension
///
/// @return vector of arrays read from file
void read_grdecl(const std::string& fname, const std::string dir, ulong nx, ulong ny, ulong nz,
	smart_ptr< h5_pool_iface > pool)
{
	// precalc sizes
	char order = 'F';
	const char* coord_name = "COORD";
	const char* zcorn_name = "ZCORN";
	const char* specgrid_name = "SPECGRID";
	t_long dims[3] = { t_long(nx), t_long(ny), t_long(nz) };

	// open file using bos_reader
	smart_ptr< bos_reader_iface > br = BS_KERNEL.create_object("bos_reader");
	if(!br) return;
	if(br->open(fname.c_str(), dir.c_str()) < 0) {
		BSERR << "Error opening file " << fname << bs_end;
		return;
	}

	// set C locale for proper numbers reading
	setlocale(LC_NUMERIC, "C");

	// start reading
	std::string buf(MAX_LINE_LEN, ' ');
	std::istringstream is;
	spv_float databuf = BS_KERNEL.create_object(v_float::bs_type());
	spv_int idatabuf = BS_KERNEL.create_object(v_int::bs_type());
	bool is_mesh_read = false, is_pool_touched = false;
	long N;
	while((N = br->read_line(const_cast< char* >(buf.c_str()), MAX_LINE_LEN)) > 0) {
		// cut trailing garbage from buf
		const std::string line = buf.substr(0, N);
		// skip /
		if(line[0] == '/' || line[N - 1] == '/')
			continue;

		if(line == specgrid_name) {
			// SPECGRID keyword can override dimensions
			N = br->read_line(const_cast< char* >(buf.c_str()), MAX_LINE_LEN);
			is.str(buf.substr(0, N));
			is.clear();

			// read dims from string
			if(is >> dims[0] >> dims[1] >> dims[2]) {
				nx = ulong(dims[0]); ny = ulong(dims[1]); nz = ulong(dims[2]);
				// reset pool dims
				pool->set_pool_dims(&dims[0], 3);
				is_pool_touched = true;
			}
			else
				continue;
			// read array order
			// skip 4-d number
			if(!(is >> dims[0])) continue;
			char new_order;
			if(is >> new_order)
				order = new_order;
		}
		else if(line == coord_name) {
			// process COORD
			if(read_array(br, coord_name, databuf, (nx + 1)*(ny + 1)*6)) {
				if(!is_pool_touched) {
					pool->set_pool_dims(&dims[0], 3);
					is_pool_touched = true;
				}
				npy_intp coord_dims[6] = { 1, 1, 1, 1, 0, 6 };
				pool->declare_fp_data(coord_name, 0, 3, &coord_dims[0], 1);
				pool->set_fp_data(coord_name, databuf);
			}
		}
		else if(line == zcorn_name) {
			// process ZCORN
			if(read_array(br, zcorn_name, databuf, nx * ny * nz * 8)) {
				if(!is_pool_touched) {
					pool->set_pool_dims(&dims[0], 3);
					is_pool_touched = true;
				}
				npy_intp zcorn_dims[6] = { 2, 0, 2, 0, 2, 0 };
				pool->declare_fp_data(zcorn_name, 0, 3, &zcorn_dims[0], 1);
				pool->set_fp_data(zcorn_name, databuf);
			}
			is_mesh_read = true;
			continue;
		}

		// skip all lines until mesh is read
		// TODO: replace it with better solution
		// take a parameter-list of known properties to read
		if(!is_mesh_read) continue;

		// if we are here then we should process unknown properties array of size nx * ny * nz
		is.str(line);
		is.clear();
		// read array name
		std::string prop_name;
		is >> prop_name;

		// read prop data
		bool read_res = false;
		if(prop_name == "ACTNUM")
			read_res = read_array(br, prop_name.c_str(), idatabuf, nx * ny * nz);
		else
			read_res = read_array(br, prop_name.c_str(), databuf, nx * ny * nz);
		if(!read_res)
			continue;
		// set dims once
		if(!is_pool_touched) {
			pool->set_pool_dims(&dims[0], 3);
			is_pool_touched = true;
		}

		// save to pool
		npy_intp prop_dims[6] = { 1, 0, 1, 0, 1, 0 };
		if(prop_name == "ACTNUM") {
			// store only C-ordered arrays
			if(order == 'F')
				idatabuf = convert_arrayF2C(idatabuf, nx, ny, nz);
			pool->declare_i_data(prop_name, 0, 3, &prop_dims[0], 1);
			pool->set_i_data(prop_name, idatabuf);
		}
		else {
			// store only C-ordered arrays
			if(order == 'F')
				databuf = convert_arrayF2C(databuf, nx, ny, nz);
			pool->declare_fp_data(prop_name, 0, 3, &prop_dims[0], 1);
			pool->set_fp_data(prop_name, databuf);
		}
	}

	// restore locale
	setlocale(LC_NUMERIC, "");
	br->close();
}

// wrapper for wide strings
//void read_grdecl_w(const std::wstring& fname, const std::wstring dir, ulong nx, ulong ny, ulong nz,
//	smart_ptr< h5_pool_iface > pool) {
//
//	// decode wide strings and call 
//	read_grdecl(
//#ifdef UNIX
//		wstr2str(fname)
//#else
//		wstr2str(fname, "ru_RU.CP1251")
//#endif
//		,
//#ifdef UNIX
//		wstr2str(dir)
//#else
//		wstr2str(dir, "ru_RU.CP1251")
//#endif
//		, nx, ny, nz, pool
//	);
//}

} /* blue_sky */

