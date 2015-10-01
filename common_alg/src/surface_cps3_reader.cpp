/// @file surface_cps3_reader.cpp
/// @brief Read surface in CPS-3 ASCII format
/// @author uentity
/// @version 1.0
/// @date 30.09.2015
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "bs_kernel.h"
#include "conf.h"

#include <fstream>
#include <sstream>
#include <limits>

namespace blue_sky {
using namespace std;

// hidden namespace
namespace {

/*-----------------------------------------------------------------
 * helper structure to get dx[i] when dim is given by one number
 *----------------------------------------------------------------*/
struct dim_subscript {
	dim_subscript(const double dim, const double offset = .0)
		: dim_(dim), offset_(offset), sum_(offset)
	{}

	void reset() { sum_ = offset_; }

	double operator[](t_ulong idx) {
		return offset_ + dim_ * idx;
	}

private:
	const double dim_;
	const double offset_;
	double sum_;
};

} /* eof hidden namespace  */

ulong read_cps3_grid(std::ifstream& f, ulong* dims, spv_float& databuf) {
	std::string linebuf;
	std::istringstream line_s;

	// [x_min, x_max, y_min, y_max, z_min, z_max]
	t_float bounds[] = {0., 0., 0., 0., 0., 0.};
	// increment in X and Y directions
	t_float D[] = {0., 0.};

	// actually read file
	ulong n_points = 0;
	bool data_started = false, overflow = false;
	ulong i = 0, j = 0;
	while(std::getline(f, linebuf)) {
		line_s.clear();

		// parse header
		if(!data_started) {
			size_t pos = linebuf.find("FSLIMI");
			if(pos != string::npos) {
				line_s.str(linebuf.substr(pos + 6));
				if(!(line_s >> bounds[0] >> bounds[1] >> bounds[2] >> bounds[3] >> bounds[4] >> bounds[5])) {
					BSERR << "Error reading limits from CPS-3 FSLIMI keyword" << bs_end;
					return 0;
				}
			}

			pos = linebuf.find("FSNROW");
			if(pos != string::npos) {
				line_s.str(linebuf.substr(pos + 6));
				if(!(line_s >> dims[1] >> dims[0])) {
					BSERR << "Error reading dimensions from CPS-3 FSNROW keyword" << bs_end;
					return 0;
				}
				// resize Z data buffer
				databuf->resize(dims[0] * dims[1] * 3);
				// TODO: enable NaN when GUI part will be fixed to support it
				//fill(databuf->begin(), databuf->end(), std::numeric_limits< double >::quiet_NaN());
				fill(databuf->begin(), databuf->end(), 1e30);
				// data is read from upper left corner, i = 0, j = max
				i = 0;
				j = dims[1] - 1;
			}

			pos = linebuf.find("FSXINC");
			if(pos != string::npos) {
				line_s.str(linebuf.substr(pos + 6));
				if(!(line_s >> D[1] >> D[0])) {
					BSERR << "Error reading deltas from CPS-3 FSXINC keyword" << bs_end;
					return 0;
				}
			}

			// check if we reached actual data
			if(linebuf.find("->MSMODL") != string::npos) {
				data_started = true;
			}
			continue;
		}

		// if we are here, then we're in data section
		line_s.str(linebuf);
		t_float val;
		while(line_s >> val) {
			// store x, y, z surface value
			const v_float::iterator dest = databuf->begin() + (dims[1] * i + j) * 3;
			if(dest >= databuf->end() - 3) {
				overflow = true;
				break;
			}
			dest[0] = bounds[0] + D[0] * i;
			dest[1] = bounds[2] + D[1] * j;
			if(val >= bounds[4] && val <= bounds[5])
				dest[2] = val;

			// step to next surface point
			if(--j > dims[1]) {
				++i;
				j = dims[1] - 1;
			}
			++n_points;
		}
	}

	return n_points;
}

} /* namespace blue_sky */
