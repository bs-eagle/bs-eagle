/// @file surface_earth_vision_reader.cpp
/// @brief Read surface in EarthVision format
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

ulong read_earth_vision_grid(std::ifstream& f, ulong* dims, spv_float& databuf) {
	std::string linebuf;
	std::istringstream line_s;
	std::size_t pos;
	t_float point[3];
	ulong nx, ny, row, col;

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
					// TODO: enable NaN when GUI will support it
					//std::fill(databuf->begin(), databuf->end(), std::numeric_limits< double >::quiet_NaN());
					std::fill(databuf->begin(), databuf->end(), 1e30);
				}
			}
			continue;
		}

		// if we are here, then start array reading
		line_s.str(linebuf);
		if(line_s >> point[0] >> point[1] >> point[2] >> col >> row) {
			//const ulong offs = ((row - 1)*ny + col - 1) * 3;
			const ulong offs = ((col - 1)*ny + row - 1) * 3;
			if(offs > databuf->size() - 3)
				continue;
			std::copy(&point[0], &point[0] + 3, databuf->begin() + offs);
			// count points actually read from file
			++n_points;
		}
	}

	return n_points;
}

} /* namespace blue_sky */

