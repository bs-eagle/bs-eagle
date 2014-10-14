/// @file wlog_tools.cpp
/// @brief Well Logs manipulation algorithms
/// @author uentity
/// @version 1.0
/// @date 27.02.2014
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#ifndef WLOG_TOOLS_TMEEHO84
#define WLOG_TOOLS_TMEEHO84

#include "bs_kernel.h"
#include "conf.h"
#include <limits>
// DEBUG
//#include <iostream>

namespace blue_sky {

// hidden implementation
namespace {

typedef v_float::iterator vf_iterator;
typedef v_float::const_iterator cvf_iterator;
typedef v_float::reverse_iterator vf_riterator;
typedef v_float::const_reverse_iterator cvf_riterator;

// depth and share same iterator type
template< class data_iterator, class grid_iterator, class res_iterator >
void projection_impl(
	//const spv_float& wlog_data, const spv_float& wlog_dept,
	const data_iterator& p_data_begin, const data_iterator& p_data_end,
	const data_iterator& p_dept_begin, const data_iterator& p_dept_end,
	grid_iterator& p_grid, const grid_iterator& p_grid_end,
	res_iterator& p_res
) {
	// setup iterators
	//const cvf_iterator p_dept_begin = wlog_dept->begin();
	//const cvf_iterator p_dept_end = p_dept_begin + wlog_sz;

	const ulong wlog_sz = std::min< ulong >(p_data_end - p_data_begin, p_dept_end - p_dept_begin);
	//const ulong wlog_sz = std::min(wlog_data->size(), wlog_dept->size());
	if(wlog_sz == 0)
		return;

	// find starting point on well log
	data_iterator p_dept = std::lower_bound(
		p_dept_begin, p_dept_end, *p_grid
	);
	if(p_dept == p_dept_end)
		// grid is fully outside of well log depths
		return;

	data_iterator p_data = p_data_begin + (p_dept - p_dept_begin);

	// DEBUG
	//std::cout << "++wlog_mean_projection:" << std::endl;
	//std::cout << "log_dept = [" << *p_dept << ", " << *(p_dept_end - 1) << "]" << std::endl;
	//std::cout << "grid = [" << *p_grid << ", " << *(p_grid_end - 1) << "]" << std::endl;
	// position dest grid to next boundary
	++p_grid;
	const res_iterator p_res_end = p_res + (p_grid_end - p_grid - 1);

	// main cycle
	t_float win_sum = 0;
	// win_sz counts only valid numbers that fits in window
	// raw_sz counts all values in window, including NaNs
	ulong win_sz = 0, raw_sz = 0;
	while(p_grid != p_grid_end) {
		// if we reached window boundary
		if(*p_dept > *p_grid) {
			if(win_sz)
				*p_res = win_sum / win_sz;
			else if(!raw_sz && (p_dept != p_dept_begin)) {
				// assign prev log value only if we had no NaNs in window
				*p_res = *(p_data - 1);
			}

			// next step on dest grid and resulting array
			if(++p_res == p_res_end)
				break;
			++p_grid;
			win_sum = 0;
			win_sz = 0;
			raw_sz = 0;
		}
		else {
			// we're inside window
			// check for NaN
			if(!(*p_data != *p_data)) {
				win_sum += *p_data;
				++win_sz;
			}
			// count all values in raw_sz
			++raw_sz;
			// next step in well log
			++p_data; ++p_dept;
			if(p_dept == p_dept_end) {
				if(win_sz)
					*p_res = win_sum / win_sz;
				break;
			}
		}
	}
}

} // eof hidden namespace

BS_API_PLUGIN spv_float wlog_mean_projection(
	spv_float wlog_data, spv_float wlog_dept, spv_float dest_grid
) {
	// sanity
	spv_float res = BS_KERNEL.create_object(v_float::bs_type());
	if(!wlog_data->size() || !wlog_dept->size() || dest_grid->size() < 2 || !res)
		return res;
	res->resize(dest_grid->size() - 1);
	// fill resulting array with NaN
	std::fill(res->begin(), res->end(), std::numeric_limits< double >::quiet_NaN());

	// check grid ordering
	if(dest_grid->ss(0) < dest_grid->ss(dest_grid->size() - 1)) {
		// normal ordering
		cvf_iterator p_grid = dest_grid->begin();
		const cvf_iterator p_grid_end = dest_grid->end();
		vf_iterator p_res = res->begin();

		// check depth ordering
		if(wlog_dept->ss(0) < wlog_dept->ss(wlog_dept->size() - 1)) {
			projection_impl(
				wlog_data->begin(), wlog_data->end(),
				wlog_dept->begin(), wlog_dept->end(),
				p_grid, p_grid_end, p_res
			);
		}
		else {
			projection_impl(
				wlog_data->rbegin(), wlog_data->rend(),
				wlog_dept->rbegin(), wlog_dept->rend(),
				p_grid, p_grid_end, p_res
			);
		}
	}
	else {
		// reversed grid ordering
		cvf_riterator p_grid = dest_grid->rbegin();
		const cvf_riterator p_grid_end = dest_grid->rend();
		vf_riterator p_res = res->rbegin();

		// check depth ordering
		if(wlog_dept->ss(0) < wlog_dept->ss(wlog_dept->size() - 1)) {
			projection_impl(
				wlog_data->begin(), wlog_data->end(),
				wlog_dept->begin(), wlog_dept->end(),
				p_grid, p_grid_end, p_res
			);
		}
		else {
			projection_impl(
				wlog_data->rbegin(), wlog_data->rend(),
				wlog_dept->rbegin(), wlog_dept->rend(),
				p_grid, p_grid_end, p_res
			);
		}
		//projection_impl(wlog_data, wlog_dept, p_grid, p_grid_end, p_res);
	}
	return res;
}

} /* namespace blue_sky */

#endif /* end of include guard: WLOG_TOOLS_TMEEHO84 */

