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

namespace blue_sky {

// hidden implementation
namespace {

typedef v_float::iterator vf_iterator;
typedef v_float::const_iterator cvf_iterator;
typedef v_float::reverse_iterator vf_riterator;
typedef v_float::const_reverse_iterator cvf_riterator;

template< class dept_iterator, class res_iterator >
void projection_impl(
	const spv_float& wlog_data, const spv_float& wlog_dept,
	dept_iterator& p_grid, const dept_iterator& p_grid_end,
	res_iterator& p_res
) {
	const ulong wlog_sz = std::min(wlog_data->size(), wlog_dept->size());
	//spv_float res = BS_KERNEL.create_object(v_float::bs_type());
	if(wlog_sz == 0)
		return;
	//res->resize(p_grid_end - p_grid - 1);

	// setup iterators
	const cvf_iterator p_dept_end = wlog_dept->begin() + wlog_sz;
	// find starting point on well log
	cvf_iterator p_dept = std::lower_bound(
		wlog_dept->begin(), wlog_dept->begin() + wlog_sz, *p_grid
	);
	cvf_iterator p_data = wlog_data->begin() + (p_dept - wlog_dept->begin());
	//vf_iterator  p_res  = res->begin();
	// position dest grid to next boundary
	++p_grid;
	// zero-fill resulting array
	// TODO: do something for out-of-bounds values, instead of simple 0
	std::fill(p_res, p_res + (p_grid_end - p_grid), 0.0);

	// main cycle
	t_float win_sum = 0;
	ulong win_sz = 0;
	while(p_grid != p_grid_end) {
		// if we reached window boundary
		if(*p_dept >= *p_grid) {
			if(win_sz)
				*p_res = win_sum / win_sz;
			else
				*p_res = *p_data;
			// next step on dest grid and resulting array
			++p_res;
			++p_grid;
			win_sum = 0;
			win_sz = 0;
		}
		else {
			++p_data; ++p_dept;
			if(p_dept == p_dept_end) {
				// here we always have win_sz >= 1
				*p_res = win_sum / win_sz;
				break;
			}
		}
		// we're inside window
		win_sum += *p_data;
		++win_sz;
	}
}

} // eof hidden namespace

BS_API spv_float wlog_mean_projection(
	spv_float wlog_data, spv_float wlog_dept, spv_float dest_grid
) {
	// sanity
	spv_float res = BS_KERNEL.create_object(v_float::bs_type());
	if(dest_grid->size() < 2 || !res)
		return res;
	res->resize(dest_grid->size() - 1);

	if(dest_grid->ss(0) < dest_grid->ss(dest_grid->size() - 1)) {
		// normal ordering
		cvf_iterator p_grid = dest_grid->begin();
		const cvf_iterator p_grid_end = dest_grid->end();
		vf_iterator p_res = res->begin();
		projection_impl(wlog_data, wlog_dept, p_grid, p_grid_end, p_res);
	}
	else {
		// reversed grid ordering
		cvf_riterator p_grid = dest_grid->rbegin();
		const cvf_riterator p_grid_end = dest_grid->rend();
		vf_riterator p_res = res->rbegin();
		projection_impl(wlog_data, wlog_dept, p_grid, p_grid_end, p_res);
	}
	return res;
}

} /* namespace blue_sky */

#endif /* end of include guard: WLOG_TOOLS_TMEEHO84 */

