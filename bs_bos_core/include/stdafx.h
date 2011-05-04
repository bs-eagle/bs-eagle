/**
 *       \file  stdafx.h
 *      \brief  Precompiled header... 
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

/**
 * \file stdafx.h
 * \brief precompiled header
 * \author Sergey Miryanov
 * */
#ifndef BS_PRECOMPILED_HEADERS_H_
#define BS_PRECOMPILED_HEADERS_H_

#ifndef UNIX
#include <windows.h>
#endif

#include <vector>
#include <map>
#include <string>
#include <cstdio>
#include <list>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <iterator>
#include <sstream>
#include <cmath>

#include <memory.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

//#include <boost/thread/condition.hpp>
//#include <boost/thread/mutex.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
//#include <boost/assert.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/array.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/tuple/to_seq.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/tuple/elem.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/date_time.hpp>
#include <boost/noncopyable.hpp>
#include <boost/type_traits.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/static_assert.hpp>
#include <boost/format.hpp>

// Handles broken standard libraries better than <iterator>
#include <boost/detail/iterator.hpp>
#include <boost/throw_exception.hpp>
// FIXES for broken compilers
#include <boost/config.hpp>

#ifdef _OPENMP
  #include <omp.h>
#endif // #ifdef _OPENMP

#pragma intrinsic (memset, memcpy)

#include "bs_common.h"
#include "smart_ptr.h"
#include "bs_kernel.h"
#include "bs_link.h"
#include "bs_object_base.h"
#include "bs_tree.h"
#include "bs_exception.h"
#include "bs_prop_base.h"

#ifdef BSPY_EXPORTING_PLUGIN
  #include <boost/python.hpp>
  #include <boost/python/module.hpp>
  #include <boost/python/class.hpp>
  #include <boost/python/def.hpp>
  #include <boost/python/manage_new_object.hpp>
  #include <boost/python/suite/indexing/vector_indexing_suite.hpp>
  #include <boost/python/suite/indexing/map_indexing_suite.hpp>
  #include <boost/python/wrapper.hpp>
  #include <boost/python/iterator.hpp>

  #include "bs_plugin_common.h"
  #include "py_bs_object_base.h"
  #include "py_bs_command.h"
  #include "py_bs_tree.h"
#endif

#include BS_FORCE_PLUGIN_IMPORT ()
//#include "force_inline.h"

#include "bs_assert.h"
#include "bos_report.h"
#include "auto_value.h"
#include "throw_exception.h"

#include "memory_macroses.h"
#include "macroses.h"
//#include "lu_decomposition.h"
#include "debug_macroses.h"
#include "matrix_macroses.h"

#include "omp_tools.h"
#include "timers.h"
#include "mv_tools.h"

#include "aligned_allocator.h"
#include "shared_vector.h"

#ifdef _MPI
  #include "mpi_type_t.h"
  #include "mpi_vector.h"
#endif  // #ifdef _MPI

#include "strategies.h"

#include "bos_map.h"
#include "constants.h"
#include "main_def.h"
#include "err_num_def.h"
#include "vector_assign.h"
//#include "save_seq_vector.h"
//#include "jacobian_matrix.h"

#ifdef _HDF5
//#include "bs_hdf5_storage.h"
#endif

#include "property_base.h"
#include "named_pbase_access.h"

//#include "linear_solvers.h"
#include "print_process_memory_info.h"
#include "write_time_to_log.h"

//#include "bs_array_map.h"
#include "interpolation_macro.h"
#include "rs_mesh_iface.h"
#include "rs_smesh_iface.h"
#include "localization.h"
#include "arrays.h"
#include "arrays_tables.h"
#include "convert_units.h"

// sergey.miryanov - to compile pvt
#include "table_iface.h"

#include "pvt_water.h"
#include "pvt_gas.h"
#include "pvt_oil.h"

// sergey.miryanov - to compile scal
#include "vartype_table_iface.h"

#include "scal_3p_iface.hpp"
//#include "scale_array_holder.h"
//#include "scal_region_info.h"
//#include "scal_region.h"
//#include "scal_2p_data_holder.h"
#include "jfunction.h"

#include "rock_grid.h"

#ifdef BSPY_EXPROTING_PLUGIN
#include "py_matrix_base.h"
#include "py_bcsr_matrix.h"
#include "py_linear_solvers.h"
#include "py_jacobian_matrix.h"
#include "py_scal_wrapper.h"
#include "py_data_class.h"
#endif

#include BS_STOP_PLUGIN_IMPORT ()

#endif  // #ifndef BS_PRECOMPILED_HEADERS_H_
