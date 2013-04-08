/**
 * \file stdafx.h
 * \brief precompiled header
 * \author Sergey Miryanov
 * */
#ifndef BS_BOS_CORE_DATA_STORAGE_PRECOMPILED_HEADERS_H_
#define BS_BOS_CORE_DATA_STORAGE_PRECOMPILED_HEADERS_H_

// Modify the following defines if you have to target a platform prior to the ones specified below.
// Refer to MSDN for the latest info on corresponding values for different platforms.
#ifndef WINVER				// Allow use of features specific to Windows XP or later.
#define WINVER 0x0501		// Change this to the appropriate value to target other versions of Windows.
#endif

#ifndef _WIN32_WINNT		// Allow use of features specific to Windows XP or later.                   
#define _WIN32_WINNT 0x0501	// Change this to the appropriate value to target other versions of Windows.
#endif						

#ifndef _WIN32_WINDOWS		// Allow use of features specific to Windows 98 or later.
#define _WIN32_WINDOWS 0x0410 // Change this to the appropriate value to target Windows Me or later.
#endif

#ifndef _WIN32_IE			// Allow use of features specific to IE 6.0 or later.
#define _WIN32_IE 0x0600	// Change this to the appropriate value to target other versions of IE.
#endif

#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers

#ifdef BSPY_EXPORTING_PLUGIN
#include <boost/python.hpp>
#endif

// Windows Header Files:
#ifndef UNIX
#include <Windows.h>
#else
#include <pthread.h>
#endif

#include <sstream>
#include <math.h>
#include <stack>
#include <set>
#include <map>
#include <list>

#include <boost/array.hpp>
#include <boost/format.hpp>

#include "bs_common.h"
#include "bs_prop_base.h"

#include "aligned_allocator.h"
#include "shared_vector.h"
#include "bs_assert.h"

#ifdef BSPY_EXPORTING_PLUGIN
#include <boost/python.hpp>
#endif

#include BS_FORCE_PLUGIN_IMPORT ()
#include "bos_report.h"
#include "force_inline.h"
#include "auto_value.h"
#include "err_num_def.h"

#include "vector_assign.h"
#include "locale_keeper.h"
#include "interpolation_macro.h"
//#include "naive_file_reader.h"
#include BS_STOP_PLUGIN_IMPORT ()

#endif // #ifndef BS_BOS_CORE_DATA_STORAGE_PRECOMPILED_HEADERS_H_
