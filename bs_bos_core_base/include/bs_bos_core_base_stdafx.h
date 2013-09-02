/**
 * \file stdafx.h
 * \brief precompiled header
 * \author Sergey Miryanov
 * */
#ifndef BS_BOS_CORE_BASE_PRECOMPILED_HEADERS_H_
#define BS_BOS_CORE_BASE_PRECOMPILED_HEADERS_H_

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

#include <vector>

#include <boost/array.hpp>

#include "bs_common.h"
#include "smart_ptr.h"
#include "bs_kernel.h"
#include "bs_link.h"
#include "bs_object_base.h"
#include "bs_tree.h"
#include "bs_exception.h"
#include "bs_assert.h"

#ifdef BSPY_EXPORTING_PLUGIN
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

#include "force_inline.h"
#include "bos_report.h"


#include "aligned_allocator.h"

#endif // #ifndef BS_BOS_CORE_BASE_PRECOMPILED_HEADERS_H_
