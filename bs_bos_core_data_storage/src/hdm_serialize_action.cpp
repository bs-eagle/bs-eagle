/// @file hdm_serialize_action.cpp
/// @brief Implementation of hdm save & load functions called by clients
/// @author uentity
/// @version 
/// @date 30.01.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "bs_bos_core_data_storage_stdafx.h"

#include "hdm_serialize.h"
#include "bs_misc.h"
#include <fstream>
#include <sstream>

// save to and load from text archives
#include <boost/archive/polymorphic_text_iarchive.hpp>
#include <boost/archive/polymorphic_text_oarchive.hpp>
#include <boost/python.hpp>

// path separator
#ifdef UNIX
#define PATHSEP '/'
#else
#define PATHSEP '\\'
#endif

// extension of file with hdm dump
#define HDM_DUMP_EXT ".hdm"

namespace boarch = boost::archive;
namespace boser = boost::serialization;
namespace bp = boost::python;

namespace blue_sky {

// helper that allows to manually set project path and name for serialization process
// invoke with empty strings to clear kernel table
// returns prev content of kernel table
// pass 3 empty strs to clear kernel table
std::vector< std::string > hdm_serialize_set_project_prop(
	const std::string& prj_path = "",
	const std::string& prj_name = "",
	const std::string& deep_copy_suffix = "",
	bool force = false
) {
	// save project path for serialization code
	kernel::idx_dt_ptr p_dt = BS_KERNEL.pert_idx_dt(hdm::bs_type());
	// save existing contents of kernel table
	ulong sz = p_dt->size< std::string >();
	std::vector< std::string > prev_tbl(sz);
	for(ulong i = 0; i < sz; ++i)
		prev_tbl[i] = p_dt->ss< std::string >(i);

	// clear kernel table to not confuse with further saves
	p_dt->clear< std::string >();

	// remember project path
	if(force || prj_path.size() > 0 || prj_name.size() > 0)
		p_dt->insert< std::string >(prj_path);
	// and project name
	if(force || prj_name.size() > 0 || deep_copy_suffix.size() > 0)
		p_dt->insert< std::string >(prj_name);
	// and deep copy suffix
	if(force || deep_copy_suffix.size() > 0) {
		p_dt->insert< std::string >(deep_copy_suffix);
	}

	return prev_tbl;
}

void hdm_serialize_restore_project_prop(const std::vector< std::string >& prev_tbl) {
	kernel::idx_dt_ptr p_dt = BS_KERNEL.pert_idx_dt(hdm::bs_type());
	p_dt->clear< std::string >();
	for(ulong i = 0; i < prev_tbl.size(); ++i)
		p_dt->insert< std::string >(prev_tbl[i]);
}

// hidden implementation details
namespace  {

template< class dst_stream >
dst_stream& hdm_serialize_save_impl(
	dst_stream& f,
	smart_ptr< hdm > t,
	const std::string& prj_path,
	const std::string& prj_name,
	const std::string& deep_copy_suffix
){
	// setup paths
	std::vector< std::string > prev_tbl = hdm_serialize_set_project_prop(
		prj_path, prj_name, deep_copy_suffix
	);

	// make archive
	try {
		boarch::polymorphic_text_oarchive oa(f);
		oa << t;
	}
	catch(const bs_exception& e) {
		BSERR << std::string("Error during hdm_serialize_save: ") + e.what() << bs_end;
	}
	catch(const std::exception& e) {
		BSERR << std::string("Error during hdm_serialize_save: ") + e.what() << bs_end;
	}
	catch(...) {
		BSERR << "Unknown error during hdm_serialize_save!" << bs_end;
		// clear kernel table to not confuse with further saves
		hdm_serialize_restore_project_prop(prev_tbl);
		throw;
	}

	// restore kernel table contents
	hdm_serialize_restore_project_prop(prev_tbl);
	return f;
}

template< class src_stream >
smart_ptr< hdm > hdm_serialize_load_impl(
	src_stream& f,
	const std::string& prj_path,
	const std::string& prj_name
){
	// setup paths
	std::vector< std::string > prev_tbl = hdm_serialize_set_project_prop(
		prj_path, prj_name
	);

	// load archive
	// actually need to rethrow all exceptions
	smart_ptr< hdm > t;
	try {
		boarch::polymorphic_text_iarchive ia(f);
		ia >> t;
	}
	catch(const bs_exception& e) {
		BSERR << std::string("Error during hdm_serialize_load: ") + e.what() << bs_end;
		// clear kernel table
		hdm_serialize_restore_project_prop(prev_tbl);
		throw;
		//t = BS_KERNEL.create_object(hdm::bs_type());
		//t.release();
	}
	catch(const std::exception& e) {
		BSERR << std::string("Error during hdm_serialize_load: ") + e.what() << bs_end;
		// clear kernel table
		hdm_serialize_restore_project_prop(prev_tbl);
		throw;
		//t = BS_KERNEL.create_object(hdm::bs_type());
		//t.release();
	}
	catch(...) {
		BSERR << "Unknown error during hdm_serialize_save!" << bs_end;
		// clear kernel table
		hdm_serialize_restore_project_prop(prev_tbl);
		throw;
	}

	// restore kernel table contents
	hdm_serialize_restore_project_prop(prev_tbl);
	return t;
}

} /* hidden namespace */

/*-----------------------------------------------------------------
 * save hdm
 *----------------------------------------------------------------*/
void hdm_serialize_save(
	smart_ptr< hdm > t,
	const std::string& prj_path,
	const std::string& prj_name,
	const std::string& deep_copy_suffix
){
//#ifdef UNIX
//	const char* enc_name = "ru_RU.UTF-8";
//#else
//	const char* enc_name = "ru_RU.CP1251";
//#endif
//	std::ofstream f(
//		(wstr2str(prj_path, enc_name) + PATHSEP +
//		wstr2str(prj_name, enc_name) + HDM_DUMP_EXT).c_str()
//	);

	std::ofstream f((prj_path + PATHSEP + prj_name + HDM_DUMP_EXT).c_str());
	hdm_serialize_save_impl(f, t, prj_path, prj_name, deep_copy_suffix);
}

std::string hdm_serialize_to_str(
	smart_ptr< hdm > t,
	const std::string& prj_path,
	const std::string& prj_name,
	const std::string& deep_copy_suffix
){
	std::ostringstream f;
	return hdm_serialize_save_impl(f, t, prj_path, prj_name, deep_copy_suffix).str();
}

/*-----------------------------------------------------------------
 * load hdm
 *----------------------------------------------------------------*/
smart_ptr< hdm > hdm_serialize_load(
	const std::string& prj_path,
	const std::string& prj_name
){
//#ifdef UNIX
//	const char* enc_name = "ru_RU.UTF-8";
//#else
//	const char* enc_name = "ru_RU.CP1251";
//#endif
//	std::ifstream f(
//		(wstr2str(prj_path, enc_name) + PATHSEP +
//		wstr2str(prj_name, enc_name) + HDM_DUMP_EXT).c_str()
//	);
	std::ifstream f((prj_path + PATHSEP + prj_name + HDM_DUMP_EXT).c_str());
	return hdm_serialize_load_impl(f, prj_path, prj_name);
}

smart_ptr< hdm > hdm_serialize_from_str(
	const std::string& hdm_dump,
	const std::string& prj_path,
	const std::string& prj_name
){
	std::istringstream f(hdm_dump);
	return hdm_serialize_load_impl(f, prj_path, prj_name);
}

/*-----------------------------------------------------------------
 * export save & load to python
 *----------------------------------------------------------------*/
BOOST_PYTHON_FUNCTION_OVERLOADS(hdm_serialize_save_overl, hdm_serialize_save, 3, 4)
//BOOST_PYTHON_FUNCTION_OVERLOADS(hdm_serialize_load_overl, hdm_serialize_load, 2, 3)
BOOST_PYTHON_FUNCTION_OVERLOADS(hdm_serialize_to_str_overl, hdm_serialize_to_str, 3, 4)
//BOOST_PYTHON_FUNCTION_OVERLOADS(hdm_serialize_from_str_overl, hdm_serialize_from_str, 3, 4)

BOOST_PYTHON_FUNCTION_OVERLOADS(hdm_serialize_setpprop_overl,
	hdm_serialize_set_project_prop, 0, 4
)

namespace python {

BS_API_PLUGIN void py_export_hdm_serialize() {
	bp::def("hdm_serialize_save", &hdm_serialize_save, hdm_serialize_save_overl());
	bp::def("hdm_serialize_load", &hdm_serialize_load);
    // register pvt to/from str serialization
	bp::def("serialize_to_str", &blue_sky::hdm_serialize_to_str, hdm_serialize_to_str_overl());
	bp::def("serialize_from_str", &blue_sky::hdm_serialize_from_str);

	// export only wide-string version
	bp::def(
		"hdm_serialize_set_project_prop", &hdm_serialize_set_project_prop,
		hdm_serialize_setpprop_overl()
	);
	bp::def(
		"hdm_serialize_restore_project_prop", &blue_sky::hdm_serialize_restore_project_prop
	);
}

} /* python */

} /* blue_sky */
