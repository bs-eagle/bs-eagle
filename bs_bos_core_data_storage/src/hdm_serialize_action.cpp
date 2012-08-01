/// @file hdm_serialize_action.cpp
/// @brief Implementation of hdm save & load functions called by clients
/// @author uentity
/// @version 
/// @date 30.01.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "bs_bos_core_data_storage_stdafx.h"

#include "hdm_serialize.h"
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

namespace  {

template< class dst_stream >
dst_stream& hdm_serialize_save_impl(
	dst_stream& f,
	smart_ptr< hdm > t,
	const std::string& prj_path,
	const std::string& prj_name,
	const std::string& deep_copy_suffix
){
	// save project path for serialization code
	kernel::idx_dt_ptr p_dt = BS_KERNEL.pert_idx_dt(hdm::bs_type());
	//std::string well_pool_filename = prj_name + "_well_pool.db";
	p_dt->insert< std::string >(prj_path);
	// and project name
	p_dt->insert< std::string >(prj_name);
	// inform to explicitly copy storage files
	if(deep_copy_suffix.size())
		p_dt->insert< std::string >(deep_copy_suffix);

	// make archive
	boarch::polymorphic_text_oarchive oa(f);
	oa << t;

	// clear kernel table to not confuse with further saves
	p_dt->clear< std::string >();
	return f;
}

template< class src_stream >
smart_ptr< hdm > hdm_serialize_load_impl(
	src_stream& f,
	const std::string& prj_path,
	const std::string& prj_name
){
	// save project path for serialization code
	kernel::idx_dt_ptr p_dt = BS_KERNEL.pert_idx_dt(hdm::bs_type());
	//std::string well_pool_filename = prj_name + "_well_pool.db";
	p_dt->insert< std::string >(prj_path);
	// and project name
	p_dt->insert< std::string >(prj_name);

	// load archive
	boarch::polymorphic_text_iarchive ia(f);
	smart_ptr< hdm > t;
	ia >> t;

	// clear kernel table
	p_dt->clear< std::string >();
	// clear dangling refs of loaded in kernel
	BS_KERNEL.tree_gc();

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
	std::string fname = prj_path + PATHSEP + prj_name + HDM_DUMP_EXT;
	std::ofstream f(fname.c_str());
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

namespace python {

BS_API_PLUGIN void py_export_hdm_serialize() {
	bp::def("hdm_serialize_save", &hdm_serialize_save, hdm_serialize_save_overl());
	bp::def("hdm_serialize_load", &hdm_serialize_load);
    // register pvt to/from str serialization
	bp::def("serialize_to_str", &blue_sky::hdm_serialize_to_str, hdm_serialize_to_str_overl());
	bp::def("serialize_from_str", &blue_sky::hdm_serialize_from_str);
}

} /* python */

} /* blue_sky */
