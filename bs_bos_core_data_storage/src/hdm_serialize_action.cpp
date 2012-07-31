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

void hdm_serialize_save(
	smart_ptr< hdm > t,
	const std::string& prj_path,
	const std::string& prj_name,
	bool deep_copy
){
	// save project path for serialization code
	kernel::idx_dt_ptr p_dt = BS_KERNEL.pert_idx_dt(hdm::bs_type());
	//std::string well_pool_filename = prj_name + "_well_pool.db";
	p_dt->insert< std::string >(prj_path);
	// and project name
	p_dt->insert< std::string >(prj_name);
	// inform to explicitly copy storage files
	if(deep_copy)
		p_dt->insert< std::string >("deep_copy");

	// make archive
	std::string fname = prj_path + PATHSEP + prj_name + HDM_DUMP_EXT;
	std::ofstream f(fname.c_str());
	boarch::polymorphic_text_oarchive oa(f);
	oa << t;

	// clear kernel table to not confuse with further saves
	p_dt->clear< std::string >();
}

smart_ptr< hdm > hdm_serialize_load(
	const std::string& prj_path,
	const std::string& prj_name,
	bool deep_copy
){
	// save project path for serialization code
	kernel::idx_dt_ptr p_dt = BS_KERNEL.pert_idx_dt(hdm::bs_type());
	//std::string well_pool_filename = prj_name + "_well_pool.db";
	p_dt->insert< std::string >(prj_path);
	// and project name
	p_dt->insert< std::string >(prj_name);
	// inform to explicitly load storage files
	if(deep_copy)
		p_dt->insert< std::string >("deep_copy");

	// load archive
	std::ifstream f((prj_path + PATHSEP + prj_name + HDM_DUMP_EXT).c_str());
	boarch::polymorphic_text_iarchive ia(f);
	smart_ptr< hdm > t;
	ia >> t;

	// clear kernel table
	p_dt->clear< std::string >();
	// clear dangling refs of loaded in kernel
	BS_KERNEL.tree_gc();

	return t;
}

/*-----------------------------------------------------------------
 * export save & load to python
 *----------------------------------------------------------------*/
BOOST_PYTHON_FUNCTION_OVERLOADS(hdm_serialize_save_overl, hdm_serialize_save, 3, 4)
BOOST_PYTHON_FUNCTION_OVERLOADS(hdm_serialize_load_overl, hdm_serialize_load, 2, 3)

namespace python {

BS_API_PLUGIN void py_export_hdm_serialize() {
	bp::def("hdm_serialize_save", &hdm_serialize_save, hdm_serialize_save_overl());
	bp::def("hdm_serialize_load", &hdm_serialize_load, hdm_serialize_load_overl());
    // register pvt to/from str serialization
	bp::def("serialize_to_str", &blue_sky::serialize_to_str< hdm >);
	bp::def("serialize_from_str", &blue_sky::serialize_from_str< hdm >);
}

} /* python */

} /* blue_sky */
