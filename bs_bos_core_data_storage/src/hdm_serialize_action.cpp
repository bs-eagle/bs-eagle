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
	const std::string& prj_name
){
	// save project path for serialization code
	kernel::idx_dt_ptr p_dt = BS_KERNEL.pert_idx_dt(hdm::bs_type());
	//std::string well_pool_filename = prj_name + "_well_pool.db";
	p_dt->insert< std::string >(prj_path);
	// and project name
	p_dt->insert< std::string >(prj_name);

	// make archive
	std::string fname = prj_path + PATHSEP + prj_name + HDM_DUMP_EXT;
	std::ofstream f(fname.c_str());
	boarch::polymorphic_text_oarchive oa(f);
	oa << t;
}

smart_ptr< hdm > hdm_serialize_load(
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
	std::ifstream f((prj_path + PATHSEP + prj_name + HDM_DUMP_EXT).c_str());
	boarch::polymorphic_text_iarchive ia(f);
	smart_ptr< hdm > t;
	ia >> t;
	return t;
}

/*-----------------------------------------------------------------
 * export save & load to python
 *----------------------------------------------------------------*/
namespace python {

BS_API_PLUGIN void py_export_hdm_serialize() {
	bp::def("hdm_serialize_save", &hdm_serialize_save);
	bp::def("hdm_serialize_load", &hdm_serialize_load);
    // register pvt to/from str serialization
	bp::def("serialize_to_str", &blue_sky::serialize_to_str< hdm >);
	bp::def("serialize_from_str", &blue_sky::serialize_from_str< hdm >);
}

} /* python */

} /* blue_sky */
