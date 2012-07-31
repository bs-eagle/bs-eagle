/// @file sql_well_serialize.cpp
/// @brief Serialization implementations for sql_well
/// @author uentity
/// @version 
/// @date 20.01.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "bs_bos_core_data_storage_stdafx.h"

#include "bs_serialize.h"
#include "sql_well.h"

#include <boost/serialization/string.hpp>

#include <boost/uuid/uuid.hpp>
//#include <boost/uuid/string_generator.hpp>
#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid_io.hpp>

#include <stdio.h>

// path separator
#ifdef UNIX
#define PATHSEP '/'
#else
#define PATHSEP '\\'
#endif

using namespace blue_sky;
namespace bu = boost::uuids;
namespace boser = boost::serialization;

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(save, blue_sky::sql_well)
	// save file_name
	ar << t.file_name;

	// check if database is open
	int cur_s, hiwtr, opres = SQLITE_OK + 1;
	sqlite3* src_db = t.db;
	if(src_db)
		opres = sqlite3_db_status(src_db, 0, &cur_s, &hiwtr, 0);
	if(opres != SQLITE_OK) {
		// probably db is closed or not open or whatever
		// try to open it using contained file_name
		if(t.file_name.size()) {
			opres = sqlite3_open(t.file_name.c_str(), &src_db);
			if(opres != SQLITE_OK) {
				sqlite3_close(src_db);
				src_db = 0;
			}
		}
	}

	// flag indicating whether we should backup db
	bool do_write_db = true;

	// if db not open now, we're done
	if(opres != SQLITE_OK) {
		do_write_db = false;
		ar << do_write_db;
		return;
	}
	ar << do_write_db;

	// check db has an associated filename
	// strangely there is no such API function though it is documented
	//std::string db_fname(sqlite3_db_filename(t.db, "main"));
	std::string db_fname, prj_path, prj_name, db_basename;
	bu::uuid db_uuid = bu::random_generator()();
	// check if hdm serializer saved filename for us
	kernel::idx_dt_ptr kdt = BS_KERNEL.pert_idx_dt(BS_KERNEL.find_type("hdm").td_);
	if(kdt && kdt->size< std::string >()) {
		prj_path = kdt->ss< std::string >(0);
		prj_name = kdt->ss< std::string >(1);
		db_basename = prj_name;
		// if we should make a deep copy of well pool - add a random salt to filename
		if(kdt->size< std::string >() > 2 && kdt->ss< std::string >(2) == "deep_copy")
			db_basename += std::string("_") + bu::to_string(db_uuid);
		db_basename += "_well_pool.db";
		db_fname = prj_path + PATHSEP + db_basename;
	}
	else
		db_basename = t.file_name;

	// generate uuid that is part of filename
	if(db_basename.empty()) {
		db_basename = bu::to_string(db_uuid) + "_well_pool.db";
		db_fname = std::string("file:") + db_basename;
	}

	// save basename first
	ar << db_basename;
	ar << db_fname;

	// open db for backup
	sqlite3* save_db = NULL;
	sqlite3_backup* save_bu = NULL;
	opres = sqlite3_open(db_fname.c_str(), &save_db);
	// and prepare backup process
	if(opres == SQLITE_OK)
		save_bu = sqlite3_backup_init(save_db, "main", src_db, "main");
	if(!save_bu) {
		// probably source & dest db are the same
		sqlite3_close(save_db);
		return;
	}

	// copy data to save_db
	while((opres = sqlite3_backup_step(save_bu, -1)) == SQLITE_OK) {}
	sqlite3_backup_finish(save_bu);
	sqlite3_close(save_db);
	// SQLITE_BUSY returned if we're trying to backup into source db
	if(opres != SQLITE_DONE && opres != SQLITE_BUSY) {
		// if error happens - kill db file
		::remove(db_fname.c_str());
	}

	// save flag that db backup was successful and filename of backup db
	//do_write_db = true;
	//ar << do_write_db;
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(load, blue_sky::sql_well)
	// load filename
	ar >> (std::string&)t.file_name;
	// and flag if we should restore db
	bool do_load_db;
	ar >> do_load_db;
	if(!do_load_db) return;

	// backup db handle
	sqlite3* bu_db;
	int opres;

	// load backup db fname
	std::string db_basename, db_fname;
	ar >> db_basename >> db_fname;

	// check if hdm serializer saved filename for us
	kernel::idx_dt_ptr kdt = BS_KERNEL.pert_idx_dt(BS_KERNEL.find_type("hdm").td_);
	if(kdt && kdt->size< std::string >()) {
		std::string prj_path = kdt->ss< std::string >(0);
		// try to open db relative to project path
		std::string db_relname = prj_path + PATHSEP + db_basename;
		opres = sqlite3_open(db_relname.c_str(), &bu_db);
		if(opres == SQLITE_OK) {
			// save handle to open db
			t.db = bu_db;
			return;
		}
		else
			sqlite3_close(bu_db);
	}

	opres = sqlite3_open(db_fname.c_str(), &bu_db);
	if(opres != SQLITE_OK) {
		sqlite3_close(bu_db);
		return;
	}
	// just save handle to db into sql_well and we're done
	t.db = bu_db;
BLUE_SKY_CLASS_SRZ_FCN_END

// generate serialize() function that uses save & load
//BLUE_SKY_CLASS_SERIALIZE_SPLIT(blue_sky::sql_well)
BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, blue_sky::sql_well)
	// register conversion to base iface
	boser::bs_void_cast_register(
		static_cast< sql_well* >(NULL),
		static_cast< well_pool_iface* >(NULL)
	);
	// split into save/load
	boser::split_free(ar, t, version);
BLUE_SKY_CLASS_SRZ_FCN_END

// instantiate serialization code
BLUE_SKY_TYPE_SERIALIZE_DECL(blue_sky::sql_well)
BLUE_SKY_TYPE_SERIALIZE_IMPL(blue_sky::sql_well)

