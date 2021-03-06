/**
 * @file py_sql_well.h
 * @brief python wraper for #sql_well
 * @author Oleg Borschuk
 * @version
 * @date 2011-07-29
 */
#ifndef PY_SQL_WELL_O4VRW03K

#define PY_SQL_WELL_O4VRW03K


#include <string>
#include "well_pool_iface.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "dummy_base.h"
#include "construct_python_object.h"
#include "make_me_happy.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include "export_python_wrapper.h"

#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky { namespace python {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(add_branch_gis_ol, add_branch_gis, 3, 5)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(get_branch_gis_ol, get_branch_gis, 2, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(get_wlog_names_ol, get_wlog_names, 2, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(delete_well_log_ol, delete_well_log, 2, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(add_branch_ol, add_branch, 2, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(open_db_ol, open_db, 1, 2)

PY_EXPORTER (py_sql_well_exporter, default_exporter)
	.def ("open_db",                            &T::open_db, open_db_ol())
	//	args ("file_name"), "Open database")
	.def ("close_db",                           &T::close_db,
		args (""), "Close database")
	.def ("fill_db",                           &T::fill_db,
		args (""), "Fill database")
	.def ("create_db_struct",                   &T::create_db_struct,
		args (""), "Create data base structure")
	.def ("add_well",                           &T::add_well,
		args ("well_name"), "Add new well to the storage")
	.def ("add_branch_gis",                     &T::add_branch_gis, add_branch_gis_ol())
	.def ("get_branch_gis",                     &T::get_branch_gis, get_branch_gis_ol())
	.def ("get_wlog_names",                     &T::get_wlog_names, get_wlog_names_ol())
	.def ("add_branch_traj",                     &T::add_branch_traj,
		args ("well_name", "branch_name", "traj"), "Add traj to the well branch")
	.def ("update_branch_traj",                     &T::add_branch_traj,
		args ("well_name", "branch_name", "traj"), "Update traj to the well branch")
	.def ("get_branch_traj",                     &T::get_branch_traj,
		args ("well_name", "branch_name"),  "Get traj of well branch")
	.def ("get_well_names",                     &T::get_well_names,
		args (""), "Return well namew")
	.def ("prepare_sql",                        &T::prepare_sql,
		args ("SQL_str"), "Prepare SQL")
	.def ("step_sql",                           &T::step_sql,
		args (""), "Step SQL")
	.def ("finalize_sql",                       &T::finalize_sql,
		args (""), "Finalize SQL")
	.def ("get_sql_int",                        &T::get_sql_int,
		args ("Column"), "Return int value of column")
	.def ("get_sql_real",                       &T::get_sql_real,
		args ("Column"), "Return real value of column")
	.def ("get_sql_bool",                       &T::get_sql_bool,
		args ("Column"), "Return bool value of column")
	.def ("get_sql_str",                        &T::get_sql_str,
		args ("Column"), "Return text value of column")
	.def ("get_sql_exist",                      &T::get_sql_exist,
		args ("Column"), "Return existance of value in column")
	.def ("get_table",                          &T::get_table,
		args ("Table_name", "Column list", "Filter"), "Return nparray of values in columns list")
	.def ("exec_sql",                           &T::exec_sql,
		args ("Sql_string"), "Execute SQL")
	.def ("merge_with_db",                      &T::merge_with_db,
		args ("DB_name"), "Merge two DBs")
	.def ("exec_sql_and_return_rowid",          &T::exec_sql_and_return_rowid,
		args ("Sql_string"), "Execute SQL and return ROWID")
	.def ("read_from_ascii_file",               &T::read_from_ascii_file,
		args ("File_name", "starting_date"), "Read data from ascii file")
	.def ("save_to_bos_ascii_file",             &T::save_to_bos_ascii_file,
		args ("File_name", "h5_pool", "hdm_prop"),                  "Save data to ascii file in BOS format")
	.def ("insert_or_update",                   &T::insert_or_update,
		args ("select", "insert", "update"), "Insert or update data")
	.def ("backup_to_file",                     &T::backup_to_file,
		args ("filename"), "Backup memory DB to disk")
	.def ("__str__",                            &T::py_str)
	.def("rename_well_log", &T::rename_well_log)
	.def("delete_well_log", &T::delete_well_log, delete_well_log_ol())
	.def ("get_branches_names",                 &T::get_branches_names,
		args (""), "Return names of branches for given well")
	.def("add_branch", &T::add_branch, add_branch_ol())
	.def("rename_branch", &T::rename_branch)
	.def("delete_branch", &T::delete_branch)
PY_EXPORTER_END;

//! export matrices to python
void py_export_sql_well ();


	} // namespace python
} // namespace blue_sky

#endif // #ifdef BSPY_EXPORTING_PLUGIN
#endif /* end of include guard: PY_sql_well_O4VRW03K */

