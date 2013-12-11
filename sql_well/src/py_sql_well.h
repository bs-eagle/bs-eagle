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
namespace blue_sky
  {
  namespace python
    {

  PY_EXPORTER (py_sql_well_exporter, default_exporter)
    .def ("open_db",                            &T::open_db,
        args ("file_name"), "Open database")
    .def ("close_db",                           &T::close_db,
        args (""), "Close database")
    .def ("fill_db",                           &T::fill_db,
        args (""), "Fill database")
    .def ("create_db_struct",                   &T::create_db_struct,
        args (""), "Create data base structure")
    .def ("add_well",                           &T::add_well,
        args ("well_name"), "Add new well to the storage")
    .def ("add_branch_gis",                     &T::add_branch_gis,
        args ("well_name", "branch_name", "gis"), "Add gis to the well branch")
    .def ("get_branch_gis",                     &T::get_branch_gis,
        args ("well_name", "branch_name"),  "Get gis of well branch")
    .def ("get_wlog_names",                     &T::get_wlog_names,
        args ("well_name", "branch_name"), "Get list of custom well logs")
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
    //.def ("d2date",                               &T::d2date,
    //  args ("d"), "convert d to list (year, month, day, hour, minute, second)")
    //.def ("date2d",                               &T::date2d,
    //  args ("year", "month", "day", "hour", "minute", "second"), "date and time to double")
    //.def ("d2str",                                &T::d2str,
    //  args ("d"), "convert date from double to str")
    //.def ("t2str",                                &T::t2str,
    //  args ("d"), "convert time from double to str")
    .def ("__str__",                            &T::py_str)
  PY_EXPORTER_END;

  //! export matrices to python
  void py_export_sql_well ();


  } // namespace python
} // namespace blue_sky

#endif // #ifdef BSPY_EXPORTING_PLUGIN
#endif /* end of include guard: PY_sql_well_O4VRW03K */

