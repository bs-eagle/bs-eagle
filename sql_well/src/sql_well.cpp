/**
 * @file sql_well.cpp
 * @brief implementation of frac storage
 * @author Oleg Borschuk
 * @version
 * @date 2011-07-29
 */

// d REAL NOT NULL REFERENCES dates(d) ON UPDATE CASCADE ON DELETE CASCADE,
//
// d REAL NOT NULL,

#include "bs_kernel.h"
#include "bs_misc.h"
#include "sql_well.h"

#include <iostream>
#include <stdio.h>
#include <ctype.h>
#include <boost/lexical_cast.hpp>
//#include <boost/format.hpp>

#ifdef BSPY_EXPORTING_PLUGIN

#include <boost/python.hpp>
using namespace boost::python;

#endif //BSPY_EXPORTING_PLUGIN

#define DB_FORMAT_PROP L"DB_format"

using namespace boost;

namespace blue_sky {

// blobs extracter
template< class T >
struct sql_well::extract_blob {
	typedef std::vector< T > res_t;

	static res_t go(sql_well& sqw, const int blob_idx) {
		typedef std::vector< T > res_t;

		// extract blob
		ulong n = sqlite3_column_bytes (sqw.stmp_sql, blob_idx);
		res_t res(n);
		const T* b = (const T*)sqlite3_column_blob (sqw.stmp_sql, blob_idx);
		std::copy(b, b + n, res.begin());
		return res;
	}
};

template< >
struct sql_well::extract_blob< std::string > {
	typedef std::string res_t;

	static std::string go(sql_well& sqw, const int blob_idx) {
		// extract blob
		std::string s;
		const int n = sqlite3_column_bytes (sqw.stmp_sql, blob_idx);
		const char *b = (const char *)sqlite3_column_blob (sqw.stmp_sql, blob_idx);
		s.assign(b, n);
		return s;
	}
};

// hidden details
namespace {
static const char* spaces = " \n\r\t";

std::string trim(const std::string& ss) {
	std::string s = ss;
	while(s.size() > 0 && strchr(spaces, s[0]) != NULL)
		s.erase(s.begin());
	while(s.size() > 0 && strchr(spaces, s[s.size() - 1]) != NULL)
		s.erase(s.size() - 1);
	return s;
}

std::string to_upper(const std::string& s) {
	std::string us(s.size(), ' ');
	std::transform(s.begin(), s.end(), us.begin(), ::toupper);
	return us;
}

std::string to_lower(const std::string& s) {
	std::string us(s.size(), ' ');
	std::transform(s.begin(), s.end(), us.begin(), ::tolower);
	return us;
}

} // eof hidden namespace

sql_well::sql_well(bs_type_ctor_param) {
	db = 0;
	stmp_sql = 0;
	fr_file = 0;

}

sql_well::sql_well(const sql_well& rhs)
	: bs_refcounter (), file_name(rhs.file_name), db(rhs.db)
{
	//*this = rhs;
	// don't copy pending statement
	stmp_sql = 0;
	fr_file = 0;
}

sql_well::~sql_well() {
	if (db)
		close_db ();
}

int sql_well::add_well(const std::string &well_name) {
	if (!db)
		return -2;
	if (stmp_sql)
		finalize_sql ();

	int rc = 0;
	char *zErrMsg = 0;
	char buf[2048];

	sprintf (buf, "INSERT INTO wells (name) VALUES('%s');", well_name.c_str ());

	rc = sqlite3_exec (db, buf, NULL, 0, &zErrMsg);
	if(rc != SQLITE_OK) {
		fprintf (stderr, "SQL error (add_well): %s\n", zErrMsg);
		sqlite3_free (zErrMsg);
		return -1;
	}
	return 0;
}

int sql_well::delete_well(const std::string& well_name) {
	// remove well logs
	// assume that all DBs are in new format and contain well_logs table
	exec_sql("DELETE FROM well_logs WHERE well_name = '" + well_name + "'");
	// delete from other HDM-related tables
	exec_sql("DELETE FROM completions WHERE well_name = '" + well_name + "'");
	exec_sql("DELETE FROM fractures WHERE well_name = '" + well_name + "'");
	exec_sql("DELETE FROM well_hist WHERE well_name = '" + well_name + "'");
	exec_sql("DELETE FROM well_res WHERE well_name = '" + well_name + "'");
	exec_sql("DELETE FROM wells_in_group WHERE well_name = '" + well_name + "'");
	// remove vranches
	exec_sql("DELETE FROM branches WHERE well_name = '" + well_name + "'");
	// remove well
	int res = exec_sql("DELETE FROM wells WHERE name = '" + well_name + "'");
	// clean deleted entries from DB
	exec_sql("VACUUM");
	return res;
}

sql_well::list_t sql_well::get_well_names() const {
	list_t lst;
	if (!db)
		return lst;

	ulong rc = 0;
	//char *zErrMsg = 0;
	const char *ttt;
	sqlite3_stmt *stmp;
	std::string sql = "SELECT name FROM wells ORDER BY name ASC";
	rc = sqlite3_prepare_v2 (db, sql.c_str (), sql.length () + 1, &stmp, &ttt);
	if(rc) {
		fprintf (stderr, "Can't make select: %s\n", sqlite3_errmsg (db));
		return lst;
	}
	while (sqlite3_step (stmp) == SQLITE_ROW) {
		// UPDATE
		lst.push_back (std::string ((const char *)sqlite3_column_text (stmp, 0)));
	}
	sqlite3_finalize (stmp);
	return lst;
}

sql_well::list_t sql_well::get_branches_names(const std::string& well_name) const {
	list_t lst;
	if (!db)
		return lst;

	// FIXME: union branch names taken from branches and from well_logs tables!
	// selection from well_logs should be removed after correct multibranch implementation is ready
	ulong rc = 0;
	const char *ttt;
	sqlite3_stmt *stmp;
	std::string sql = "SELECT branch_name FROM branches WHERE well_name = '" + well_name +
		"' UNION SELECT DISTINCT branch_name FROM well_logs WHERE well_name = '" + well_name + "'";
	rc = sqlite3_prepare_v2 (db, sql.c_str (), sql.length () + 1, &stmp, &ttt);
	if(rc) {
		fprintf (stderr, "Can't make select: %s\n", sqlite3_errmsg (db));
		return lst;
	}
	while (sqlite3_step (stmp) == SQLITE_ROW) {
		// UPDATE
		lst.push_back (std::string ((const char *)sqlite3_column_text (stmp, 0)));
	}
	sqlite3_finalize (stmp);
	return lst;
}

int sql_well::add_branch(
	const std::string &wname, const std::string &branch, const std::string &parent, t_double md
) {
	std::string sql = "INSERT OR REPLACE INTO branches (well_name, branch_name, parent, md)";
	sql += " VALUES ('" + wname + "', '" + branch + "', '" + parent + "', " +
		boost::lexical_cast< std::string >(md) + ")";
	return exec_sql(sql);
}

int sql_well::add_branch_gis (const std::string &wname, const std::string &branch,
	sp_gis_t g, std::string wlog_name, uint wlog_type,
	bool replace_existing
) {
	// helper
	// bind well log data blob and execute PREPARED sql statement
	struct dump_wlog_data {
		static int go(sql_well& sqw, const sp_gis_t& g) {
			// bind well log data
			const std::string log_data = g->to_str();
			if (sqlite3_bind_blob(sqw.stmp_sql, 1, &log_data.c_str()[0], int(log_data.size()), SQLITE_STATIC)) {
				fprintf (stderr, "Can't make select: %s\n", sqlite3_errmsg (sqw.db));
				sqw.finalize_sql();
				return -3;
			}

			// exec query
			sqw.step_sql();
			sqw.finalize_sql();
			return 0;
		}
	};

	if (!db || !g)
		return -1;
	if (stmp_sql)
		finalize_sql ();

	std::string q;
	if(wlog_name.size() == 0) {
		// 0. Check wlog property indicating if it has been converted to new format before
		sp_prop_t log_prop = g->get_prop();
		std::vector< std::wstring > names = log_prop->get_names_i();
		// ensure DB format prop exists
		if(std::find(names.begin(), names.end(), DB_FORMAT_PROP) == names.end()) {
			log_prop->add_property_i(0, DB_FORMAT_PROP, L"");
		}
		// mark that we store logs in new format
		log_prop->set_i(DB_FORMAT_PROP, 1);

		// 1. convert old representation to new
		// split logs table into sequence of DEPT-LOG data tables
		sp_table_t log_data = g->get_table();
		names = log_data->get_col_names();
		// process only tables with >= 2 columns
		if(names.size() < 2)
			return -1;

		// search DEPTH values
		// take first column by default
		ulong dept_idx = 0;
		for(ulong i = 0; i < names.size(); ++i) {
			const std::string cur_col = trim(to_lower(wstr2str(names[i], "utf-8")));
			if(cur_col.find("dept", 0, 4) != std::string::npos) {
				dept_idx = i;
				break;
			}
		}
		// extract depth vector
		spv_double dept_data = BS_KERNEL.create_object(v_double::bs_type());
		dept_data->init(log_data->get_col_vector(dept_idx));

		// prepare string representation of properties from parent gis
		// TODO: find better way to copy props from one object to another
		const std::string log_prop_dump = g->get_prop()->to_str();

		// loop over all well log columns
		spv_double cur_log = BS_KERNEL.create_object(v_double::bs_type());
		sp_gis_t cur_gis = BS_KERNEL.create_object("gis");
		for(ulong i = 0; i < names.size(); ++i) {
			// skip depth
			if(i == dept_idx) continue;
			// create and fill new gis object
			cur_gis->get_prop()->from_str(log_prop_dump);
			sp_table_t cur_data = cur_gis->get_table();
			cur_data->init(0, 2);
			cur_data->add_col_vector(0, names[dept_idx], dept_data);
			cur_log->init(log_data->get_col_vector(i));
			cur_data->add_col_vector(1, names[i], cur_log);

			// write it to DB
			add_branch_gis(
				wname, branch, cur_gis, wstr2str(names[i], "utf-8"), wlog_type, replace_existing
			);
		}

		// 2. Write gis in old format if no wlog_name specified
		// TODO: deprecate and disable writing to branches table at all
		q = "UPDATE ";
		//if(!replace_existing)
		//  q += "OR IGNORE ";
		q += "branches SET well_log = ?1 WHERE well_name = '";
		q += wname + "' AND branch_name = '" + branch + "'";

		if(prepare_sql(q) < 0)
			return -1;
		// exec query
		return dump_wlog_data::go(*this, g);
	}

	// new implementation writes to separate well logs table
	// format query
	q = "INSERT OR ";
	if(replace_existing)
		q += "REPLACE";
	else
		q += "IGNORE";
	q += " INTO well_logs (well_name, branch_name, wlog_name, wlog_type, wlog_data) VALUES ('";
	q += wname + "', '" + branch + "', '" + wlog_name + "', " +
		boost::lexical_cast< std::string >(wlog_type) + ", ?1)";

	// prepare query
	if(prepare_sql(q) < 0)
		return -1;
	// exec query
	return dump_wlog_data::go(*this, g);
}

std::vector< std::string > sql_well::get_wlog_names(
	const std::string &wname, const std::string &branch, uint wlog_type
) {
	std::vector< std::string > res;
	if (!db)
		return res;

	std::string q = "SELECT wlog_name FROM well_logs WHERE well_name = '" + wname +
		"' AND branch_name = '" + branch + "' AND wlog_type = " +
		boost::lexical_cast< std::string >(wlog_type);

	// exec sql
	if(prepare_sql(q) < 0) {
		finalize_sql();
		return res;
	}

	while(step_sql() == 0) {
		res.push_back(get_sql_str(0));
	}

	finalize_sql();
	return res;
}

sql_well::sp_gis_t sql_well::get_branch_gis(
	const std::string &wname, const std::string &branch,
	std::string wlog_name, uint wlog_type
) {
	// result
	sp_gis_t sp_gis;

	if (!db)
		return sp_gis;
	finalize_sql();

	std::string q;
	const std::string select_filter = " WHERE well_name = '" + wname +
		"' AND branch_name = '" + branch + "'";
	// what blobs can we extract from result?
	bool has_data = false;

	if(wlog_name.size()) {
		// check well_logs table
		// format query
		q = "SELECT wlog_data FROM well_logs" + select_filter +
			" AND wlog_name = '" + wlog_name + "' AND wlog_type = " +
			boost::lexical_cast< std::string >(wlog_type);
		// exec sql
		if(prepare_sql(q) == 0 && step_sql() == 0) {
			has_data = true;
		}
	}
	else {
		// make old-fashioned query to branches table
		// format query
		q = "SELECT well_log FROM branches" + select_filter;
		// exec sql
		if(prepare_sql(q) == 0 && step_sql() == 0) {
			has_data = true;
		}
	}

	// 3. read traj BLOB
	if(has_data) {
		// leave this for debugging purposes
		std::cout << "READ WELL LOG: " << wname;

		// extract well log data
		q = extract_blob< std::string >::go(*this, 0);
		if(!q.empty()) {
			// we have some to read
			sp_gis = BS_KERNEL.create_object ("gis");
			sp_gis->from_str(q);
		}
		std::cout << ", DATA = " << q.size() << std::endl;
	}
	finalize_sql();

	return sp_gis;
}

bool sql_well::rename_well_log(
	const std::string &wname, const std::string &branch,
	const std::string& old_name, const std::string& new_name
) {
	if (!db)
		return false;
	if (stmp_sql)
		finalize_sql ();

	const std::string q = "UPDATE well_logs SET wlog_name = '" + new_name +
		"' WHERE well_name = '" + wname + "' AND branch_name = '" + branch +
		"' AND wlog_name = '" + old_name + "'";

	bool res = (exec_sql(q) == 0);

	// try to rename also in old-fashioned logs table inside branches
	sp_gis_t g = get_branch_gis(wname, branch);
	if(!g) return res;

	sp_table_t wlog_data = g->get_table();
	const std::wstring wold_name = str2wstr(old_name, "utf-8");
	const ulong n_wlogs(wlog_data->get_n_cols());
	for(ulong i = 0; i < n_wlogs; ++i) {
		if(wlog_data->get_col_name(i) != wold_name) continue;
		// we found a matching column
		// rename it
		wlog_data->set_col_name(i, str2wstr(new_name, "utf-8"));
		// and write result to DB
		add_branch_gis(wname, branch, g);
		return true;
	}

	return res;
}

bool sql_well::delete_well_log(
	const std::string &wname, const std::string &branch, std::string wlog_name
) {
	if (!db)
		return false;
	if (stmp_sql)
		finalize_sql ();

	std::string q;
	const std::string select_filter = " WHERE well_name = '" + wname +
		"' AND branch_name = '" + branch + "'";
	bool res = false;

	if(wlog_name.size() == 0) {
		// old-fashioned query
		q = "UPDATE branches SET well_log = NULL" + select_filter;
		res = (exec_sql(q) == 0);
	}
	else {
		// delete from well_logs table
		q = "DELETE FROM well_logs" + select_filter +
			" AND wlog_name = '" + wlog_name + "'";
		res = (exec_sql(q) == 0);
	}

	// try to delete also from old-fashioned logs table inside branches
	sp_gis_t g = get_branch_gis(wname, branch);
	if(!g) return res;

	sp_table_t wlog_data = g->get_table();
	// delete corresponding column
	const std::wstring w_name = str2wstr(wlog_name, "utf-8");
	const ulong n_wlogs(wlog_data->get_n_cols());
	for(ulong i = 0; i < n_wlogs; ++i) {
		if(wlog_data->get_col_name(i) != w_name) continue;
		// we found a matching column
		// delete it
		wlog_data->remove_col(i);
		// and write result to DB
		add_branch_gis(wname, branch, g);
		break;
	}

	// clean deleted entries from DB
	exec_sql("VACUUM");
	return res;
}

int sql_well::add_branch_traj(
	const std::string &wname, const std::string &branch, sp_traj_t t
) {
	if (!db || !t)
		return -1;
	if (stmp_sql)
		finalize_sql ();

	const std::string q = "UPDATE branches SET traj = ?1 WHERE well_name = '" +
		wname + "' AND branch_name = '" + branch + "'";

	// prepare query
	if(prepare_sql(q) < 0)
		return -1;

	// bind well log data
	const std::string traj_data = t->to_str();
	if (sqlite3_bind_blob (stmp_sql, 1, &traj_data.c_str()[0], int(traj_data.size()), SQLITE_STATIC)) {
		fprintf (stderr, "Can't make select: %s\n", sqlite3_errmsg (db));
		finalize_sql();
		return -3;
	}

	// exec query
	step_sql();
	finalize_sql();

	return 0;
}

sql_well::sp_traj_t sql_well::get_branch_traj(
	const std::string &wname, const std::string &branch
) {
	sp_traj_t sp_traj;
	if (!db)
		return sp_traj;
	finalize_sql();

	// format sql
	std::string q = "SELECT traj, parent, md FROM branches WHERE well_name = '" + wname +
		"' AND branch_name = '" + branch + "'";

	// exec sql
	if(prepare_sql(q) == 0 && step_sql() == 0) {
		// leave this for debugging purposes
		std::cout << "READ WELL TRAJ: (" << wname << ", " << branch << ')';
		// extract trajectory, parent, md
		q = extract_blob< std::string >::go(*this, 0);
		std::string parent = get_sql_str(1);
		t_double md = get_sql_real(2);
		finalize_sql();

		if(!q.empty()) {
			sp_traj = BS_KERNEL.create_object ("traj");
			sp_traj->from_str(q);
			std::cout << ", ";
		}
		// check if branch contains reference to parent
		else if(!parent.empty()) {
			std::cout << " -> PARENT: ";
			sp_traj = get_branch_traj(wname, parent);
			if(md <= 0)
				return sp_traj;

			// if positive MD is set, return parent trajectory after given MD
			// obtain parent traj data
			sp_table_t ptraj = sp_traj->get_table();
			// check if MD column exists
			std::vector< std::wstring > col_names = ptraj->get_col_names();
			int md_col = -1;
			for(ulong i = 0; i < col_names.size(); ++i) {
				if(col_names[i] == L"MD" || col_names[i] == L"md") {
					md_col = int(i);
					break;
				}
			}
			// cut trajectory below md
			// DEBUG
			//std::cout << "MD col = " << md_col << std::endl;
			if(md_col >= 0) {
				const t_double* p_md = ptraj->get_col_ptr(md_col);
				for(ulong i = ptraj->get_n_rows() - 1; i < ptraj->get_n_rows(); --i) {
					if(p_md[i] < md) {
						//std::cout << "cut row " << i << ", MD = " << p_md[i] << std::endl;
						ptraj->remove_row(i);
					}
				}
			}
		}
		std::cout << "DATA = " << sp_traj->get_table()->get_n_rows() << std::endl;
	}
	else
		finalize_sql();

	return sp_traj;
}

bool sql_well::delete_branch(
	const std::string &wname, const std::string &branch
) {
	if (!db)
		return false;
	if (stmp_sql)
		finalize_sql ();

	std::string q;
	const std::string select_filter = " WHERE well_name = '" + wname +
		"' AND branch_name = '" + branch + "'";
	bool res = false;

	// remove branch-related data in misc tables
	// well logs
	q = "DELETE FROM well_logs" + select_filter;
	exec_sql(q);
	// TODO: get rid of following tables completely
	// completions
	q = "DELETE FROM completions" + select_filter;
	exec_sql(q);
	// fractures
	q = "DELETE FROM fractures" + select_filter;
	exec_sql(q);

	// and finally delete branch
	q = "DELETE FROM branches" + select_filter;
	res = (exec_sql(q) == 0);

	// clean deleted entries from DB
	exec_sql("VACUUM");
	return res;
}

bool sql_well::rename_branch(
	const std::string& wname, const std::string& old_branch, const std::string& new_branch
) {
	if (!db)
		return false;
	if (stmp_sql)
		finalize_sql ();

	std::string q;
	const std::string update_filter = " SET branch_name = '" + new_branch +
		"' WHERE well_name = '" + wname + "' AND branch_name = '" + old_branch + "'";

	// rename branch in misc tables
	// well logs
	q = "UPDATE well_logs" + update_filter;
	exec_sql(q);
	// TODO: get rid of following tables completely
	// completions
	q = "UPDATE completions" + update_filter;
	exec_sql(q);
	// fractures
	q = "UPDATE fractures" + update_filter;
	exec_sql(q);

	// and finally delete branch
	q = "UPDATE branches" + update_filter;
	return (exec_sql(q) == 0);
}

#ifdef BSPY_EXPORTING_PLUGIN

std::string sql_well::py_str() const {
	std::stringstream s;
	s << file_name << "\n";
	return s.str ();
}

#endif //BSPY_EXPORTING_PLUGIN

/////////////////////////////////BS Register
/////////////////////////////////Stuff//////////////////////////

  BLUE_SKY_TYPE_STD_CREATE (sql_well);
  BLUE_SKY_TYPE_STD_COPY (sql_well);

  BLUE_SKY_TYPE_IMPL(sql_well,  well_pool_iface, "sql_well", "sql_well storage", "realization of well sql_well storage");

}  // blue_sky namespace

