/// @file h5_pool_serialize.cpp
/// @brief h5_pool serialization implementation
/// @author uentity
/// @version
/// @date 24.01.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

#include "bs_serialize.h"
#include "h5_pool.hpp"
#include "bs_prop_base.h"

#include <boost/serialization/array.hpp>
#include <boost/serialization/map.hpp>
#include <boost/lexical_cast.hpp>

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid_io.hpp>

#include <algorithm>
#include <iostream>
#include <fstream>

// path separator
#ifdef UNIX
#define PATHSEP '/'
#else
#define PATHSEP '\\'
#endif

using namespace blue_sky;
namespace boser = boost::serialization;
namespace bu = boost::uuids;

// helper functions
namespace {
using namespace std;

hid_t demand_h5_group(hid_t file_id, const std::string& g_name, bool* is_created = NULL) {
	if(is_created)
		*is_created = false;
	// try to open group
	hid_t g_id = -1;
	if(H5Lexists(file_id, g_name.c_str(), 0))
		g_id = H5Gopen(file_id, g_name.c_str());
	if(g_id < 0) {
		// try to create group
		g_id = H5Gcreate(file_id, g_name.c_str (), -1);
		if(g_id < 0)
			bs_throw_exception(
				boost::format ("h5_pool_serialize: Can't create group: %s") % g_name
			);
		if(is_created)
			*is_created = true;
	}
	return g_id;
}

std::string basename(const std::string& fname) {
	std::size_t sep_pos = fname.rfind(PATHSEP);
	if(sep_pos == std::string::npos)
		return fname;
	else
		return fname.substr(sep_pos + 1, std::string::npos);
}

// helper struct to find group name that match best with given name
struct grp_match {
	// max difference between existing and given group names
	const float tol_;

	// file handle
	hid_t f_;
	// group names contained in file
	set< string > gnames_;
	typedef set< string >::const_iterator gn_iterator;

	// resolves group name -> ID
	typedef pair< hid_t, const string& > g_info;
	typedef map< string, g_info > resolver_t;
	typedef resolver_t::iterator r_iterator;
	resolver_t r_;

	grp_match(float diff_tol) : tol_(diff_tol) {};

	void init(hid_t file_id) {
		f_ = file_id;
		discover_groups();
	}

	g_info& operator[](const string& gkey) {
		typedef set< string >::const_iterator n_iterator;
		typedef pair< string::const_iterator, string::const_iterator > mis_pair;

		r_iterator phid = r_.find(gkey);
		if(phid != r_.end())
			return phid->second;

		// if we are here then group name wasn't resolved yet
		// first try to find exact match
		string res_name;
		n_iterator gn = gnames_.find(gkey);
		if(gn != gnames_.end()) {
			res_name = *gn;
		}
		else {
			// try to find group with best match
			size_t match_len = 0;
			for(set< string >::const_iterator gn = gnames_.begin(), end = gnames_.end(); gn != end; ++gn) {
				size_t len = min(gkey.size(), gn->size());
				mis_pair mp = mismatch(gn->begin(), gn->begin() + len, gkey.begin());
				size_t cur_match_len = mp.first - gn->begin();
				if(cur_match_len > match_len) {
					res_name = *gn;
					match_len = cur_match_len;
				}
			}

			// if difference is more than tolerance - then try to create group with key name
			if(gkey.size() - match_len > gkey.size() * tol_)
				res_name = gkey;
		}

		// open/create group
		bool is_created;
		hid_t res_hid = demand_h5_group(f_, res_name.c_str(), &is_created);
		// if new group was created - update dict
		if(is_created)
			gnames_.insert(res_name);

		return r_.insert(make_pair(gkey, g_info(res_hid, *gnames_.find(res_name)))).first->second;
	}

	void discover_groups() {
		gnames_.clear();
		hid_t root_gid = H5Gopen(f_, "/");
		H5Literate(
			root_gid, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, &(grp_match::enumerator), this
		);
	}

	static herr_t enumerator(hid_t root_gid, const char* gname, const H5L_info_t*, void* m) {
		grp_match& self = *static_cast< grp_match* >(m);
		self.gnames_.insert(string(gname));
		return 0;
	}
};

}

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, h5_pair)
	ar & t.n_dims & t.size;
	ar & t.py_dims & t.h5_dims & t.src_dims;
	ar & t.var_dims;
	// save hid to remember if corresponding item was open
	ar & t.dset; ar & t.dtype; ar & t.plist;
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(save, h5_pool)
	// save flag if h5 is open
	bool is_open = (t.file_id >= 0);
	ar << is_open;
	//ar << t.fname;

	std::string h5_fname = t.fname;
	std::string h5_basename = (const std::string&)basename(t.fname);

	// check if we should explicitly copy h5 storage instead of referencing existing file
	kernel::idx_dt_ptr kdt = BS_KERNEL.pert_idx_dt(BS_KERNEL.find_type("hdm").td_);
	const bool do_deep_copy = is_open && kdt && kdt->size< std::string >() > 2;
	if(do_deep_copy) {
		// randomly generate h5 filename
		// bu::uuid cpy_uuid = bu::random_generator()();
		// filename = project path + basename
		h5_basename = kdt->ss< std::string >(1) + kdt->ss< std::string >(2) + ".h5";
		h5_fname = kdt->ss< std::string >(0) + PATHSEP + h5_basename;
	}

	// save full filename & basename
	// NOTE: use UTF-8 encoding for filenames
	ar << (const std::string&)str2ustr(h5_fname);
	ar << (const std::string&)str2ustr(h5_basename);

	// we should only dump group names
	typedef h5_pool::map_hid_t::const_iterator g_iterator;
	const std::size_t sz = t.group_id.size();
	ar << sz;
	for(g_iterator pg = t.group_id.begin(); pg != t.group_id.end(); ++pg)
		ar << pg->first;

	// dump nodes
	typedef h5_pool::map_t::const_iterator m_iterator;
	// name of group that will contain committed dtypes
	const std::string dtype_grp = "_dtype";
	hid_t dtype_grp_hid = -1;
	std::string g_name, dtype_name;
	t_ulong dtype_idx = 0;
	vector< char > dtype_ex_name(16);

	ar << (const std::size_t&)t.h5_map.size();
	for(m_iterator pm = t.h5_map.begin(); pm != t.h5_map.end(); ++pm) {
		ar << pm->first;
		ar << pm->second;

		// find corresponding group name by id
		g_name.clear();
		for(g_iterator pg = t.group_id.begin(); pg != t.group_id.end(); ++pg) {
			if(pg->second == pm->second.group_id) {
				g_name = pg->first;
				break;
			}
		}
		// and save it
		ar << g_name;

		// save dataset type info
		hid_t dtype = pm->second.dtype;
		if(is_open && dtype >= 0) {
			// create group to hide all dtypes
			if(dtype_grp_hid < 0) {
				try {
					dtype_grp_hid = demand_h5_group(t.file_id, dtype_grp);
				}
				catch(bs_exception& e) {
					BSERR
						<< "h5 pool serialize: WARNING! can't create group for storing data types: "
						<< e.what() << bs_end;
				}
				catch(...) {
					throw;
				}
			}

			dtype_name.clear();
			// check if dtype already has associated name
			if(H5Tcommitted(dtype)) {
				ssize_t sz = H5Iget_name(dtype, NULL, 0);
				if(dtype_ex_name.size() < std::size_t(sz + 1))
					dtype_ex_name.resize(sz + 1);
				H5Iget_name(dtype, &dtype_ex_name[0], sz + 1);
				dtype_name = &dtype_ex_name[0];
			}
			else if(dtype_grp_hid >= 0) {
				while(1 > 0) {
					// commit type into H5 as named type
					dtype_name = dtype_grp + "/_dtype" + boost::lexical_cast< std::string >(dtype_idx++);
					if(H5Lexists(t.file_id, dtype_name.c_str(), 0)) {
						continue;
						//H5Ldelete(t.file_id, dtype_name.c_str(), 0);
					}
					H5Tcommit(t.file_id, dtype_name.c_str(), dtype);
					break;
				}
			}
			ar << dtype_name;

			//ar << (const H5T_class_t&)H5Tget_class(dtype);
			//ar << (const hid_t&)H5Tget_size(dtype);
		}
	}

	// other data
	ar << t.pool_dims;
	//ar << boser::make_array(&t.pool_dims, 3);
	ar << t.n_pool_dims;
	// flush buffers
	const_cast< h5_pool& >(t).flush();

	// explicitly copy src h5 if we're doing deep copy
	if(do_deep_copy) {
		// copy file via memory -- probably bad for big h5
		// first create new h5
		//hid_t cpy_id = H5Fcreate(cpy_fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT);
		//if(cpy_id < 0)
		//	bs_throw_exception(
		//		boost::format ("h5_pool_serialize: Can't create h5 file %s to store copy of source %s") % cpy_fname % t.fname
		//	);

		// make a copy of h5 file
		std::ifstream src_h5(t.fname.c_str(), std::ios::binary);
		std::ofstream dst_h5(h5_fname.c_str(), std::ios::trunc | std::ios::binary);
		dst_h5 << src_h5.rdbuf();
		dst_h5.flush();
	}
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(load, h5_pool)
	bool do_open;
	ar >> do_open;
	ar >> t.fname;
	// NOTE: expect filenames in unicode, convert to native locale
	t.fname = ustr2str(t.fname);

	// load base name
	std::string h5_basename;
	ar >> h5_basename;
	h5_basename = ustr2str(h5_basename);

	// open file
	if(do_open) {
		t.file_id = -1;
		// first try to open base name relative to project path
		kernel::idx_dt_ptr kdt = BS_KERNEL.pert_idx_dt(BS_KERNEL.find_type("hdm").td_);
		if(kdt && kdt->size< std::string >()) {
			std::string prj_path = kdt->ss< std::string >(0);
			h5_basename = prj_path + PATHSEP + h5_basename;
			t.file_id = H5Fopen(h5_basename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
			if(t.file_id >= 0)
				t.fname = h5_basename;
		}

		if(t.file_id < 0)
			t.file_id = H5Fopen(t.fname.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
		if(t.file_id < 0) // try to create file
			t.open_file(t.fname.c_str());
		// something really bad
		if(t.file_id < 0)
			bs_throw_exception(
				boost::format("Can't open and even create h5 file %s")
				% t.fname
			);
	}

	// create matcher object
	// note tolerance value = 0.5 by default
	// if set to zero only absolute matches accepted
	grp_match grp_dict(0.5);
	if(do_open)
		grp_dict.init(t.file_id);

	// build groups map
	std::size_t g_size;
	std::string g_name;
	ar >> g_size;
	for(std::size_t i = 0; i < g_size; ++i) {
		ar >> g_name;
		// add group to map
		if(do_open) {
			grp_match::g_info& gi = grp_dict[g_name];
			t.fill_map(gi.first);
			t.group_id[gi.second] = gi.first;
		}
		else
			t.group_id[g_name] = -1;
	}

	// read and build h5_map
	hid_t g_id;
	h5_pair p;
	std::string d_name, dtype_name;
	// type info
	//H5T_class_t type_class;
	//size_t type_size;
	// fill h5_map
	ar >> g_size;
	for(std::size_t i = 0; i < g_size; ++i) {
		ar >> d_name; // dataset name
		ar >> p;      // h5_pair
		ar >> g_name; // corresponding group name

		// if file is not open then just insert current node and continue
		if(!do_open) {
			p.group_id = p.dset = p.dspace = p.dtype = p.plist = -1;
			t.h5_map[d_name] = p;
			continue;
		}

		// read type info
		if(p.dtype >= 0) {
			ar >> dtype_name;
			//ar >> type_class;
			//ar >> type_size;
		}

		// init group id to first group in list
		// and first group should be always 'actual'
		g_id = -1;
		if(t.group_id.size())
			g_id = t.group_id.begin()->second;

		// try to open dataset
		// if we are lucky, then all hids are filled from open dataset
		hid_t dset_tmp = -1;
		if(p.dset >= 0 && g_name.size()) {
			// find corresponding group id
			g_id = grp_dict[g_name].first;

			// try to open dataset
			dset_tmp = H5Dopen(g_id, d_name.c_str());
			if(dset_tmp >= 0) {
				p.dset = dset_tmp;
				// open other attributes
				p.dspace = H5Dget_space(dset_tmp);
				if(p.dspace < 0)
					bs_throw_exception(
						boost::format("Can't get dataspace for dataset %s in group %d")
						% d_name % g_name
					);

				p.dtype = H5Dget_type(dset_tmp);
				if(p.dtype < 0)
					bs_throw_exception(
						boost::format("Can't get datatype for dataset %s in group %d")
						% d_name % g_name
					);
				// type open from dataset, so delete unused named entry
				if(dtype_name.size())
					H5Ldelete(t.file_id, dtype_name.c_str(), 0);

				p.plist = H5Dget_create_plist(dset_tmp);
				if(p.plist < 0)
					bs_throw_exception(
						boost::format("Can't get property list for dataset %s in group %d")
						% d_name % g_name
					);
			}
		}

		// if hids cannot be loaded from open dataset
		// try to create plist & dtype manually
		if(dset_tmp < 0) {
			p.plist = H5Pcreate(H5P_DATASET_CREATE);
			if(p.plist < 0)
				bs_throw_exception(
					boost::format("Can't create property for dataset %s in group %d")
					% d_name % g_name
				);

			// if we have to restore data type - try to open it
			if(p.dtype >= 0) {
				if(dtype_name.size()) {
					p.dtype = H5Topen(t.file_id, dtype_name.c_str());
					if(p.dtype < 0)
						BSERR << "h5_pool serialization: WARNING! Can't open data type "
							<< dtype_name << bs_end;
				}
				else {
					p.dtype = -1;
					BSERR << "h5_pool serialization: WARNING! Can't open unknown data type" << bs_end;
				}
			}

			//if(
			//	type_class == H5T_COMPOUND ||
			//	type_class == H5T_OPAQUE ||
			//	type_class == H5T_ENUM ||
			//	type_class == H5T_STRING
			//)
			//	p.dtype = H5Tcreate(type_class, type_size);
			//else
			//	p.dtype = H5Tcopy(p.dtype);

				//bs_throw_exception(
				//	boost::format("Can't create datatype of class %d and size %d")
				//	% type_class % type_size
				//);

			// if during saving dataset was open but we failed to open it now
			// try to create it
			if(p.dset >= 0 && g_name.size()) {
				// try to create dataspace
				p.dspace = H5Screate_simple(p.n_dims, p.h5_dims, NULL);
				if(p.dspace < 0)
					bs_throw_exception(
						boost::format("Can't create simple dataspace for dataset %s in group %d")
						% d_name % g_name
					);

				// and finally dataset
				p.dset = H5Dcreate(g_id, d_name.c_str(), p.dtype, p.dspace, p.plist);
				if(p.dset < 0)
					bs_throw_exception(
						boost::format("Can't create dataset %s in group %d")
						% d_name % g_name
					);
			}
			else
				p.dspace = -1;
		}

		// insert dataset node
		p.group_id = g_id;
		t.h5_map[d_name] = p;
	}

	ar >> t.pool_dims;
	ar >> t.n_pool_dims;
	t.flush();
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, h5_pool)
	// register conversion to base interface
	boser::bs_void_cast_register< h5_pool, h5_pool_iface >(
		static_cast< h5_pool* >(NULL),
		static_cast< h5_pool_iface* >(NULL)
	);
	boser::split_free(ar, t, version);
BLUE_SKY_CLASS_SRZ_FCN_END

BOOST_SERIALIZATION_ASSUME_ABSTRACT(h5_pool_iface)
BLUE_SKY_TYPE_SERIALIZE_DECL(h5_pool)
BLUE_SKY_TYPE_SERIALIZE_IMPL(h5_pool)

#ifdef BSPY_EXPORTING_PLUGIN
namespace blue_sky { namespace python {

BS_API_PLUGIN void py_export_h5_pool_serialize() {
	using namespace boost::python;

	std::string (*s2s_hpi)(const smart_ptr< h5_pool_iface, true >&) = &blue_sky::serialize_to_str< h5_pool_iface >;
	std::string (*s2s_pi)(const smart_ptr< h5_pool, true >&) =
		&blue_sky::serialize_to_str_indirect< h5_pool, h5_pool_iface >;
	smart_ptr< h5_pool_iface, true > (*sfs_hpi)(const std::string&) = &blue_sky::serialize_from_str< h5_pool_iface >;
	smart_ptr< h5_pool, true > (*sfs_pi)(const std::string&) =
		&blue_sky::serialize_from_str_indirect< h5_pool, h5_pool_iface >;

	def("serialize_to_str", s2s_hpi);
	def("serialize_to_str", s2s_pi);
	def("serialize_from_str", sfs_hpi);
	def("serialize_from_str", sfs_pi);
}

}} /* blue_sky::python */
#endif // #ifdef BSPY_EXPORTING_PLUGIN

