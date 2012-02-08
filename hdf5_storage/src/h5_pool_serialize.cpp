/// @file h5_pool_serialize.cpp
/// @brief h5_pool serialization implementation
/// @author uentity
/// @version 
/// @date 24.01.2012
/// @copyright This source code is released under the terms of
///            the BSD License. See LICENSE for more details.

//#include "bs_bos_core_data_storage_stdafx.h"

#include "h5_pool_serialize.h"
#include "bs_prop_base.h"

#include <boost/serialization/array.hpp>
#include <boost/serialization/map.hpp>

// path separator
#ifdef UNIX
#define PATHSEP '/'
#else
#define PATHSEP '\\'
#endif

using namespace blue_sky;
namespace boser = boost::serialization;

// helper functions
namespace {

hid_t demand_h5_group(hid_t file_id, const std::string& g_name) {
	// try to open group
	hid_t g_id = H5Gopen(file_id, g_name.c_str());
	if(g_id < 0) {
		// try to create group
		g_id = H5Gcreate(file_id, g_name.c_str (), -1);
		if(g_id < 0)
			bs_throw_exception(
				boost::format ("h5_pool_serialize: Can't create group: %s") % g_name
			);
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

}

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(serialize, h5_pair)
	ar & t.n_dims & t.size;
	ar & t.py_dims & t.h5_dims & t.src_dims;
	ar & t.var_dims & t.diff_from_base;
	// save hid to remember if corresponding item was open
	ar & t.dset; ar & t.dtype; ar & t.plist;
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(save, h5_pool)
	// save flag if h5 is open
	bool is_open = (t.file_id != 0);
	ar << is_open;
	ar << t.fname;
	// also save base name
	ar << (const std::string&)basename(t.fname);

	// we should only dump group names
	typedef h5_pool::map_hid_t::const_iterator g_iterator;
	ar << (const std::size_t&)t.group_id.size();
	for(g_iterator pg = t.group_id.begin(); pg != t.group_id.end(); ++pg)
		ar << pg->first;

	// dumb nodes
	typedef h5_pool::map_t::const_iterator m_iterator;
	ar << (const std::size_t&)t.h5_map.size();
	std::string g_name;
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
			ar << (const H5T_class_t&)H5Tget_class(dtype);
			ar << (const hid_t&)H5Tget_size(dtype);
		}
	}

	// other data
	ar << t.pool_dims;
	//ar << boser::make_array(&t.pool_dims, 3);
	ar << t.n_pool_dims;
	// flush buffers
	const_cast< h5_pool& >(t).flush();
BLUE_SKY_CLASS_SRZ_FCN_END

BLUE_SKY_CLASS_SRZ_FCN_BEGIN(load, h5_pool)
	bool do_open;
	ar >> do_open;
	ar >> t.fname;

	// load base name
	std::string h5_basename;
	ar >> h5_basename;

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

	// build groups map
	std::size_t g_size;
	std::string g_name;
	hid_t g_id;
	ar >> g_size;
	for(std::size_t i = 0; i < g_size; ++i) {
		ar >> g_name;
		g_id = -1;
		if(do_open) {
			g_id = demand_h5_group(t.file_id, g_name);
			t.fill_map(g_id);
		}
		// add group to map
		t.group_id[g_name] = g_id;
	}

	// read and build h5_map
	ar >> g_size;
	h5_pair p;
	std::string d_name;
	// type info
	H5T_class_t type_class;
	size_t type_size;
	// fill h5_map
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
			ar >> type_class;
			ar >> type_size;
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
			// find corresponding group
			g_id = t.group_id[g_name];
			if(g_id < 0)
				g_id = demand_h5_group(t.file_id, g_name);

			// try to open dataset
			dset_tmp = H5Dopen(g_id, d_name.c_str());
			if(dset_tmp >= 0) {
				p.dset = dset_tmp;
				// open other attributes
				p.dspace = H5Dget_space(dset_tmp);
				if(p.dspace < 0) {
					bs_throw_exception(
						boost::format("Can't get dataspace for dataset %s in group %d")
						% d_name % g_name
					);
				}

				p.dtype = H5Dget_type(dset_tmp);
				if(p.dtype < 0) {
					bs_throw_exception(
						boost::format("Can't get datatype for dataset %s in group %d")
						% d_name % g_name
					);
				}

				p.plist = H5Dget_create_plist(dset_tmp);
				if(p.plist < 0) {
					bs_throw_exception(
						boost::format("Can't get property list for dataset %s in group %d")
						% d_name % g_name
					);
				}
			}
		}

		// if hids cannot be loaded from open dataset
		// try to create plist & dtype manually
		if(dset_tmp < 0) {
			p.plist = H5Pcreate(H5P_DATASET_CREATE);
			if(p.plist < 0) {
				bs_throw_exception(
					boost::format("Can't create property for dataset %s in group %d")
					% d_name % g_name
				);
			}

			// create type
			if(
				type_class == H5T_COMPOUND ||
				type_class == H5T_OPAQUE ||
				type_class == H5T_ENUM ||
				type_class == H5T_STRING
			)
				p.dtype = H5Tcreate(type_class, type_size);
			else
				p.dtype = H5Tcopy(p.dtype);
			if(p.dtype < 0) {
				bs_throw_exception(
					boost::format("Can't create datatype of class %d and size %d")
					% type_class % type_size
				);
			}

			// if during saving dataset was open but we failed to open it now
			// try to create it
			if(p.dset >= 0 && g_name.size()) {
				// try to create dataspace
				p.dspace = H5Screate_simple(p.n_dims, p.h5_dims, NULL);
				if(p.dspace < 0) {
					bs_throw_exception(
						boost::format("Can't create simple dataspace for dataset %s in group %d")
						% d_name % g_name
					);
				}

				// and finally dataset
				p.dset = H5Dcreate(g_id, d_name.c_str(), p.dtype, p.dspace, p.plist);
				if(p.dset < 0) {
					bs_throw_exception(
						boost::format("Can't create dataset %s in group %d")
						% d_name % g_name
					);
				}
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
BLUE_SKY_TYPE_SERIALIZE_IMPL(h5_pool)

