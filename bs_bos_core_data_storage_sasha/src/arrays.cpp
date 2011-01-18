#include "arrays.h"
#include "pool_treeish.h"
#include <boost/assign/list_of.hpp> // for 'map_list_of()'

namespace blue_sky { namespace pool {
using namespace std;
using namespace boost::assign;

const char* array_name_i[] = {
	"MPFANUM",
	"EQLNUM",
	"SATNUM",
	"PVTNUM",
	"ACTNUM",
	"FIPNUM",
	"BNDNUM",
	"EOSNUM",
	"ROCKNU"
};

const char* array_name_fp[] = {
	"PERMXY",
	"PERMYZ",
	"PERMZX",
	"PERMX",
	"PERMY",
	"PERMZ",
	"PORO",
	"NTG",
	"SGL",
	"SGU",
	"SOGCR",
	"SGCR",
	"SWL",
	"SWU",
	"SOWCR",
	"SWCR",
	"PBUB",
	"RS",
	"PCW",
	"SWATINIT",
	"SOIL",
	"SWAT",
	"PRESSURE",
	"DX",
	"DY",
	"DZ",
	"MULTX",
	"MULTY",
	"MULTZ",
	"TOPS",
	"MULTPV",
	"SGAS"
};

const char* array_name_comp[] = {
	"ACF",
	"MW",
	"PCRIT",
	"TCRIT",
	"VCRIT",
	"ZCRIT",
	"RTEMP",
	"ACFS",
	"MWS",
	"PCRITS",
	"TCRITS",
	"VCRITS",
	"ZCRITS",
	"XMF",
	"YMF",
	"ZMF",
	"ZI"
};

// make instance of array_initv_fp
BS_API_PLUGIN const map< string, float >& array_initv_fp() {
	static map< string, float > data = map_list_of
		("PERMXY"   ,0)
		("PERMYZ"   ,0)
		("PERMZX"   ,0)
		("PERMX"    ,0)
		("PERMY"    ,0)
		("PERMZ"    ,0)
		("PORO"     ,0)
		("NTG"      ,1)
		("SGL"      ,0)
		("SGU"      ,0)
		("SOGCR"    ,0)
		("SGCR"     ,0)
		("SWL"      ,0)
		("SWU"      ,0)
		("SOWCR"    ,0)
		("SWCR"     ,0)
		("PBUB"     ,0)
		("RS"       ,0)
		("PCW"      ,0)
		("SWATINIT" ,0)
		("SOIL"     ,0)
		("SWAT"     ,0)
		("PRESSURE" ,0)
		("DX"       ,0)
		("DY"       ,0)
		("DZ"       ,0)
		("MULTX"    ,1)
		("MULTY"    ,1)
		("MULTZ"    ,1)
		("TOPS"     ,0)
		("MULTPV"   ,1)
		("SGAS"     ,0)
		;
	return data;
}

BS_API_PLUGIN const map< string, int >& array_initv_i() {
	static map< string, int > data = map_list_of
		("MPFANUM" ,0)
		("EQLNUM"  ,0)
		("SATNUM"  ,0)
		("PVTNUM"  ,0)
		("ACTNUM"  ,1)
		("FIPNUM"  ,0)
		("BNDNUM"  ,0)
		("EOSNUM"  ,0)
		("ROCKNU"  ,0)
		;
	return data;
}

BS_API_PLUGIN const map< string, float >& array_initv_c() {
	static map< string, float > data = map_list_of
		("ACF"    ,0)
		("MW"     ,0)
		("PCRIT"  ,0)
		("TCRIT"  ,0)
		("VCRIT"  ,0)
		("ZCRIT"  ,0)
		("RTEMP"  ,0)
		("ACFS"   ,0)
		("MWS"    ,0)
		("PCRITS" ,0)
		("TCRITS" ,0)
		("VCRITS" ,0)
		("ZCRITS" ,0)
		("XMF"    ,0)
		("YMF"    ,0)
		("ZMF"    ,0)
		("ZI"     ,0)
		;
	return data;
}

BS_API_PLUGIN const map< string, array_size >& array_sizes_fp() {
	static map< string, array_size > data = map_list_of
		("PERMXY"   ,array_size(1, 1, 1, 0, 0, 0))
		("PERMYZ"   ,array_size(1, 1, 1, 0, 0, 0))
		("PERMZX"   ,array_size(1, 1, 1, 0, 0, 0))
		("PERMX"    ,array_size(1, 1, 1, 0, 0, 0))
		("PERMY"    ,array_size(1, 1, 1, 0, 0, 0))
		("PERMZ"    ,array_size(1, 1, 1, 0, 0, 0))
		("PORO"     ,array_size(1, 1, 1, 0, 0, 0))
		("NTG"      ,array_size(1, 1, 1, 0, 0, 0))
		("SGL"      ,array_size(1, 1, 1, 0, 0, 0))
		("SGU"      ,array_size(1, 1, 1, 0, 0, 0))
		("SOGCR"    ,array_size(1, 1, 1, 0, 0, 0))
		("SGCR"     ,array_size(1, 1, 1, 0, 0, 0))
		("SWL"      ,array_size(1, 1, 1, 0, 0, 0))
		("SWU"      ,array_size(1, 1, 1, 0, 0, 0))
		("SOWCR"    ,array_size(1, 1, 1, 0, 0, 0))
		("SWCR"     ,array_size(1, 1, 1, 0, 0, 0))
		("PBUB"     ,array_size(1, 1, 1, 0, 0, 0))
		("RS"       ,array_size(1, 1, 1, 0, 0, 0))
		("PCW"      ,array_size(1, 1, 1, 0, 0, 0))
		("SWATINIT" ,array_size(1, 1, 1, 0, 0, 0))
		("SOIL"     ,array_size(1, 1, 1, 0, 0, 0))
		("SWAT"     ,array_size(1, 1, 1, 0, 0, 0))
		("PRESSURE" ,array_size(1, 1, 1, 0, 0, 0))
		("DX"       ,array_size(1, 1, 1, 0, 0, 0))
		("DY"       ,array_size(1, 1, 1, 0, 0, 0))
		("DZ"       ,array_size(1, 1, 1, 0, 0, 0))
		("MULTX"    ,array_size(1, 1, 1, 0, 0, 0))
		("MULTY"    ,array_size(1, 1, 1, 0, 0, 0))
		("MULTZ"    ,array_size(1, 1, 1, 0, 0, 0))
		("TOPS"     ,array_size(1, 1, 0, 0, 0, 1))
		("MULTPV"   ,array_size(1, 1, 1, 0, 0, 0))
		("SGAS"     ,array_size(1, 1, 1, 0, 0, 0))
		;
	return data;
}

BS_API_PLUGIN const map< string, array_size >& array_sizes_i() {
	static map< string, array_size > data = map_list_of
		("MPFANUM" ,array_size(1, 1, 1, 0, 0, 0))
		("EQLNUM"  ,array_size(1, 1, 1, 0, 0, 0))
		("SATNUM"  ,array_size(1, 1, 1, 0, 0, 0))
		("PVTNUM"  ,array_size(1, 1, 1, 0, 0, 0))
		("ACTNUM"  ,array_size(1, 1, 1, 0, 0, 0))
		("FIPNUM"  ,array_size(1, 1, 1, 0, 0, 0))
		("BNDNUM"  ,array_size(1, 1, 1, 0, 0, 0))
		("EOSNUM"  ,array_size(1, 1, 1, 0, 0, 0))
		("ROCKNU"  ,array_size(1, 1, 1, 0, 0, 0))
		;
	return data;
}

} } 	// namespace blue_sky { namespace pool

