/**
	\file data_class.cpp
	\brief implimenataion of idata class methods
	\author Nikonov Maxim
*/
#include "bs_bos_core_data_storage_stdafx.h"

#include "data_class.h"
#include "arrays.h"
#include "constants.h"

//! Default minimal pore volume allowed for active cells
//#define DEFAULT_MINIMAL_PORE_VOLUME     1.e-1
#define DEFAULT_MINIMAL_PORE_VOLUME     1.e-7
//! Default minimal pore volume allowed for active cells to splice with other cells
// miryanov: set MINSV equal to MINPV for SPE10 model
#define DEFAULT_MINIMAL_SPLICE_VOLUME     1.e-7
//! Default maximum thickness allowed between active cells to be coupled
#define DEFAULT_MAXIMUM_SPLICE_THICKNESS     1.e-2

using namespace std;

namespace blue_sky {
using namespace pool;

template< class strategy_t >
struct idata< strategy_t >::idata_traits : public bs_node::restrict_types {
	const char *sort_name() const {
		return "idata trait";
	}

	bool accepts(const sp_link &l) const {
		if (l->name().find("_in_pool",0,l->name().size()))
			return smart_ptr< this_t >(l->data(), bs_dynamic_cast());
		return false;
	}
};

template <class strategy_t>
idata<strategy_t>::~idata ()
{

}

template< class strategy_t >
idata< strategy_t >::idata(bs_type_ctor_param param)
{
	init_section = 0;
	nx=0;
	ny=0;
	nz=0;
	pvt_region = 1;
	sat_region = 1;
	eql_region = 1;
	fip_region = 1;
	rock_region = 0;
	minimal_pore_volume = DEFAULT_MINIMAL_PORE_VOLUME;
	minimal_splice_volume = DEFAULT_MINIMAL_SPLICE_VOLUME;
	maximum_splice_thickness = DEFAULT_MAXIMUM_SPLICE_THICKNESS;
	fi_n_phase = 0;
	fi_phases = 0;
	rpo_model = 0;//RPO_DEFAULT_MODEL;
	build_operator_list();
	init();
}

//template <class strategy_t>
//idata<strategy_t>::idata(const this_t &src)
//  :bs_pool_node(src)
//      d_map(give_kernel::Instance().create_object_copy(src.d_map))
//{
//	*this = src;
//}

template <class strategy_t>
void idata<strategy_t>::init()
{
	//depth.resize((nx+1) * (ny+1) * (nz+1));
	//ahelper.init_names_maps ();
}

//template <class strategy_t>
//idata<strategy_t> &idata<strategy_t>::operator=(const this_t &src)
//{
//	i_map = src.i_map;
//	d_map = src.d_map;
//
//	ahelper = src.ahelper;
//	return *this;
//}

template< class strategy_t >
typename idata< strategy_t >::vval_vs_depth &idata< strategy_t >::get_prvd() {
	return prvd;
}

template< class strategy_t >
typename idata< strategy_t >::vval_vs_depth &idata< strategy_t >::get_rsvd() {
	return rsvd;
}

template< class strategy_t >
typename idata< strategy_t >::vval_vs_depth &idata< strategy_t >::get_pbvd() {
	return pbvd;
}

/*!
\brief updating dx, dy, dz
			 add values in nonactive blocks for mesh generation algorithm
\return if success  0
*/

template< class strategy_t >
void idata< strategy_t >::set_defaults_in_pool() {
//	using bs_pool_subnode::array_size::dimens_idx;

	//create_array(array_name_fp[MULTX], at(array_sizes_fp(), "MULTX"),  at(array_initv_fp(), "MULTX"));
	//create_array(array_name_fp[MULTY], at(array_sizes_fp(), "MULTY"),  at(array_initv_fp(), "MULTY"));
	//create_array(array_name_fp[MULTZ], at(array_sizes_fp(), "MULTZ"),  at(array_initv_fp(), "MULTZ"));
	//create_array(array_name_fp[NTG], at(array_sizes_fp(), "NTG"), at(array_initv_fp(), "NTG"));
	//create_array(array_name_fp[MULTPV], at(array_sizes_fp(), "MULTPV"), at(array_initv_fp(), "MULTPV"));
	create_array(array_name_i[ACTNUM], at(array_sizes_i(), "ACTNUM"), at(array_initv_i(), "ACTNUM"));
}


template< class strategy_t >
void idata< strategy_t >::set_region (int r_pvt,int r_sat, int r_eql, int r_fip) {
	this->pvt_region = r_pvt;
	this->fip_region = r_fip;
	this->sat_region = r_sat;
	this->eql_region = r_eql;

	// check
	if (r_pvt <= 0 || r_sat <= 0 || r_eql <= 0 || r_fip <= 0)
	  {
		bs_throw_exception ("One of init parameters <= 0");
	  }

	assign (this->rock, r_pvt, -1);
	assign (this->p_ref, r_pvt, -1);

	this->equil.resize(EQUIL_TOTAL * eql_region); //!TODO: EQUIL_TOTAL instead of 3

	this->pvto.resize(r_pvt);
	this->pvtdo.resize(r_pvt);
	this->pvtg.resize(r_pvt);
	this->pvtw.resize(r_pvt);
}

template< class strategy_t >
void idata< strategy_t >::set_density (const std::vector <item_t> &density) {
	if ((density.size() % 3 != 0) || (density.size()<3))
	{
		bs_throw_exception ("Not enough valid arguments yet");
	}

	set_density_internal(&density[0]);
}

template <class strategy_t>
void idata<strategy_t>::set_density_internal (const item_t *density) {
	std::ostringstream out_s;

	if (!pvto.size())
		throw bs_exception("idata.set_density()","PVT table for oil has not been initialized yet");

	if (!pvtw.size())
		throw bs_exception("idata.set_density()","PVT table for water has not been initialized yet");

	if (!pvtg.size())
		throw bs_exception("idata.set_density()","PVT table for gas has not been initialized yet");

	if (!rock.size())
		throw bs_exception("idata.set_density()","Rock properties table has not been initialized yet");

	if (!equil.size())
		throw bs_exception("idata.set_density()","EQUIL table has not been initialized yet");

	for (index_t i = 0; i < pvt_region; ++i) {
		typename idata<strategy_t>::pvt_info &pvto__ = pvto[i];
		typename idata<strategy_t>::pvt_info &pvtw__ = pvtw[i];
		typename idata<strategy_t>::pvt_info &pvtg__ = pvtg[i];
		typename idata<strategy_t>::pvt_info &pvdo__ = pvtdo[i];

		pvto__.has_density_   = true;
		pvto__.density_       = density[i*3];
		pvto__.molar_density_ = density[i*3];

		pvdo__.has_density_   = true;
		pvdo__.density_       = density[i*3];
		pvdo__.molar_density_ = density[i*3];

		pvtw__.has_density_   = true;
		pvtw__.density_       = density[i*3+1];
		pvtw__.molar_density_ = density[i*3+1];

		pvtg__.has_density_   = true;
		pvtg__.density_       = density[i*3+2];
		pvtg__.molar_density_ = density[i*3+2];
		}
}
  
/* 
		ARRAYS_HELPER methods
*/

template< typename strategy_t >
bs_array_i& idata< strategy_t >::get_int_array (const string& array_name) {
	//if (!i_pool()->contain (array_name)) {
	//	bs_throw_exception ("Array not initialized yet");
	//}
	return get_non_empty_i(array_name);
}

template < typename strategy_t >
bs_array_fp& idata< strategy_t >::get_float_array (const string& array_name) {
	//if (!d_pool()->contain (array_name)) {
	//	bs_throw_exception ("Array not initialized yet");
	//}
	return get_non_empty_fp(array_name);
}

// create object
BLUE_SKY_TYPE_STD_CREATE_T_DEF (idata,(class));
BLUE_SKY_TYPE_STD_COPY_T_DEF (idata,(class));


// array map implementation
BLUE_SKY_TYPE_IMPL_T_EXT (1, (idata<base_strategy_fif>), 1, (bs_pool_node), "idata float", "Initial data storage", "", false);
BLUE_SKY_TYPE_IMPL_T_EXT (1, (idata<base_strategy_did>), 1, (bs_pool_node), "idata double", "Initial data storage", "", false);
BLUE_SKY_TYPE_IMPL_T_EXT (1, (idata<base_strategy_dif>), 1, (bs_pool_node), "idata mixi", "Initial data storage", "", false);

/*
\brief register plugin
*/
bool register_idata (const plugin_descriptor& pd) {
	bool res;
	res = BLUE_SKY_REGISTER_TYPE (pd, idata<base_strategy_fif>); 
	res &= BLUE_SKY_REGISTER_TYPE (pd, idata<base_strategy_did>);
	res &= BLUE_SKY_REGISTER_TYPE (pd, idata<base_strategy_dif>);
	return res;
}

}	// namespace blue_sky

