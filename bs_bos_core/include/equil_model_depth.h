#ifndef EQUIL_MODEL_DEPTH_H
#define EQUIL_MODEL_DEPTH_H

#include "data_class.h"
#include "pvt_3p_iface.h"
#include "scal_3p_iface.hpp"
#include "pvt_3p_iface.h"
#include "scal_3p.h"
//Needed for scal_3p.h
#include "scale_array_holder.h"
#include "scal_region_info.h"
#include "scal_region.h"
#include "scal_2p_data_holder.h"
//--------------------
#include "scal_dummy_iface.h"
#include "pvt_dummy_iface.h"

namespace blue_sky
{
	typedef boost::array <t_long, FI_PHASE_TOT>			phase_d_t;
    typedef boost::array <t_long, FI_PHASE_TOT>			sat_d_t;
	typedef BS_SP(jfunction)							sp_jfunction;

	void
	calc_equil_vs_depth (bool is_o, bool is_g, bool is_w,
						 const phase_d_t &phase_d, const sat_d_t &sat_d,
						 t_long n_phases, t_long phases, t_long sat_counter,
						 t_long n_depth,
                         const t_long n_eql,               //!< number of equil regions      
                         const stdv_long &sat_regions,      //!< (n_eql) scal region number for current equil region 
                         const stdv_long &pvt_regions,      //!< (n_eql) pvt region number for current equil region 
                         BS_SP (pvt_3p_iface) pvt_prop,    //!< PVT properties
                         BS_SP (scal_3p_iface) scal_prop,  //!< SCAL properties
                         spv_float equil,                 //!< main equil data (from EQUIL keyword)
                         idata::vval_vs_depth *rsvd,
                         idata::vval_vs_depth *pbvd,
                         const stdv_float &min_depth,     //!< n_eql
                         const stdv_float &max_depth,     //!< n_eql
                         const stdv_float &perm,          //!< (n_eql) permeability
                         const stdv_float &poro,          //!< (n_eql) poro
                         stdv_double &pressure,          //!< (n_phases * n_eql * n_depth)
                         stdv_double &saturation         //!< (n_phases * n_eql * n_depth)
                        );

	class BS_API_PLUGIN equil_model_depth : public objbase
	{
	public:
		void
		py_calc_equil(bool is_o, bool is_g, bool is_w,
		              BS_SP(scal_3p_iface) scal_props, 
					  BS_SP(pvt_3p_iface) pvt_props,
					  spv_float equil,
					  const stdv_float &min_depth,
					  const stdv_float &max_depth,
					  const stdv_float &perm,
					  const stdv_float &poro,
					  t_long n_depth,
					  //stdv_double density,
					  sp_jfunction jfunc_water,
					  sp_jfunction jfunc_oil);

		spv_double get_pressure ();
		spv_double get_saturation ();

	public:
		BLUE_SKY_TYPE_DECL (equil_model_depth);
	private:
		spv_double pressure;
		spv_double saturation;
	};
}

#endif
