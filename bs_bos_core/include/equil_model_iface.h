#ifndef __EQUIL_MODEL_IFACE_H
#define __EQUIL_MODEL_IFACE_H


namespace blue_sky
{
  class jfunction;
  class table_iface;
  class pvt_3p_iface;
  class scal_3p_iface;

	class BS_API_PLUGIN equil_model_iface : public objbase
	{
	public:

  	typedef BS_SP(jfunction)							sp_jfunction;
	  typedef std::vector <spv_double>           stdv_spv_double;
	  typedef BS_SP (table_iface)                sp_table_t;
	  typedef std::vector <BS_SP (table_iface)>  sp_table_array_t;

    virtual ~equil_model_iface () {}

		virtual void
		py_calc_equil(bool is_o, bool is_g, bool is_w,
		              BS_SP(scal_3p_iface) scal_props,
					  BS_SP(pvt_3p_iface) pvt_props,
					  BS_SP (table_iface) equil,
					  const stdv_float &min_depth,
					  const stdv_float &max_depth,
					  const stdv_float &perm,
					  const stdv_float &poro,
					  t_long n_depth,
					  //stdv_double density,
					  sp_jfunction jfunc_water,
					  sp_jfunction jfunc_oil) = 0;

		virtual
		spv_double get_pressure () = 0;

		virtual
		spv_double get_saturation () = 0;

    virtual void init_equil_model (const t_long n_equil_regions_,
                                   t_int n_phases) = 0;

    virtual BS_SP (table_iface)
    get_equil_region_data (const t_long region) const = 0;

    virtual std::list <BS_SP( table_iface)>
    get_equil_data () const = 0;

    virtual BS_SP (table_iface)
    get_pbvd_region_data (const t_long region) const = 0;

    virtual std::list <BS_SP( table_iface)>
    get_pbvd_data () const = 0;

    virtual BS_SP (table_iface)
    get_rsvd_region_data (const t_long region) const = 0;

    virtual std::list <BS_SP( table_iface)>
    get_rsvd_data () const = 0;

    virtual t_long
    get_n_equil_regions () const = 0;
	};
}


#endif // __EQUIL_MODEL_IFACE_H
