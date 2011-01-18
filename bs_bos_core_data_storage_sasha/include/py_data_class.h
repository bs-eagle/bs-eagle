#ifndef PY_DATA_CLASS_H
#define PY_DATA_CLASS_H

#include "data_class.h"
#include "shared_array.h"

namespace blue_sky { namespace python {

template<class strategy_t>
class BS_API_PLUGIN py_idata : public py_bs_node
{
public:
	typedef idata<strategy_t>             wrapped_t;
	typedef smart_ptr<wrapped_t, true>    sp_idata_t;

	typedef typename strategy_t::fp_type_t   item_t;
	typedef typename strategy_t::i_type_t  index_t;
	typedef typename wrapped_t::vec_i     vec_i;

	py_idata();
	py_idata(const sp_idata_t&);
	py_idata(const py_idata&);

	void init();

	int get_rpo_model () const;
	void set_rpo_model (int t);

	double get_minimal_pore_volume () const;
	void set_minimal_pore_volume (double t);

	int get_restart () const;
	void set_restart (int t);

	index_t get_nx () const;
	index_t get_ny () const;
	index_t get_nz () const;
	void set_nx (index_t t);
	void set_ny (index_t t);
	void set_nz (index_t t);

	index_t get_pvt_region () const;
	index_t get_sat_region () const;
	index_t get_eql_region () const;
	index_t get_fip_region () const;
	void set_pvt_region (index_t t);
	void set_sat_region (index_t t);
	void set_eql_region (index_t t);
	void set_fip_region (index_t t);

	int get_fi_phases () const;
	index_t get_fi_n_phase () const;
	index_t get_rock_region () const;
	void set_fi_phases (int t);
	void set_fi_n_phase (index_t t);
	void set_rock_region (index_t t);

	const vec_i &get_equil_regions () const;

	int get_units_in () const;
	void set_units_in (int t);
	int get_units_out () const;
	void set_units_out (int t);

	std::string get_title () const;
	void set_title (const std::string &t);

	array_uchar_t   get_int_array (const std::string& array_name);
	array_float_t get_float_array (const std::string& array_name);

	void set_int_array (const std::string& array_name, const boost::python::object &obj);
	void set_float_array (const std::string& array_name, const boost::python::object &obj);

};

void py_export_idata();

} }	// namesapce blue_sky::python

#endif // PY_DATA_CLASS_H

