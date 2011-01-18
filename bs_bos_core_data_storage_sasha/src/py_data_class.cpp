#include "bs_bos_core_data_storage_stdafx.h"

#include "py_data_class.h"
#include "numpy_tools.h"
#include "py_object_handler.h"

#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>

using namespace boost::python;
namespace bp = boost::python;

using namespace blue_sky::pool;

namespace blue_sky {

typedef float_t float_t;
typedef uchar_t uchar_t;

namespace python {

template < class strategy_t >
py_idata<strategy_t>::py_idata()
: py_bs_node(give_kernel::Instance().create_object(wrapped_t::bs_type()))
{}

template < class strategy_t >
py_idata<strategy_t>::py_idata(const sp_idata_t &src)
: py_bs_node(sp_obj(src))
{}

template < class strategy_t >
py_idata<strategy_t>::py_idata(const py_idata<strategy_t> &src)
: py_bs_node(src.sp)
{}

template < class strategy_t >
void py_idata<strategy_t>::init()
{
}

template < class strategy_t >
int py_idata<strategy_t>::get_rpo_model () const
{
	return this->get_spx<wrapped_t> ()->rpo_model;
}

template < class strategy_t >
void py_idata<strategy_t>::set_rpo_model (int t)
{
	this->get_spx<wrapped_t> ()->rpo_model = t;
}


template < class strategy_t >
double py_idata<strategy_t>::get_minimal_pore_volume () const
{
	return this->get_spx<wrapped_t> ()->minimal_pore_volume;
}

template < class strategy_t >
void py_idata<strategy_t>::set_minimal_pore_volume (double t)
{
	this->get_spx<wrapped_t> ()->minimal_pore_volume = t;
}


template < class strategy_t >
int py_idata<strategy_t>::get_restart () const
{
	return this->get_spx<wrapped_t> ()->restart;
}

	template < class strategy_t >
void py_idata<strategy_t>::set_restart (int t)
{
	this->get_spx<wrapped_t> ()->restart = t;
}

template < class strategy_t >
typename py_idata<strategy_t>::index_t py_idata<strategy_t>::get_nx () const
{
	return this->get_spx<wrapped_t> ()->nx;
}

template < class strategy_t >
void py_idata<strategy_t>::set_nx (index_t t)
{
	this->get_spx<wrapped_t> ()->nx = t;
}

template < class strategy_t >
typename py_idata<strategy_t>::index_t py_idata<strategy_t>::get_ny () const
{
	return this->get_spx<wrapped_t> ()->ny;
}

	template < class strategy_t >
void py_idata<strategy_t>::set_ny (index_t t)
{
	this->get_spx<wrapped_t> ()->ny = t;
}

template < class strategy_t >
typename py_idata<strategy_t>::index_t py_idata<strategy_t>::get_nz () const
{
	return this->get_spx<wrapped_t> ()->nz;
}

template < class strategy_t >
void py_idata<strategy_t>::set_nz (index_t t)
{
	this->get_spx<wrapped_t> ()->nz = t;
}


template < class strategy_t >
typename py_idata<strategy_t>::index_t py_idata<strategy_t>::get_pvt_region () const
{
	return this->get_spx<wrapped_t> ()->pvt_region;
}

template < class strategy_t >
void py_idata<strategy_t>::set_pvt_region (index_t t)
{
	this->get_spx<wrapped_t> ()->pvt_region = t;
}

template < class strategy_t >
typename py_idata<strategy_t>::index_t py_idata<strategy_t>::get_sat_region () const
{
	return this->get_spx<wrapped_t> ()->sat_region;
}

template < class strategy_t >
void py_idata<strategy_t>::set_sat_region (index_t t)
{
	this->get_spx<wrapped_t> ()->sat_region = t;
}

template < class strategy_t >
typename py_idata<strategy_t>::index_t py_idata<strategy_t>::get_eql_region () const
{
	return this->get_spx<wrapped_t> ()->eql_region;
}

template < class strategy_t >
void py_idata<strategy_t>::set_eql_region (index_t t)
{
	this->get_spx<wrapped_t> ()->eql_region = t;
}

template < class strategy_t >
typename py_idata<strategy_t>::index_t py_idata<strategy_t>::get_fip_region () const
{
	return this->get_spx<wrapped_t> ()->fip_region;
}

template < class strategy_t >
void py_idata<strategy_t>::set_fip_region (index_t t)
{
	this->get_spx<wrapped_t> ()->fip_region = t;
}


template < class strategy_t >
typename py_idata<strategy_t>::index_t py_idata<strategy_t>::get_fi_n_phase () const
{
	return this->get_spx<wrapped_t> ()->fi_n_phase;
}

template < class strategy_t >
void py_idata<strategy_t>::set_fi_n_phase (index_t t)
{
	this->get_spx<wrapped_t> ()->fi_n_phase = t;
}


template < class strategy_t >
int py_idata<strategy_t>::get_fi_phases () const
{
	return this->get_spx<wrapped_t> ()->fi_phases;
}

template < class strategy_t >
void py_idata<strategy_t>::set_fi_phases (int t)
{
	this->get_spx<wrapped_t> ()->fi_phases = t;
}

template < class strategy_t >
typename py_idata<strategy_t>::index_t py_idata<strategy_t>::get_rock_region () const
{
	return this->get_spx<wrapped_t> ()->rock_region;
}

template < class strategy_t >
void py_idata<strategy_t>::set_rock_region (index_t t)
{
	this->get_spx<wrapped_t> ()->rock_region = t;
}

template < class strategy_t >
const typename py_idata<strategy_t>::vec_i &
py_idata<strategy_t>::get_equil_regions () const
{
	return this->get_spx<wrapped_t> ()->equil_regions;
}

template < class strategy_t >
int py_idata<strategy_t>::get_units_in () const
{
	return this->get_spx<wrapped_t> ()->units_in;
}

template < class strategy_t >
void py_idata<strategy_t>::set_units_in (int t)
{
	this->get_spx<wrapped_t> ()->units_in = t;
}

template < class strategy_t >
int py_idata<strategy_t>::get_units_out () const
{
	return this->get_spx<wrapped_t> ()->units_out;
}

template < class strategy_t >
void py_idata<strategy_t>::set_units_out (int t)
{
	this->get_spx<wrapped_t> ()->units_out = t;
}

template < class strategy_t >
std::string py_idata<strategy_t>::get_title () const
{
	return this->get_spx<wrapped_t> ()->title;
}

template < class strategy_t >
void py_idata<strategy_t>::set_title (const std::string &t)
{
	this->get_spx<wrapped_t> ()->title = t;
}

template <typename strategy_t>
array_uchar_t py_idata <strategy_t>::get_int_array (const std::string& array_name)
{
	bs_array_i& arr = this->get_spx(this)->get_int_array(array_name);
	return array_uchar_t(&arr[0], arr.size());
}

template <typename strategy_t>
array_float_t py_idata <strategy_t>::get_float_array (const std::string& array_name)
{
	bs_array_fp& arr = this->get_spx(this)->get_float_array(array_name);
	return array_float_t(&arr[0], arr.size());
}

template <typename strategy_t>
void py_idata <strategy_t>::set_int_array (const std::string& array_name, const boost::python::object &obj)
{
	bp::ssize_t len = numpy_tools::get_len (obj);
	uchar_t *buffer = numpy_tools::get_buffer <uchar_t> (obj, len);

	bs_array_i& arr = this->get_spx(this)->get_int_array(array_name);
	if(arr.size() == len)
		std::copy(buffer, buffer + len, arr.begin());

	//(*this->get_spx (this)->i_map)[cur_index].array = array_uchar_t (buffer, len);
	//this->get_spx (this)->subscribe (objbase::on_delete, new tools::py_object_handler (obj.ptr ()));

}

template <typename strategy_t>
void py_idata <strategy_t>::set_float_array (const std::string& array_name, const boost::python::object &obj)
{
	bp::ssize_t len   = numpy_tools::get_len (obj);
	float_t *buffer = numpy_tools::get_buffer <float_t> (obj, len);

	bs_array_fp& arr = this->get_spx(this)->get_float_array(array_name);
	if(arr.size() == len)
		std::copy(buffer, buffer + len, arr.begin());

	//(*this->get_spx (this)->d_map)[cur_index].array = array_float_t (buffer, len);
	//this->get_spx (this)->subscribe (objbase::on_delete, new tools::py_object_handler (obj.ptr ()));
}

template < typename strategy_t, typename strategy_t::i_type_t n >
struct idata_int_array_getter
{
	array_uchar_t get_array (py_idata <strategy_t> *self) {
		return self->get_int_array(array_name_i[n]);
	}

	void set_array (py_idata <strategy_t> *self, const boost::python::object &obj) {
		self->set_int_array (array_name_i[n], obj);
	}
};

template <typename strategy_t, typename strategy_t::i_type_t n >
struct idata_float_array_getter
{
	array_float_t get_array (py_idata <strategy_t> *self) {
		return self->get_float_array (array_name_fp[n]);
	}

	void set_array (py_idata <strategy_t> *self, const boost::python::object &obj) {
		self->set_float_array (array_name_fp[n], obj);
	}

	std::string name_;
};

#define REGISTER_PROPERTIES_I(z, n_, t_)                                              \
	BOOST_PP_TUPLE_ELEM (4, 0, t_).add_property (BOOST_PP_TUPLE_ELEM (4, 1, t_)[n_],    \
			&BOOST_PP_TUPLE_ELEM (4, 2, t_)< BOOST_PP_TUPLE_ELEM (4, 3, t_), n_ >::get_array,  \
			&BOOST_PP_TUPLE_ELEM (4, 2, t_)< BOOST_PP_TUPLE_ELEM (4, 3, t_), n_ >::set_array);

#define REGISTER_PROPERTIES(from_, to_, t_) \
	BOOST_PP_REPEAT_FROM_TO (from_, to_, REGISTER_PROPERTIES_I, t_)

	template <class strategy_t>
void py_export_idata_t(const char *name)
{
	typedef py_idata <strategy_t> T;

	class_<T, bases<py_bs_node> > idata_export(name);

	idata_export
		.def("init",&T::init)
		.add_property ("rpo_model",           &T::get_rpo_model, &T::set_rpo_model)
		.add_property ("minimal_pore_volume", &T::get_minimal_pore_volume, &T::set_minimal_pore_volume)
		.add_property ("restart",             &T::get_restart, &T::set_restart)
		.add_property ("nx",                  &T::get_nx, &T::set_nx)
		.add_property ("ny",                  &T::get_ny, &T::set_ny)
		.add_property ("nz",                  &T::get_nz, &T::set_nz)
		.add_property ("pvt_region",          &T::get_pvt_region, &T::set_pvt_region)
		.add_property ("sat_region",          &T::get_sat_region, &T::set_sat_region)
		.add_property ("eql_region",          &T::get_eql_region, &T::set_eql_region)
		.add_property ("fip_region",          &T::get_fip_region, &T::set_fip_region)
		.add_property ("fi_n_phase",          &T::get_fi_n_phase, &T::set_fi_n_phase)
		.add_property ("fi_phases",           &T::get_fi_phases, &T::set_fi_phases)
		.add_property ("rock_region",         &T::get_rock_region, &T::set_rock_region)
		.add_property ("equil_regions",       make_function (&T::get_equil_regions, return_value_policy<copy_const_reference> ()))
		.add_property ("units_in",            &T::get_units_in, &T::set_units_in)
		.add_property ("units_out",           &T::get_units_out, &T::set_units_out)
		.add_property ("title",               &T::get_title, &T::set_title)
		.add_property ("int_array",           &T::get_int_array, &T::set_int_array)
		.add_property ("float_array",         &T::get_float_array, &T::set_float_array)
		.def ("get_int_array",                &T::get_int_array)
		.def ("get_float_array",              &T::get_float_array)
		.def ("set_int_array",                &T::set_int_array)
		.def ("set_float_array",              &T::set_float_array)
		;

	REGISTER_PROPERTIES (0, 9,  (idata_export, array_name_i,    idata_int_array_getter,   strategy_t));
	REGISTER_PROPERTIES (0, 16, (idata_export, array_name_fp, idata_float_array_getter, strategy_t))
}

template <typename T>
long get_index (T *c, typename T::size_type i, const char *context)
{
	long index = long (i);
	if (index < 0)
	{
          index += (long)c->size ();
	}
	if (index >= long (c->size ()) || index < 0)
	{
		bs_throw_exception ((boost::format ("[%s] Index out of range: %d -> %d") % context % c->size () % index).str ());
	}

	return index;
}

template <typename T>
typename T::value_type get_item_array_ext (T *c, typename T::size_type i)
{
	return (*c)[get_index (c, i, "get_item_array_ext")];
}

template <typename T>
void set_item_array_ext (T *c, typename T::size_type i, const typename T::value_type &v)
{
	(*c)[get_index (c, i, "set_item_array_ext")] = v;
}

template <typename T>
void py_export_array_ext (const char *name)
{
	class_ <T> (name)
		.def ("__len__", make_function (&T::size))
		.def ("__iter__", bp::iterator <T> ())
		.def ("__getitem__", get_item_array_ext <T>)
		.def ("__setitem__", set_item_array_ext <T>)
		;
}

void py_export_idata()
{
	py_export_array_ext <array_float_t> ("array_ext_f");
	py_export_array_ext <array_uchar_t>   ("array_ext_i");

	py_export_idata_t< base_strategy_did >("idata_d");
	py_export_idata_t< base_strategy_fif >("idata_f");
	py_export_idata_t< base_strategy_dif >("idata_m");
}

template class py_idata< base_strategy_fif >;
template class py_idata< base_strategy_did >;
template class py_idata< base_strategy_dif >;

} }	// namespace blue_sky { python

