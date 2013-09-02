/**
 *       \file  py_calc_well.h
 *      \brief  Python wrapper for calc_well, see calc_well.h
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  23.06.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  A bit outdate
 * */
#ifndef PY_BS_CALC_WELL_H_
#define PY_BS_CALC_WELL_H_

#ifdef BSPY_EXPORTING_PLUGIN
#include "calc_well.h"
#include "well_connection.h"
#include "calc_model.h"
#include "py_calc_well_detail.h"
#include "facility_manager.h"

#include "export_python_wrapper.h"


namespace blue_sky {
namespace python {

  CLASS_WRAPPER (wells::connection, py_connection)
  {
  public:

    typedef t_double  item_t;
    typedef t_float   rhs_item_t;

    typedef wells::connection  wrapped_t;

  public:

    CLASS_WRAPPER_DECL (py_connection);

    WRAPPER_METHOD (clear_data, void, 0, (empty_arg__));
    WRAPPER_METHOD_R (get_rw_value, shared_vector <item_t>, 0, (empty_arg__));
    WRAPPER_METHOD_R (get_wr_value, shared_vector <item_t>, 0, (empty_arg__));
    WRAPPER_METHOD_R (get_rr_value, shared_vector <item_t>, 0, (empty_arg__));
    WRAPPER_METHOD_R (get_ps_value, shared_vector <item_t>, 0, (empty_arg__));
    WRAPPER_METHOD_R (get_rate_value, shared_vector <rhs_item_t>, 0, (empty_arg__));
  };

  CLASS_WRAPPER (well, py_well)
  {
  public:
    typedef smart_ptr <calc_model, true >         sp_calc_model_t;
    typedef smart_ptr <rs_mesh_iface, true >      sp_mesh_iface_t;

    typedef well                                  wrapped_t;

  public:

    CLASS_WRAPPER_DECL (py_well);

    WRAPPER_METHOD (process_impl,     void, 5, (bool, double, const sp_calc_model_t &, const sp_mesh_iface_t &, smart_ptr <jacobian> &));
    WRAPPER_METHOD (restore_solution, void, 4, (double, const spv_double &, const spv_double &, index_t));
    WRAPPER_METHOD (clear_data,       void, 0, (empty_arg__));
  };

  PY_EXPORTER (connection_exporter, default_exporter)
    .def ("clear_data",                 &T::clear_data)
    .add_property ("rw",                &T::get_rw_value)
    .add_property ("wr",                &T::get_wr_value)
    .add_property ("rr",                &T::get_rr_value)
    .add_property ("ps",                &T::get_ps_value)
    .add_property ("rate",              &T::get_rate_value)
    .add_property ("head_term",         &T::get_head_term, &T::set_head_term)
    .add_property ("cur_bhp",           &T::get_cur_bhp, &T::set_cur_bhp)
    .add_property ("connection_depth",  &T::get_connection_depth, &T::set_connection_depth)
    .add_property ("density",           &T::get_density, make_function (detail::set_connection_density<T>))
    .add_property ("bulkp",             &T::get_bulkp, &T::set_bulkp)
    .add_property ("n_block",           &T::n_block)
    .add_property ("i",                 &T::i_coord)
    .add_property ("j",                 &T::j_coord)
    .add_property ("k",                 &T::k_coord)
    .add_property ("fact",              &T::get_fact)
    .add_property ("gw",                &T::get_fact)
    .add_property ("is_shut",           &T::is_shut)
  PY_EXPORTER_END;

  PY_EXPORTER (well_exporter, default_exporter)
    .def ("process",                      &T::process)
    .def ("process_impl",                 &T::process_impl)
    .def ("update",                       &T::restore_solution)
    .def ("clear_data",                   &T::clear_data)
    .add_property ("name",                make_function (&T::get_name, return_value_policy <copy_const_reference> ()), &T::set_name)
    .add_property ("i_coord",             make_function (detail::get_well_i_coord<T>), make_function (detail::set_well_i_coord<T>))
    .add_property ("j_coord",             make_function (detail::get_well_j_coord<T>), make_function (detail::set_well_j_coord<T>))
    .add_property ("bhp_depth",           &T::get_bhp_depth, &T::set_bhp_depth)
    .add_property ("exploitation_factor", make_function (detail::get_well_wefac<T>), &T::set_exploitation_factor)
    .add_property ("state",               make_function (detail::get_well_state<T>))
    .add_property ("bhp_rate",            &T::bhp)
    .add_property ("is_work",             make_function (detail::get_well_is_work<T>), make_function (detail::set_well_is_work<T>))
    .add_property ("is_production",       make_function (detail::well_is_production<T>))
  PY_EXPORTER_END;

  void
  py_export_calc_well ();

} // namespace python
} // namespace blue_sky

#endif // BSPY_EXPORTING_PLUGIN
#endif  // #ifndef PY_BS_CALC_WELL_H_
