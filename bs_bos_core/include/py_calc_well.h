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

  STRATEGY_CLASS_WRAPPER (wells::connection, py_connection)
  {
  public:

    typedef typename strategy_t::item_t     item_t;
    typedef typename strategy_t::rhs_item_t rhs_item_t;

    typedef wells::connection <strategy_t>  wrapped_t;

  public:

    STRATEGY_CLASS_WRAPPER_DECL (py_connection);

    WRAPPER_METHOD (clear_data, void, 0, (empty_arg__));
    WRAPPER_METHOD_R (get_rw_value, array_ext <item_t>, 0, (empty_arg__));
    WRAPPER_METHOD_R (get_wr_value, array_ext <item_t>, 0, (empty_arg__));
    WRAPPER_METHOD_R (get_rr_value, array_ext <item_t>, 0, (empty_arg__));
    WRAPPER_METHOD_R (get_ps_value, array_ext <item_t>, 0, (empty_arg__));
    WRAPPER_METHOD_R (get_rate_value, array_ext <rhs_item_t>, 0, (empty_arg__));
  };

  STRATEGY_CLASS_WRAPPER (well, py_well)
  {
  public:

    typedef typename strategy_t::index_t                      index_t;
    typedef typename strategy_t::item_t                       item_t;
    typedef typename strategy_t::rhs_item_t                   rhs_item_t;
    typedef typename strategy_t::item_array_t                 item_array_t;
    typedef typename strategy_t::rhs_item_array_t             rhs_item_array_t;

    typedef smart_ptr <calc_model <strategy_t>, true >        sp_calc_model_t;
    typedef smart_ptr <rs_mesh_iface <strategy_t>, true >           sp_mesh_iface_t;
    typedef smart_ptr <jacobian_matrix <strategy_t>, true >   sp_jmatrix_t;

    typedef well <strategy_t>                                 wrapped_t;

  public:

    STRATEGY_CLASS_WRAPPER_DECL (py_well);

    WRAPPER_METHOD_R (get_ww_value,   array_ext <item_t>, 0, (empty_arg__));
    WRAPPER_METHOD_R (get_bw_value,   array_ext <item_t>, 0, (empty_arg__));
    WRAPPER_METHOD (eliminate,        void, 5, (rhs_item_t *, index_t, index_t, item_t, index_t));
    WRAPPER_METHOD (process,          void, 5, (int, double, const sp_calc_model_t &, const sp_mesh_iface_t &, sp_jmatrix_t &));
    WRAPPER_METHOD (process_newton,   void, 4, (int, const sp_calc_model_t &, const sp_mesh_iface_t &, sp_jmatrix_t &));
    WRAPPER_METHOD (restore_solution, void, 4, (double, const item_array_t &, const item_array_t &, index_t));
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
    //.def ("eliminate",                    &T::eliminate)
    .def ("process",                      &T::process)
    .def ("update",                       &T::restore_solution)
    .def ("clear_data",                   &T::clear_data)
    //.def ("get_connection",               &T::get_connection)
    .add_property ("ww",                  &T::get_ww_value)
    .add_property ("bw",                  &T::get_bw_value)
    .add_property ("name",                make_function (&T::get_name, return_value_policy <copy_const_reference> ()), &T::set_name)
    .add_property ("i_coord",             make_function (detail::get_well_i_coord<T>), make_function (detail::set_well_i_coord<T>))
    .add_property ("j_coord",             make_function (detail::get_well_j_coord<T>), make_function (detail::set_well_j_coord<T>))
    .add_property ("bhp_depth",           &T::get_bhp_depth, &T::set_bhp_depth)
    .add_property ("exploitation_factor", make_function (detail::get_well_wefac<T>), &T::set_exploitation_factor)
    .add_property ("state",               make_function (detail::get_well_state<T>))
    .add_property ("bhp_rate",            &T::bhp)
    .add_property ("is_work",             make_function (detail::get_well_is_work<T>), make_function (detail::set_well_is_work<T>))
    .add_property ("connections",         &T::get_connections_count)
    .add_property ("is_production",       make_function (detail::well_is_production<T>))
  PY_EXPORTER_END;

  //template <typename strategy_t>
  //class py_well_iterator
  //      : public std::iterator<
  //      std::forward_iterator_tag,
  //      py_well<strategy_t>, ptrdiff_t,
  //      py_well<strategy_t>, py_well<strategy_t> >
  //  {
  //    friend class py_well <strategy_t>;

  //  public:
  //    typedef py_well <strategy_t>              py_well_t;

  //    typedef py_well_iterator <strategy_t>     this_t;

  //    typedef std::iterator<
  //    std::forward_iterator_tag,
  //    //std::bidirectional_iterator_tag,
  //    py_well<strategy_t>, ptrdiff_t,
  //    py_well<strategy_t>, py_well<strategy_t> >  base_t;

  //    typedef typename base_t::reference  reference;
  //    typedef typename base_t::pointer    pointer;

  //    typedef facility_manager <strategy_t>                       facility_manager_t;
  //    typedef smart_ptr <facility_manager_t>                      sp_facility_t;

  //    typedef this_t                                              wrapped_t;
  //    typedef typename facility_manager_t::well_const_iterator_t  well_const_iterator_t;
  //    typedef well <strategy_t>                                   well_t;
  //    typedef smart_ptr <well_t, true>                            sp_well_t;

  //    sp_facility_t   mgr;
  //    well_const_iterator_t ins;

  //  public:
  //    py_well_iterator(const py_well_iterator &src);
  //    py_well_iterator(const sp_facility_t &mgr);
  //    ~py_well_iterator();

  //    reference operator*() const;
  //    pointer operator->() const;

  //    py_well_iterator& operator++();
  //    py_well_iterator operator++(int);

  //    py_well_t next ();

  //    //py_well_iterator& operator--();
  //    //py_well_iterator operator--(int);

  //    bool operator ==(const py_well_iterator &ritr) const;
  //    bool operator !=(const py_well_iterator &ritr) const;
  //    const this_t &operator =(const this_t &ritr);
  //  };

  void
  py_export_calc_well ();

} // namespace python
} // namespace blue_sky

#endif // BSPY_EXPORTING_PLUGIN
#endif  // #ifndef PY_BS_CALC_WELL_H_
