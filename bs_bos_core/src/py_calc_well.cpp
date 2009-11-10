/**
 * \file py_calc_well.cpp
 * \brief impl of
 * \author Sergey Miryanov
 * \date 23.06.2008
 * */
#include "stdafx.h"

#ifdef BSPY_EXPORTING_PLUGIN

#include "py_calc_well.h"
#include "reservoir.h"
#include "facility_manager.h"

namespace blue_sky {
namespace python {

  //////////////////////////////////////////////////////////////////////////
  //template <typename strategy_t>
  //py_well_iterator<strategy_t>::py_well_iterator(const this_t &src)
  //    : mgr(src.mgr), ins(src.ins) //,endi(src.endi)
  //{
  //}

  //template <typename strategy_t>
  //py_well_iterator<strategy_t>::py_well_iterator(const sp_facility_t &mgr) //const well_iterator_t &src, get_end fn)
  //    : mgr(mgr), ins(mgr->wells_begin ()) //ins(src),endi(fn)
  //{
  //}

  //template <typename strategy_t>
  //py_well_iterator<strategy_t>::~py_well_iterator() {}

  //template <typename strategy_t>
  //typename py_well_iterator<strategy_t>::reference
  //py_well_iterator<strategy_t>::operator*() const
  //  {
  //    if (ins != mgr->wells_end ())
  //      {
  //        py_well_t w(sp_well_t (ins->second, bs_dynamic_cast ()));
  //        return w;
  //      }
  //    return py_well_t(NULL);
  //  }

  //template <typename strategy_t>
  //typename py_well_iterator<strategy_t>::pointer
  //py_well_iterator<strategy_t>::operator->() const
  //  {
  //    if (ins != mgr->wells_end ())
  //      {
  //        py_well_t w(sp_well_t (ins->second, bs_dynamic_cast ()));
  //        return w;
  //      }
  //    return py_well_t(NULL);
  //  }

  //template <typename strategy_t>
  //typename py_well_iterator<strategy_t>::this_t&
  //py_well_iterator<strategy_t>::operator++()
  //{
  //  ++ins;
  //  return *this;
  //}

  //template <typename strategy_t>
  //typename py_well_iterator<strategy_t>::this_t
  //py_well_iterator<strategy_t>::operator++(int)
  //{
  //  py_well_iterator tmp(*this);
  //  ++ins;
  //  return tmp;
  //}

  //template <typename strategy_t>
  //typename py_well_iterator<strategy_t>::py_well_t
  //py_well_iterator<strategy_t>::next()
  //{
  //  if (ins == mgr->wells_end ())
  //    {
  //      PyErr_SetString(PyExc_StopIteration, "No more data.");
  //      boost::python::throw_error_already_set();
  //    }

  //  py_well_iterator tmp(*this);

  //  (*this)++;

  //  return *tmp;
  //}

  //template <typename strategy_t>
  //bool
  //py_well_iterator<strategy_t>::operator ==(const this_t &ritr) const
  //  {
  //    return (ins == ritr.ins);
  //  }

  //template <typename strategy_t>
  //bool
  //py_well_iterator<strategy_t>::operator !=(const this_t &ritr) const
  //  {
  //    return (!(ins == ritr.ins));
  //  }

  //template <typename strategy_t>
  //const typename py_well_iterator<strategy_t>::this_t&
  //py_well_iterator<strategy_t>::operator =(const this_t &ritr)
  //{
  //  ins = ritr.ins;
  //  mgr = ritr.mgr;
  //  return *this;
  //}

  /*********************************************************************/
  template <typename strategy_t>
  void py_export_well_iterator (const char *name)
  {
    using namespace boost::python;

    class_ <typename facility_manager <strategy_t>::well_const_iterator_t> (name, no_init)
    .def ("__iter__",pass_through)
    //.def ("next", &py_well_iterator <strategy_t>::next)
    ;
  }

  void
  py_export_calc_well ()
  {
    strategy_exporter::export_base <well, well_exporter> ("well_seq");
    strategy_exporter::export_base <wells::connection, connection_exporter> ("connection_seq");

    strategy_exporter::export_class <py_well, well, well_exporter> ("py_well_seq");
    strategy_exporter::export_class <py_connection, wells::connection, connection_exporter> ("py_connection_seq");

    //using namespace boost::python;
    //class_ < py_well_storage<base_strategy_fi>::py_well_list_t> ("vector_well_seq_fi")
    //  .def (vector_indexing_suite < py_well_storage<base_strategy_fi>::py_well_list_t> ())
    //  ;
    //using namespace boost::python;
    //class_ < py_well_storage<base_strategy_di>::py_well_list_t> ("vector_well_seq_di")
    //  .def (vector_indexing_suite < py_well_storage<base_strategy_di>::py_well_list_t> ())
    //  ;

    py_export_well_iterator <base_strategy_fi> ("well_iterator_f");
    py_export_well_iterator <base_strategy_di> ("well_iterator_d");
    py_export_well_iterator <base_strategy_mixi> ("well_iterator_m");
  }

  //template class py_well_iterator <base_strategy_fi>;
  //template class py_well_iterator <base_strategy_di>;
  //template class py_well_iterator <base_strategy_mixi>;

#endif  // #ifdef BSPY_EXPORTING_PLUGIN
} // namespace python
} // namespace blue_sky

