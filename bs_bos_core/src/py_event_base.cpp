/**
 *       \file  py_event_base.cpp
 *      \brief  Python wrappers for event_base, python events iterator,
 *              for event_base see event_base.h
 *     \author  Nikonov Max
 *       \date  17.10.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  A bit outdate
 * */

#include "stdafx.h"

#include "calc_model.h"
#include "reservoir.h"

#include "py_event_base.h"

//#include "py_reservoir.h"
//#include "py_calc_model.h"

//#include "facility_manager.h"
//#include "calc_model_types.h"
//#include "well_connection.h"

//#include BS_FORCE_PLUGIN_IMPORT ()
//#include "py_rs_mesh.h"
//#include BS_STOP_PLUGIN_IMPORT ()

using namespace boost::python;

namespace blue_sky
  {
  namespace python
    {

    template <typename strategy_t>
    py_event_base_iterator<strategy_t>::py_event_base_iterator (const elist_t &events1)
    {
      events.resize (0);
      events.assign (events1.begin (), events1.end ());
      iter = events.begin ();
    }

    template <typename strategy_t>
    py_event_base_iterator<strategy_t>::py_event_base_iterator (const this_t &iter1)
//        : events (iter.events)
//        , iter (iter.iter)
    {
      events.resize (0);
      events.assign (iter1.events.begin (), iter1.events.end ());
      for (const_iterator_t i = events.begin (); i != events.end (); ++i)
        if (*iter1.iter == *i)
          {
            iter = i;
            break;
          }
    }

    template <typename strategy_t>
    typename py_event_base_iterator<strategy_t>::py_event_base_t *
    py_event_base_iterator<strategy_t>::next ()
    {
      if (!events.size () || iter == events.end ())
        {
          PyErr_SetString(PyExc_StopIteration, "No more data.");
          boost::python::throw_error_already_set();
        }

      py_event_base_t *tmp = new py_event_base_t (**this);

      (*this)++;

      return tmp;
    }

    template <typename strategy_t>
    typename py_event_base_iterator<strategy_t>::reference
    py_event_base_iterator<strategy_t>::operator*() const
      {
        if (iter == events.end () || !events.size ())
          {
            BOSWARN (section::schedule, level::debug) << "end of the list" << bs_end;
            return py_event_base_t (NULL);
          }
        return py_event_base_t (*iter);
      }

    template <typename strategy_t>
    typename py_event_base_iterator<strategy_t>::this_t&
    py_event_base_iterator<strategy_t>::operator++()
    {
      ++iter;
      return *this;
    }

    template <typename strategy_t>
    typename py_event_base_iterator<strategy_t>::this_t
    py_event_base_iterator<strategy_t>::operator++(int)
    {
      this_t tmp (*this);
      ++iter;
      return tmp;
    }

    template <typename strategy_t>
    typename py_event_base_iterator<strategy_t>::this_t&
    py_event_base_iterator<strategy_t>::operator--()
    {
      --iter;
      return *this;
    }

    template <typename strategy_t>
    typename py_event_base_iterator<strategy_t>::this_t
    py_event_base_iterator<strategy_t>::operator--(int)
    {
      this_t tmp (*this);
      --iter;
      return tmp;
    }

    template <typename strategy_t>
    bool
    py_event_base_iterator<strategy_t>::operator ==(const this_t &ritr) const
      {
        return (iter == ritr.iter);
      }

    template <typename strategy_t>
    bool
    py_event_base_iterator<strategy_t>::operator !=(const this_t &ritr) const
      {
        return (iter != ritr.iter);
      }

    template <typename strategy_t>
    const typename py_event_base_iterator<strategy_t>::this_t&
    py_event_base_iterator<strategy_t>::operator =(const this_t &ritr)
    {
      events = ritr.events;
      iter = ritr.iter;
      return *this;
    }

    template <typename strategy_t>
    py_el_pair<strategy_t>::py_el_pair ()
    {}

    template <typename strategy_t>
    py_el_pair<strategy_t>::py_el_pair (const const_iterator_t &iter)
        : first (iter->first)
        , second (iter->second)
    {}

    template <typename strategy_t>
    py_el_pair<strategy_t>::py_el_pair (const this_t &iter)
        : first (iter.first)
        , second (iter.second)
    {}

    template <typename strategy_t>
    typename py_el_pair<strategy_t>::py_event_base_iterator_t
    py_el_pair<strategy_t>::list_begin ()
    {
      return py_event_base_iterator_t (second);
    }


    template <typename strategy_t>
    event_list_iterator <strategy_t>::event_list_iterator (const sp_event_manager_t &evm) //iterator_t &iter)
        : evm (evm)
        , iter (evm->event_list.begin ())
        //, pair_ (iter)
    {
    }

    template <typename strategy_t>
    event_list_iterator <strategy_t>::event_list_iterator (const this_t &iter1)
        : evm (iter1.evm)
        , iter (iter1.iter)
        //, pair_ (iter.iter)
    {
      /*for (const_iterator_t i = evm->event_list.begin (); i != evm->event_list.end (); ++i)
        if (*i == *(iter1.iter)) {
          iter = const_iterator_t (i);
          break;
        }*/
    }

    template <typename strategy_t>
    typename event_list_iterator <strategy_t>::reference
    event_list_iterator <strategy_t>::operator*() const
      {
        py_el_pair_t elp(iter);
        return elp;
      }

    template <typename strategy_t>
    typename event_list_iterator <strategy_t>::this_t &
    event_list_iterator <strategy_t>::operator++()
    {
      ++iter;
      return *this;
    }

    template <typename strategy_t>
    typename event_list_iterator <strategy_t>::this_t
    event_list_iterator <strategy_t>::operator++(int)
    {
      this_t tmp (*this);
      ++iter;
      return tmp;
    }

    template <typename strategy_t>
    typename event_list_iterator <strategy_t>::py_el_pair_t
    event_list_iterator <strategy_t>::next ()
    {
      if (iter == evm->event_list.end ())
        {
          BOSWARN (section::schedule, level::debug) << "Stop map" << bs_end;
          PyErr_SetString(PyExc_StopIteration, "No more data.");
          boost::python::throw_error_already_set();
        }

      py_el_pair_t tmp(iter);

      (*this)++;
      //++iter;

      return tmp;

    }

    template <typename strategy_t>
    typename event_list_iterator <strategy_t>::this_t &
    event_list_iterator <strategy_t>::operator--()
    {
      --iter;
      return *this;
    }

    template <typename strategy_t>
    typename event_list_iterator <strategy_t>::this_t
    event_list_iterator <strategy_t>::operator--(int)
    {
      this_t tmp (*this);
      --iter;
      return *this;
    }

    template <typename strategy_t>
    bool
    event_list_iterator <strategy_t>::operator ==(const this_t &ritr) const
      {
        return (iter == ritr.iter);
      }

    template <typename strategy_t>
    bool
    event_list_iterator <strategy_t>::operator !=(const this_t &ritr) const
      {
        return (iter != ritr.iter);
      }

    template <typename strategy_t>
    const typename event_list_iterator <strategy_t>::this_t &
    event_list_iterator <strategy_t>::operator =(const this_t &ritr)
    {
      iter = ritr.iter;
      evm = ritr.evm;
      return *this;
    }


    template <typename strategy_t>
    void export_py_event_base (const char *name)
    {
      class_< py_el_pair <strategy_t> > (std::string (std::string (name) + "_pair").c_str (), init <> ())
      .def ("list_begin", &py_el_pair<strategy_t>::list_begin)
      ;

      class_< event_list_iterator <strategy_t> > (std::string (std::string (name) + "_iterator").c_str (), no_init)
      .def ("__iter__",pass_through)
      .def ("next", &event_list_iterator <strategy_t>::next)
      ;

      class_< py_event_base_iterator <strategy_t> > (std::string (std::string (name) + "_iterator2").c_str (), no_init)
      .def ("__iter__",pass_through)
      .def ("next", &py_event_base_iterator <strategy_t>::next, return_value_policy <manage_new_object> ())
      ;
    }

    //PYTHON_EXPORT_WRAPPER (py_export_event_base, event_base, py_event_base, 1, (apply));

    PY_EXPORTER (event_base_exporter, default_exporter)
      .def ("apply", &T::apply)
    PY_EXPORTER_END;

    void py_export_events ()
    {
      strategy_exporter::export_class <py_event_base, event_base, event_base_exporter> ("event_base");

      export_py_event_base <base_strategy_fi> ("event_base_f");
      export_py_event_base <base_strategy_di> ("event_base_d");
      export_py_event_base <base_strategy_mixi> ("event_base_m");
    }

    template class py_el_pair <base_strategy_fi>;
    template class py_el_pair <base_strategy_di>;
    template class py_el_pair <base_strategy_mixi>;

    template class event_list_iterator <base_strategy_fi>;
    template class event_list_iterator <base_strategy_di>;
    template class event_list_iterator <base_strategy_mixi>;

    template class py_event_base_iterator <base_strategy_fi>;
    template class py_event_base_iterator <base_strategy_di>;
    template class py_event_base_iterator <base_strategy_mixi>;

} // namespace python
} // namespace blue_sky
