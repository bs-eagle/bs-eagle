/**
 *       \file  py_event_manager.cpp
 *      \brief  Python wrapper for event_manager, 
 *              see event_manager.h
 *     \author  Morozov Andrey
 *       \date  07.06.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"

#include "event_manager.h"
#include "py_event_manager.h"

#include "reservoir.h"
#include "facility_manager.h"
#include "well_connection.h"
#include "calc_model.h"

#ifdef BSPY_EXPORTING_PLUGIN
#include "export_python_wrapper.h"

#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <datetime.h>

using namespace boost;
using namespace boost::posix_time;

namespace blue_sky {
namespace python {

  PY_EXPORTER (event_manager_exporter, default_exporter)
    .add_property ("event_list", &T::event_list)
  PY_EXPORTER_END;

  template <typename T>
  struct event_list_iterator 
  {
    event_list_iterator (T *t)
    : list_ (t)
    , iterator_ (t->begin ())
    , iterator_end_ (t->end ())
    {
    }

    typename T::value_type
    next ()
    {
#ifdef _DEBUG
      if (iterator_end_ != list_->end ())
        {
          bs_throw_exception ("Event list iterator not more valid");
        }
#endif

      if (iterator_ == iterator_end_)
        {
          PyErr_SetString (PyExc_StopIteration, "No more data");
          boost::python::throw_error_already_set ();
        }

      return *(iterator_++);
    }

    T                     *list_;
    typename T::iterator  iterator_;
    typename T::iterator  iterator_end_;
  };

  template <typename T>
  event_list_iterator <T>
  get_event_list_iterator (T *t)
  {
    return event_list_iterator <T> (t);
  }
  template <typename T>
  size_t
  get_event_list_size (T *t)
  {
    return t->size ();
  }
  template <typename T>
  typename T::value_type
  get_event_list_item (T *list_, size_t idx)
  {
    if (idx >= list_->size ())
      {
        bs_throw_exception (boost::format ("Index out of bound (idx: %ld, size: %ld)") 
          % idx % list_->size ());
      }

    typename T::iterator it = list_->begin ();
    std::advance (it, idx);
    return *it;
  }

  template <typename T>
  void
  export_event_list_iterator (const char *name)
  {
    using namespace boost::python;

    class_ <event_list_iterator <T> > (name, init <T *> ())
      .def ("__iter__",     pass_through)
      .def ("next",         &event_list_iterator <T>::next)
      ;
  }

  template <typename T>
  void
  export_event_list (const char *name)
  {
    using namespace boost::python;

    class_ <T> (name)
      .def ("__iter__",     get_event_list_iterator <T>)
      .def ("__len__",      get_event_list_size <T>)
      .def ("__getitem__",  get_event_list_item <T>)
      ;

    export_event_list_iterator <T> (std::string (std::string (name) + "_iterator").c_str ());
  }

  template <typename T>
  void
  export_event_map (const char *name)
  {
    using namespace boost::python;

    class_ <T> (name)
      .def (map_indexing_suite <T> ())
      ;
  }

  struct boost_ptime_to_python_datetime
  {
    static PyObject *
    convert (const boost::posix_time::ptime &pt)
    {
      boost::gregorian::date date = pt.date ();
      boost::posix_time::time_duration td = pt.time_of_day ();

      int year    = date.year ();
      int month   = date.month ();
      int day     = date.day ();
      int hour    = td.hours ();
      int minute  = td.minutes ();
      int second  = td.seconds ();

      return PyDateTime_FromDateAndTime (year, month, day, hour, minute, second, 0);
    }
  };

  struct boost_ptime_from_python_datetime
  {
    boost_ptime_from_python_datetime ()
    {
      boost::python::converter::registry::push_back (
        &convertible, &construct,
        boost::python::type_id <boost::posix_time::ptime> ());
    }

    static void *
    convertible (PyObject *obj_ptr)
    {
      return PyDateTime_Check (obj_ptr) ? obj_ptr : 0;
    }

    static void
    construct (PyObject *obj_ptr, 
      boost::python::converter::rvalue_from_python_stage1_data *data)
    {
      const PyDateTime_DateTime *pydate = reinterpret_cast <PyDateTime_DateTime *> (obj_ptr);

      int year      = PyDateTime_GET_YEAR (pydate);
      int month     = PyDateTime_GET_MONTH (pydate);
      int day       = PyDateTime_GET_DAY (pydate);
      int hour      = PyDateTime_DATE_GET_HOUR (pydate);
      int minute    = PyDateTime_DATE_GET_MINUTE (pydate);
      int second    = PyDateTime_DATE_GET_SECOND (pydate);
      int mcsecond  = PyDateTime_DATE_GET_MICROSECOND (pydate);

      boost::gregorian::date d (year, month, day);
      boost::posix_time::time_duration td (hour, minute, second, 0);
      td += boost::posix_time::microseconds (mcsecond);

      using namespace boost::python;
      using namespace boost::posix_time;

      void *storage = ((converter::rvalue_from_python_storage <ptime>*)data)->storage.bytes;

      new (storage) boost::posix_time::ptime (d, td);

      data->convertible = storage;
    }
  };


  void
  py_export_event_manager ()
  {
    using namespace boost::python;

    PyDateTime_IMPORT;

    boost_ptime_from_python_datetime ();
    to_python_converter <boost::posix_time::ptime, boost_ptime_to_python_datetime> ();

    typedef std::list <smart_ptr <event_base> > event_list_t;

    export_event_list <event_list_t> ("event_list");
    export_event_map <std::map <boost::posix_time::ptime, event_list_t> > ("event_map");
    base_exporter <event_manager, event_manager_exporter>::export_class ("event_manager");
  }

} // namespace python
} // namespace blue_sky
#endif
