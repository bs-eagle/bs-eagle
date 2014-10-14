/**
 *       \file  py_facility_manager.cpp
 *      \brief  Python wrappers for facility_manager,
 *              see facility_manager.h
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  17.10.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */

#ifdef BSPY_EXPORTING_PLUGIN
#include "py_facility_manager.h"
#include "reservoir.h"
#include "facility_manager.h"
#include "export_python_wrapper.h"
#include "py_bs_iterator.h"

using namespace boost::python;

namespace blue_sky {
namespace python {

  template <typename facility_manager_t>
  struct well_iterator
  {
    typedef typename facility_manager_t::sp_well_t              sp_well_t;
    typedef typename facility_manager_t::well_const_iterator_t  iterator_t;

    well_iterator (facility_manager_t *facility_manager_)
    : facility_manager_ (facility_manager_)
    , iterator_ (facility_manager_->wells_begin ())
    , iterator_end_ (facility_manager_->wells_end ())
    {
    }

    sp_well_t
    next () 
    {
#ifdef _DEBUG
      if (iterator_end_ != facility_manager_->wells_end ())
        {
          bs_throw_exception ("Well iterator not more valid");
        }
#endif
      while (iterator_ != iterator_end_)
        {
          sp_well_t well (iterator_->second, bs_dynamic_cast ());
          ++iterator_;

          if (well)
            return well;
        }
      if (iterator_ == iterator_end_)
        {
          PyErr_SetString (PyExc_StopIteration, "No more data");
          boost::python::throw_error_already_set ();
        }
      bs_throw_exception ("Not reachable state reachabled");
    }

    smart_ptr <facility_manager_t>  facility_manager_;
    iterator_t                      iterator_;
    iterator_t                      iterator_end_;
  };

  template <typename T>
  well_iterator <T>
  get_wells (T *t)
  {
    return well_iterator <T> (t);
  }

  PY_EXPORTER (facility_manager_exporter, default_exporter)
    .def ("add_well",       &T::add_well)
    .add_property ("wells", get_wells <T>)
  PY_EXPORTER_END;

  void
  export_well_iterator (const char *name)
  {
    typedef well_iterator <facility_manager> T;

    using namespace boost::python;
    class_ <T> (name, init <facility_manager *> ())
      .def ("next", &T::next)
      .def ("__iter__", pass_through)
      ;
  }

  void py_export_facility_manager()
  {
    base_exporter <facility_manager, facility_manager_exporter>::export_class ("facility_manager");
    export_well_iterator ("well_iterator");
  }

} // namespace python
} // namespace blue_sky
#endif
