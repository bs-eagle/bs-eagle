/**
 *       \file  py_two_stage_preconditioner.cpp
 *      \brief  Python wrapper for two_stage_preconditioner
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  16.04.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#include "stdafx.h"

#include "py_two_stage_preconditioner.h"

#ifdef BSPY_EXPORTING_PLUGIN

namespace blue_sky
  {
  namespace python
    {

    template <class solver_t>
    struct two_stage_prec_exporter
    {
      template <typename class_t>
      static class_t &
      export_class (class_t &class__)
      {
        using namespace boost::python;

        solver_exporter <solver_t>::export_class (class__)
          .def ("set_prec_1",   &solver_t::set_prec_1)
          .def ("set_prec_2",   &solver_t::set_prec_2)
          .add_property ("prec_1", &solver_t::get_prec_1, &solver_t::set_prec_1)
          .add_property ("prec_2", &solver_t::get_prec_2, &solver_t::set_prec_2)
          ;

        return class__;
      }
    };

    //! export classes to python
    void py_export_two_stage_prec ()
    {
      strategy_exporter::export_class <two_stage_preconditioner, linear_solver_base, two_stage_prec_exporter> ("two_stage_prec_seq");
    }


  } // namespace python
} // namespace blue_sky

#endif // #ifdef BSPY_EXPORTING_PLUGIN
