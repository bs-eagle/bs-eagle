/**
 *       \file  py_well_factory.h
 *      \brief  Python wrappers for well factories (well, well_controller,
 *              well_rate_controller, well_limit_operation factories),
 *              see calc_well.h and related files
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  21.05.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_BOS_CORE_PY_WELL_FACTORY_H_
#define BS_BOS_CORE_PY_WELL_FACTORY_H_

#ifdef BSPY_EXPORTING_PLUGIN
#include "calc_well.h"
#include "well_controller.h"
#include "well_limit_operation.h"
#include "calc_model.h"
#include "export_python_wrapper.h"


namespace blue_sky {
namespace python {

  STRATEGY_CLASS_WRAPPER (well_factory, py_well_factory)
  {
  public:
    typedef smart_ptr <well <strategy_t>, true>               sp_well_t;
    typedef smart_ptr <wells::connection <strategy_t>, true>  sp_connection_t;
    typedef well_factory <strategy_t>                         wrapped_t;
  public:
    STRATEGY_CLASS_WRAPPER_DECL (py_well_factory);
    WRAPPER_METHOD_R_CONST (create_well, sp_well_t, 2, (std::string, std::string));
    WRAPPER_METHOD_R_CONST (create_connection, sp_connection_t, 0, (empty_arg__));
  };

  using namespace wells;
  STRATEGY_CLASS_WRAPPER (well_controller_factory, py_well_controller_factory)
  {
  public:
    typedef smart_ptr <well_controller <strategy_t>, true>    sp_well_controller_t;
    typedef smart_ptr <well_rate_control <strategy_t>, true>  sp_well_control_t;
    typedef smart_ptr <calc_model <strategy_t>, true>         sp_calc_model_t;
    typedef well_controller_factory <strategy_t>              wrapped_t;

  public:
    STRATEGY_CLASS_WRAPPER_DECL (py_well_controller_factory);
    WRAPPER_METHOD_R_CONST (create_controller, sp_well_controller_t, 0, (empty_arg__));
  };

  CLASS_WRAPPER (well_limit_operation_factory, py_well_limit_operation_factory)
  {
  public:
    typedef smart_ptr <well_limit_operation, true>  sp_well_limit_operation_t;
    typedef well_limit_operation_factory            wrapped_t;
  public:
    CLASS_WRAPPER_DECL (py_well_limit_operation_factory);
    WRAPPER_METHOD_R (create_limit, sp_well_limit_operation_t, 1, (limit_operation_type));
  };

  /**
   * \brief  Exports wrappers to python
   * */
  void 
  py_export_well_factories ();

} // namespace python
} // namespace blue_sky

#endif
#endif  // #ifndef BS_BOS_CORE_PY_WELL_FACTORY_H_