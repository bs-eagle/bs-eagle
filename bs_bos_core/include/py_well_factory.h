/**
 * \file py_well_factory.h
 * \brief 
 * \author Sergey Miryanov
 * \date 21.05.2009
 * */
#ifndef BS_BOS_CORE_PY_WELL_FACTORY_H_
#define BS_BOS_CORE_PY_WELL_FACTORY_H_

#include "calc_well.h"
#include "well_controller.h"
#include "well_limit_operation.h"
#include "calc_model.h"
#include "export_python_wrapper.h"
#include "well_rate_control_interface.h"


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
    WRAPPER_METHOD_R_CONST (create_control, sp_well_control_t, 3, (rate_control_type, bool, const sp_calc_model_t &));
  };

  STRATEGY_CLASS_WRAPPER (well_rate_control_factory, py_well_rate_control_factory)
  {
  public:
    typedef well_rate_control_interface <strategy_t>          well_rate_control_t;
    typedef smart_ptr <calc_model <strategy_t>, true>         sp_calc_model_t;
    typedef smart_ptr <well_rate_control_t, true>             sp_well_rate_control_t;
    typedef well_rate_control_factory <strategy_t>            wrapped_t;

  public:
    STRATEGY_CLASS_WRAPPER_DECL (py_well_rate_control_factory);
    WRAPPER_METHOD_R (create_control, sp_well_rate_control_t, 4, (rate_control_type, bool, bool, const sp_calc_model_t &));
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

  void py_export_well_factories ();

} // namespace python
} // namespace blue_sky


#endif  // #ifndef BS_BOS_CORE_PY_WELL_FACTORY_H_

