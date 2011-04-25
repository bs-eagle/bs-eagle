#ifndef PY_LSOLVER_IFACE_H_
#define PY_LSOLVER_IFACE_H_
/**
 * @file py_lsolver_iface.h
 * @brief
 * @author
 * @date 2009-11-25
 */
#include <string>
#include "lsolver_iface.h"
#include "amg_smoother_iface.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "dummy_base.h"
#include "construct_python_object.h"
#include "make_me_happy.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include "export_python_wrapper.h"

#ifdef BSPY_EXPORTING_PLUGIN
namespace blue_sky
  {
  namespace python
    {

    // Sergey Miryanov at 07.04.2008
    // Refactored at 16.10.2009
    // for casting into python from child classes to parent class
    // we should export base class
#if 0
  /**
   * @brief python wrapper for lsolver_iface class
   *
   * @param name_           -- should be matrix_iface
   * @param py_name_        -- name of the class in python
   *
   * @return nothing
   */
  CLASS_WRAPPER_T (1, (strategy_t), lsolver_iface, py_lsolver_iface)
  {
  public:

    CLASS_WRAPPER_DECL_T (1, (strategy_t), py_lsolver_iface)
  public:

    typedef typename strategy_t::matrix_t                       matrix_t;
    typedef typename strategy_t::t_double                      t_double;
    typedef typename strategy_t::i_type_t                       i_type_t;
    typedef typename strategy_t::fp_vector_type                 fp_vector_type_t;
    typedef lsolver_iface <strategy_t>                          this_t;
    typedef prop_iface<t_double, i_type_t, std::string, bool>  prop_t;
    typedef smart_ptr<this_t, true>                             sp_this_t;
    typedef smart_ptr<prop_t, true>                             sp_prop_t;

    WRAP_PURE_METHOD_R  (solve, int, 3, (const matrix_t&, const fp_vector_type_t&, fp_vector_type_t&));
    WRAP_PURE_METHOD_R  (solve_prec, int, 3, (const matrix_t&, const fp_vector_type_t&, fp_vector_type_t&));
    WRAP_PURE_METHOD_R  (setup, int, 1, (matrix_t&));
    WRAP_PURE_METHOD    (set_prec, void, 1, (sp_this_t&));
    WRAP_PURE_METHOD    (set_prop, void, 1, (sp_prop_t&));
    WRAP_PURE_METHOD_R  (get_prop, sp_prop_t&, 0, (empty_arg__));
    WRAP_PURE_METHOD_R_CONST  (get_final_residual, t_double, 0, (empty_arg__));
    WRAP_PURE_METHOD_R_CONST  (get_niters, int, 0, (empty_arg__));
    WRAP_PURE_METHOD_R_CONST  (py_str, std::string, 0, (empty_arg__));

  };

  /**
   * @brief python wrapper for amg_smoother_iface class
   *
   * @param name_           -- should be matrix_iface
   * @param py_name_        -- name of the class in python
   *
   * @return nothing
   */
  CLASS_WRAPPER_T (1, (strategy_t), amg_smoother_iface, py_amg_smoother_iface)
  {
  public:

    CLASS_WRAPPER_DECL_T (1, (strategy_t), py_amg_smoother_iface)
  public:

    typedef typename strategy_t::matrix_t                       matrix_t;
    typedef typename strategy_t::t_double                      t_double;
    typedef typename strategy_t::i_type_t                       i_type_t;
    typedef typename strategy_t::fp_vector_type                 fp_vector_type_t;
    typedef typename strategy_t::fp_storage_vector_type         fp_storage_vector_type_t;
    typedef typename strategy_t::i_vector_type                  i_vector_type_t;
    typedef bcsr_matrix_iface<fp_vector_type_t, i_vector_type_t, fp_storage_vector_type_t> bcsr_t;
    typedef lsolver_iface <strategy_t>                          this_t;
    typedef prop_iface<t_double, i_type_t, std::string, bool>  prop_t;
    typedef smart_ptr<this_t, true>                             sp_this_t;
    typedef smart_ptr<prop_t, true>                             sp_prop_t;

    WRAP_PURE_METHOD_R  (solve, int, 3, (const matrix_t&, const fp_vector_type_t&, fp_vector_type_t&));
    WRAP_PURE_METHOD_R  (solve_prec, int, 3, (const matrix_t&, const fp_vector_type_t&, fp_vector_type_t&));
    WRAP_PURE_METHOD_R  (setup, int, 1, (matrix_t&));
    WRAP_PURE_METHOD_R  (smooth, int, 5, (const matrix_t&, const i_vector_type_t &, const i_type_t, const fp_vector_type_t&, fp_vector_type_t&));
    WRAP_PURE_METHOD    (set_prec, void, 1, (sp_this_t&));
    WRAP_PURE_METHOD    (set_prop, void, 1, (sp_prop_t&));
    WRAP_PURE_METHOD_R  (get_prop, sp_prop_t&, 0, (empty_arg__));
    WRAP_PURE_METHOD_R_CONST  (get_final_residual, t_double, 0, (empty_arg__));
    WRAP_PURE_METHOD_R_CONST  (get_niters, int, 0, (empty_arg__));
    WRAP_PURE_METHOD_R_CONST  (py_str, std::string, 0, (empty_arg__));

  };
#endif //0

  PY_EXPORTER (py_lsolver_iface_exporter, default_exporter)
    .def ("solve",                &T::solve, args ("matrix", "rhs", "solution"), "Solve linear system")
    .def ("solve_prec",           &T::solve_prec, args ("matrix", "rhs", "solution"), "Make only one iteration of method (for exact methods equil to solve)")
    .def ("setup",                &T::setup, args ("matrix"), "Prepare method for solution (should be done before calling the solve method")
    .def ("set_prec",             &T::set_prec, args ("sp_prec"), "Set preconditioner")
    .def ("set_prop",             &T::set_prop, args ("sp_prop"), "Set properties")
    .def ("get_prop",             &T::get_prop, args (""), "Return smart pointer to the properties")
    .def ("get_final_residual",   &T::get_final_residual, args (""), "Return final residual (for exact methods 0)")
    .def ("get_niters",           &T::get_niters, args (""), "Return number of iterations (for exact methods 1)")
    .def ("__str__",              &T::py_str)
  PY_EXPORTER_END;

  PY_EXPORTER (py_amg_smoother_iface_exporter, py_lsolver_iface_exporter)
    .def ("smooth",              &T::smooth,
        args ("matrix", "cf_markers", "iter_num", "rhs", "solution"),
        "Smooth solution (make only one iteration)")
  PY_EXPORTER_END;

  //! export matrices to python
  void py_export_lsolvers ();


  } // namespace python
} // namespace blue_sky
#endif // #ifdef BSPY_EXPORTING_PLUGIN
#endif // #ifndef PY_LSOLVER_IFACE_H_
