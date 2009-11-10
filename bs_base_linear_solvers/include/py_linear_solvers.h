/**
 * \file   py_linear_solvers.h
 * \brief  Python wrapper for linear solvers
 * \author Miryanov Sergey
 * \date 2008-04-04
 */

#ifndef BS_LINEAR_SOLVERS_PYTHON_H
#define BS_LINEAR_SOLVERS_PYTHON_H

#include "strategies.h"
#include "linear_solvers.h"
#include "cgs.h"
#include "tfqmr.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "py_bcsr_matrix.h"
#include "write_time_to_log.h"
#include "dummy_base.h"
#include "construct_python_object.h"
#include "make_me_happy.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include "export_python_wrapper.h"

#ifdef BSPY_EXPORTING_PLUGIN
namespace blue_sky {
namespace python {

  /**
   * \brief python wrapper for most linear solvers
   */
  STRATEGY_CLASS_WRAPPER (linear_solver_base, py_linear_solver)
  {
  public:
    typedef linear_solver_base <strategy_t>                   linear_solver_base_t;

    typedef typename strategy_t::matrix_t                     matrix_t;     ///< short name for matrix type from wrapped type
    typedef typename strategy_t::item_array_t                 item_array_t;      ///< short name for array type from wrapped type
    typedef typename strategy_t::index_array_t                index_array_t;
    typedef typename strategy_t::rhs_item_t	                  rhs_item_t;     ///< short name for type of rhs
    typedef typename strategy_t::rhs_item_array_t             rhs_item_array_t;   ///< short name for rhs array type

    typedef py_matrix_base <rhs_item_array_t, index_array_t>  py_matrix_t;

    typedef smart_ptr<linear_solver_base_t, true>             sp_linear_solver_base_t;
    typedef smart_ptr<linear_solver_base_t, true>             sp_prec_t;      ///< short name to smart pointer to this class
    typedef smart_ptr<linear_solver_prop, true>               sp_prop_t;      ///< short name to smart pointer to properties holder class

    typedef linear_solver_base_t                              wrapped_t;
    typedef py_linear_solver_base <strategy_t>                base_t;

  public:

    STRATEGY_CLASS_WRAPPER_DECL (py_linear_solver);

    WRAPPER_METHOD_R (setup, int, 1, (matrix_t *));
    WRAPPER_METHOD_R (solve, int, 3, (matrix_t *, rhs_item_array_t &, item_array_t &));
    WRAPPER_METHOD_R (solve_prec, int, 3, (matrix_t *, item_array_t &, item_array_t &));

    void
    init_object (linear_solver_base_t *solver)
    {
      solver_ = sp_linear_solver_base_t (solver);
    }

  private:

    sp_linear_solver_base_t solver_;
  };

#define WRAP_STD_SOLVER(name_, py_name_)                                                  \
  STRATEGY_CLASS_WRAPPER (name_, py_name_)                                                \
  {                                                                                       \
  public:                                                                                 \
                                                                                          \
    typedef typename strategy_t::matrix_t                 matrix_t;                       \
    typedef typename strategy_t::item_array_t             item_array_t;                   \
    typedef typename strategy_t::rhs_item_array_t         rhs_item_array_t;               \
                                                                                          \
    typedef name_ <strategy_t>                            wrapped_t;                      \
    typedef BOOST_PP_CAT (py_name_, _base) <strategy_t>   base_t;                         \
                                                                                          \
  public:                                                                                 \
                                                                                          \
    STRATEGY_CLASS_WRAPPER_DECL (py_name_);                                               \
                                                                                          \
    WRAPPER_METHOD_R (setup, int, 1, (matrix_t *));                                       \
    WRAPPER_METHOD_R (solve, int, 3, (matrix_t *, rhs_item_array_t &, item_array_t &));   \
    WRAPPER_METHOD_R (solve_prec, int, 3, (matrix_t *, item_array_t &, item_array_t &));  \
  };

  WRAP_STD_SOLVER (gmres_solver2, py_gmres_solver);
  WRAP_STD_SOLVER (bicgstab_solver, py_bicgstab_solver);
  WRAP_STD_SOLVER (cgs_solver, py_cgs_solver);
  WRAP_STD_SOLVER (tfqmr_solver, py_tfqmr_solver);

  PY_EXPORTER (solver_exporter, default_exporter)
    .def ("solve",      &T::solve)
    .def ("solve_prec", &T::solve_prec)
    .def ("setup",      &T::setup)
    .def ("set_prec",   &T::set_prec)
    .def ("set_prop",   &T::set_prop)
  PY_EXPORTER_END;

  //! export linear_solver_prop to python
  void py_export_linear_solver_prop ();

  //! export linear solvers to python
  void py_export_linear_solvers ();

} // namespace python
}	// namespace blue_sky

#endif // #ifdef BSPY_EXPORTING_PLUGIN
#endif //#ifndef BS_LINEAR_SOLVERS_PYTHON_H
