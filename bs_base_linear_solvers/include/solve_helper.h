/**
 * \file solve_helper.h
 * \brief impl of
 * \author Elmira Salimgareeva
 * \date 31.03.2009
 * */

#ifndef BS_SOLVE_HELPER_H_
#define BS_SOLVE_HELPER_H_
namespace blue_sky
{

  template <class solver_t>
  int
  solve_helper (solver_t *solver, typename solver_t::matrix_t *mx, seq_vector <float> &rhs, typename solver_t::item_array_t &sol)
    {
      solver->solve (mx, rhs, sol);
      return 0;
    }

  template <class solver_t>
  int
  solve_helper (solver_t *solver, typename solver_t::matrix_t *mx, seq_vector <double> &rhs, typename solver_t::item_array_t &sol)
    {
      solver->solve_prec (mx, rhs, sol);
      return 0;
    }

}
#endif //#ifndef BS_SOLVE_HELPER_H_
