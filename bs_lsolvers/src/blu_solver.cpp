/** 
 * @file blu_solver.cpp
 * @brief block solver for linear systems using LU decomposition
 * @date 2009-12-09
 */

#include "bs_kernel.h"
#include "bs_assert.h"

#include "blu_solver.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "strategies.h"
#include "matrix_iface.h"
#include "matrix_macroses.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {
  /*!
   * \brief constructor
   */
  template <class strategy_t>
  blu_solver<strategy_t>::blu_solver (bs_type_ctor_param /*param*/)
  {
      prop = BS_KERNEL.create_object ("prop");
      if (!prop)
        {
          bs_throw_exception ("Type (prop) not registered");
        }
      init_prop ();
  }

  template <class strat_t>
  blu_solver<strat_t>::blu_solver(const blu_solver& solver)
      : bs_refcounter (solver) //, lsolver_iface <strat_t> ()
  {
    if (&solver != this)
      *this = solver;
  }

  //! set solver's properties
  template <class strat_t>
  void blu_solver<strat_t>::set_prop (sp_prop_t prop_)
  {
    prop = prop_;

    init_prop ();
  }

  template <class strat_t> void
  blu_solver<strat_t>::init_prop ()
    {
      block_size_idx = prop->get_index_i (std::string ("block_size"));
      if (block_size_idx < 0)
        block_size_idx = prop->add_property_i (true, std::string ("block_size"), 
                                               std::string ("Size of calculation block"));

      if (block_size_idx < 0)
        {
          bs_throw_exception ("Can not regidter some properties");
        }
    }

  template <class strategy_t>
  int blu_solver<strategy_t>::solve (sp_matrix_t matrix, sp_fp_array_t sp_rhs, sp_fp_array_t sp_sol)
  {
    BS_ASSERT (matrix);
    BS_ASSERT (sp_rhs->size ());
    BS_ASSERT (sp_sol->size ());
    BS_ASSERT (sp_rhs->size () == sp_sol->size ()) (sp_rhs->size ()) (sp_sol->size ());

    i_type_t n;     
    i_type_t nb;    

    sp_dens_matrix_t ilu;
    if (dynamic_cast<dens_matrix_iface_t *> (matrix.lock ()))
      {
        ilu = matrix;
        BS_ASSERT (ilu);
      }


    sp_fp_storage_array_t sp_values = ilu->get_values   ();
    fp_storage_type_t *values = &(*sp_values)[0];
    fp_type_t *sol = &(*sp_sol)[0];
    fp_type_t *rhs = &(*sp_rhs)[0];

    
    n         = ilu->get_n_rows ();
    nb        = ilu->get_calc_block_size ();
    if (nb < 2)
      nb = n;

    BS_ASSERT (n == (i_type_t)sp_sol->size ());

    memcpy (sol, rhs, n * sizeof (fp_type_t));
     
    return impl.lu_find_L_roots (n, values, sol, nb) 
           || impl.lu_find_U_roots (n, values, sol, nb);
  }

  template <class strategy_t>
  int blu_solver<strategy_t>::solve_prec(sp_matrix_t matrix, sp_fp_array_t rhs, sp_fp_array_t sol)
  {
     return solve (matrix, rhs, sol);
  }

  /**
   * \brief setup preconditioner (merge matrices if needed)
   *
   * \param matrix Various matrix
   * \return 0 if success
   */
  template <class strategy_t> int
  blu_solver<strategy_t>::setup (sp_matrix_t matrix)
  {
    BS_ASSERT (matrix);

    i_type_t n;     
    i_type_t nb;    

    sp_dens_matrix_t ilu;
    if (dynamic_cast<dens_matrix_iface_t *> (matrix.lock ()))
      {
        ilu = matrix;
        BS_ASSERT (ilu);
      }

    sp_fp_storage_array_t values = ilu->get_values   ();
    
    n         = ilu->get_n_rows ();
    nb        = ilu->get_calc_block_size ();
    if (nb < 2)
      nb = n;

    return impl.lu_decomposition (n, &(*values)[0], nb);
  }


  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE_T_DEF(blu_solver, (class));
  BLUE_SKY_TYPE_STD_COPY_T_DEF(blu_solver, (class));

  BLUE_SKY_TYPE_IMPL_T_EXT(1, (blu_solver<base_strategy_fif>) , 1, (lsolver_iface<base_strategy_fif>),  "blu_solver_fif",     "BLU solver", "Block LU solver for dense matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (blu_solver<base_strategy_did>) , 1, (lsolver_iface<base_strategy_did>),  "blu_solver_did",     "BLU solver", "Block LU solver for dense matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (blu_solver<base_strategy_dif>) , 1, (lsolver_iface<base_strategy_dif>),  "blu_solver_dif",     "BLU solver", "Block LU solver for dense matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (blu_solver<base_strategy_flf>) , 1, (lsolver_iface<base_strategy_flf>),  "blu_solver_flf",     "BLU solver", "Block LU solver for dense matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (blu_solver<base_strategy_dld>) , 1, (lsolver_iface<base_strategy_dld>),  "blu_solver_dld",     "BLU solver", "Block LU solver for dense matricies", false);
  BLUE_SKY_TYPE_IMPL_T_EXT(1, (blu_solver<base_strategy_dlf>) , 1, (lsolver_iface<base_strategy_dlf>),  "blu_solver_dlf",     "BLU solver", "Block LU solver for dense matricies", false);

} // namespace blue_sky

