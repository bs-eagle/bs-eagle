/**
 *       \file  jacobian.h
 *      \brief  Implementation of Jacobian matrix
 *     \author  Morozov Andrey
 *       \date  29.05.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef JACOBIAN_H_
#define JACOBIAN_H_

#include "fi_params.h"
#include BS_FORCE_PLUGIN_IMPORT ()
#include "mbcsr_matrix_iface.h"
#include "bdiag_matrix_iface.h"
#include "lsolver_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
{

  /**
   * \class jacobian
   * \brief Incapsulates solving of jacobian matrix
   * */
  class BS_API_PLUGIN jacobian : public bs_node
    {
    public:

      //! initialize jacobian
      void
      init (t_long elements, t_long phases, t_long secondary);

      /**
       * \brief returns matrix by name
       * \return returns matrix or throw exception if no matrix with name
       * */
      BS_SP (bcsr_matrix_iface)
      get_matrix (std::string const &name) const;

      /**
       * \brief returns mbcsr matrix
       * */
      BS_SP (mbcsr_matrix_iface)
      get_matrix () const;

      //! set solver
      void 
      set_solver(const BS_SP (lsolver_iface) &s)
      {
        solver = s;
      }
      //! set preconditioner
      void 
      set_prec (const BS_SP (lsolver_iface) &p)
      {
        preconditioner = p;
      }

      // FIXME: where this vector should be?
      virtual spv_double
      get_solution () const;

      // FIXME: where this vector should be?
      virtual spv_double
      get_sec_solution () const;

      //! get solver
      const BS_SP (lsolver_iface) &
      get_solver () const;
      
      //! get preconditioner
      const BS_SP (lsolver_iface) &
      get_prec () const;

      //! set up params for solver
      int 
      setup_solver_params (well_model_type model_type, int n_phases, const BS_SP (fi_params) &ts_params);

      //! run solver setup
      int 
      setup ();

      //! solve
      t_double 
      solve (t_long &n_lin_iters);

      // FIXME:
      void clear_solution ();
      void summ_rhs ();
      void mult_flux_part (t_double dt_mult);
      void prepare_matrix ();
      void restore_sec_solution ();

      BS_SP (bdiag_matrix_iface) get_accum_matrix ();
      spv_float get_ss_diagonal ();
      spv_float get_sp_diagonal ();
      spv_float get_sec_rhs ();
      spv_float get_rhs ();
      spv_double get_solution ();
      spv_double get_cfl_vector ();

      BS_SP (bcsr_matrix_iface) get_prepared_matrix ();

    public:
      virtual ~jacobian();

      //! blue-sky type description
      BLUE_SKY_TYPE_DECL (jacobian);

    private:
      /**
       * \brief  Setups solver params
       * \param  ts_params Params holder
       * */
      void 
      setup_solver (const BS_SP (fi_params) &ts_params);

      //! \todo Obsolete
      void 
      setup_preconditioner (const BS_SP (fi_params) &ts_params);

    protected:

      BS_SP (mbcsr_matrix_iface)  matrix;
      BS_SP (lsolver_iface)       solver;                   //!< solver
      BS_SP (lsolver_iface)       preconditioner;           //!< preconditioner
      auto_value <int>            solver_is_gmres_flag;     //!< if != 0 solver point to the GMRES solver
      auto_value <int>            prec_is_cpr_flag;         //!< if != 0 if using CPR precondition
    };

} // namespace blue_sky
#endif // JACOBIAN_H_
