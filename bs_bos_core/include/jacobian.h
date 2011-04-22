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

      spv_float get_ss_diagonal ();
      spv_float get_sp_diagonal ();
      spv_float get_sec_rhs ();
      spv_float get_rhs ();
      spv_float get_rhs_flux ();

      spv_double get_cfl_vector ();
      spv_double get_solution ();
      spv_double get_sec_solution ();

      spv_long get_boundary ();
      spv_long get_m_array ();
      spv_long get_p_array ();

      BS_SP (flux_connections_iface)
      get_flux_connections ();

      /**
      * @brief restore solution for secondary variables
      *        Xs = Dss * Bs - Dss * Asp * Xp
      *        Xs -- (OUTPUT) solution vector for secondary variables
      *        Dss -- (Ass)^(-1)
      *        Bs -- rhs vector for secondary variables
      *        Xp -- solution vector for primary variables
      */
      void 
      restore_sec_solution ();

      void clear_solution ();
      void summ_rhs ();
      void mult_flux_part (t_double dt_mult);

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

      BS_SP (mbcsr_matrix_iface)      matrix;
      BS_SP (flux_connections_iface)  flux_conn;
      spv_float                       ss_diagonal;
      spv_float                       sp_diagonal;
      spv_float                       sec_rhs;
      spv_float                       rhs;
      spv_float                       rhs_flux;
      spv_double                      cfl_vector;
      spv_double                      solution;
      spv_double                      sec_solution;
      spv_long                        boundary;
      BS_SP (lsolver_iface)           solver;                   //!< solver
      BS_SP (lsolver_iface)           preconditioner;           //!< preconditioner
      auto_value <int>                solver_is_gmres_flag;     //!< if != 0 solver point to the GMRES solver
      auto_value <int>                prec_is_cpr_flag;         //!< if != 0 if using CPR precondition

      // FIXME: we shoudl store number of secondary 
      // and number of phases variables to properly 
      // do restore_sec_solution
      t_long                          phases;
      t_long                          secondary;
    };

} // namespace blue_sky
#endif // JACOBIAN_H_
