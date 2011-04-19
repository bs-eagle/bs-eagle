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
#include "jac_matrix_iface.h"
#include "lsolver_iface.h"
#include "rs_mesh_iface.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include "make_me_happy.h"

namespace blue_sky
{

  /**
   * \class jacobian
   * \brief Incapsulates solving of jacobian matrix
   * */
  class BS_API_PLUGIN jacobian : public bs_node
    {
    public:
      typedef matrix_iface                          matrix_t;             //!< short name to matrix type
      typedef spv_double                            item_array_t;         //!< short name to array type
      typedef spv_float                             rhs_item_array_t;
      typedef t_double                              item_t;               //!< short name to array item type
      typedef t_long                                index_t;              //!< short name to matrix's index type
      typedef spv_long                              index_array_t;
      typedef rs_mesh_iface                         mesh_iface_t;
      typedef bcsr_matrix_iface                     bcsr_matrix_t;
      typedef lsolver_iface                         linear_solver_base_t;

      typedef lsolver_iface::sp_matrix_t            sp_matrix_t;
      typedef smart_ptr <bcsr_matrix_t, true>       sp_bcsr_matrix_t;
      typedef smart_ptr <jac_matrix_iface, true>    sp_jmatrix;           //!< matrix type
      typedef smart_ptr <lsolver_iface, true>       sp_lsolver;           //!< solver & preconditioner type
      typedef smart_ptr <fi_params, true>           sp_fi_params;

      typedef smart_ptr <mesh_iface_t, true>        sp_mesh_iface_t;

      typedef jacobian                              this_t;

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

      //! set solver
      void 
      set_solver(const sp_lsolver &s)
      {
        solver = s;
      }
      //! set preconditioner
      void 
      set_prec (const sp_lsolver &p)
      {
        preconditioner = p;
      }
      //! set jacobian matrix
      void
      set_jmatrix (const sp_jmatrix &mx)
      {
        jm = mx;
      }

      //! get jacobian matrix
      const sp_jmatrix &
      get_jmatrix() const
      {
        return jm;
      }

      // TODO: where this vector should be?
      virtual const item_array_t &
      get_solution () const;

      // TODO: where this vector should be?
      virtual const item_array_t &
      get_sec_solution () const;

      //! get solver
      const sp_lsolver &
      get_solver () const;
      
      //! get preconditioner
      const sp_lsolver &
      get_prec () const;

      //! set up params for solver
      int 
      setup_solver_params (well_model_type model_type, int n_phases, const sp_fi_params &ts_params);

      //! run solver setup
      int 
      setup ();

      //! \todo Obsolete
      void 
      begin ();

      //! \todo Obsolete
      void 
      end ();

      //! solve
      item_t 
      solve (index_t &n_lin_iters);

      /**
       * \class jacob_traits
       * \brief For sorting jacobian children
       * */
      struct jacob_traits : bs_node::sort_traits
        {
          struct jacob_key : bs_node::sort_traits::key_type
            {
              virtual bool sort_order (const key_ptr & ) const
                {
                  return true;
                }
            };

          virtual const char * sort_name () const
            {
              return "jacob trait";
            };

          virtual key_ptr key_generator (const sp_link& /*l*/) const
            {
              return new jacob_key ();
            }

          virtual bool accepts (const sp_link& l)
          {
            return smart_ptr< this_t, true >(l->data(), bs_dynamic_cast());
          }
        };

    public:
      /**
       * \brief  dtor
       * */
      virtual ~jacobian()
        {
        }

      //! blue-sky type description
      BLUE_SKY_TYPE_DECL (jacobian);

    private:
      /**
       * \brief  Setups solver params
       * \param  ts_params Params holder
       * */
      void 
      setup_solver (const sp_fi_params &ts_params);

      //! \todo Obsolete
      void 
      setup_preconditioner (const sp_fi_params &ts_params);

    protected:

      sp_jmatrix                  jm;                       //!< jacobian matrix
      sp_lsolver                  solver;                   //!< solver
      sp_lsolver                  preconditioner;           //!< preconditioner
      auto_value <int>            solver_is_gmres_flag;     //!< if != 0 solver point to the GMRES solver
      auto_value <int>            prec_is_cpr_flag;         //!< if != 0 if using CPR precondition
    };

} // namespace blue_sky
#endif // JACOBIAN_H_
