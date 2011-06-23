/**
 *       \file  two_stage_preconditioner.h
 *      \brief  Class declaration for two stage preconditioner
 *     \author  Borschuk Oleg
 *       \date  26.07.2006
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef TWO_STAGE_PREC_H__
#define TWO_STAGE_PREC_H__

#ifdef _MPI
#include "mpi_csr_matrix.h"
#endif //_MPI
#include "memory_macroses.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "linear_solvers.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {

  /**
   * \class two_stage_preconditioner
   * \brief Class declaration for two stage preconditioner
   * */
  class BS_API_PLUGIN two_stage_preconditioner: public linear_solver_base
    {
    public:
      typedef linear_solver_base          base_t;             ///< short name for base class
      typedef two_stage_preconditioner    this_t;             ///< short name for this class
      typedef strategy_t::matrix_t           matrix_t;           ///< short name for matrix type
      typedef strategy_t::item_array_t       item_array_t;       ///< short name for array type
      typedef strategy_t::item_t             item_t;             ///< short name for array item type
      typedef strategy_t::rhs_item_t         rhs_item_t;         ///< short name for type of rhs
      typedef strategy_t::rhs_item_array_t   rhs_item_array_t;   ///< short name for rhs array type
      typedef strategy_t::index_t            index_t;            ///< short name for index type
      typedef item_t                                  fp_type;            ///< short name for array item type
      typedef index_t                                 i_type;             ///< short name for index type
      typedef smart_ptr <base_t, true>                sp_base_t;          ///< short name for smart pointer on base class
      typedef smart_ptr <matrix_t, true>              sp_matrix_t;        ///< short name for smart pointer on matrix
      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
    public:
      //! destructor
      virtual ~two_stage_preconditioner ();

      /**
       * \brief  Solves Jacobian matrix
       * \param  matrix Jacobian matrix
       * \param  rhs 
       * \param sol
       * \return 0 on success
       * */
      virtual int 
      solve (matrix_t *matrix, rhs_item_array_t &rhs, item_array_t &sol);

      /**
       * \brief  Solves matrix from sub-preconditioner
       * \param  matrix
       * \param  rhs
       * \param  sol
       * \return 0 on success
       * */
      virtual int 
      solve_prec (matrix_t *matrix, item_array_t &rhs, item_array_t &sol);

      /**
       * \brief  Solves matrix
       * \param  matrix
       * \param  rhs
       * \param  sol
       * \return 0 on success
       * */
      template <class rhs_t>
      int
      templ_solve (matrix_t *matrix, rhs_t &rhs, item_array_t &sol);


      /**
       * \brief  Setups preconditioner to solve matrix matrix
       * \param  matrix
       * \return 0 on success
       * */
      virtual int 
      setup (matrix_t *matrix);

      /**
       * \brief  Sets first preconditioner
       * \param  prec Instance of new preconditioner
       * */
      void 
      set_prec_1 (const base_t *prec)
      {
        base_t::prec = prec;
        set_subnode_in_tree ("prec", sp_base_t (prec));
      }

      /**
       * \brief  Sets second preconditioner
       * \param  prec Instance of new preconditioner
       * */
      void 
      set_prec_2 (const base_t *prec)
      {
        prec_2 = prec;
        set_subnode_in_tree ("prec_2", sp_base_t (prec));
      }

      /**
       * \brief  Returns iterations count from first
       *         preconditioner
       * \return Iterations count
       * */
      int 
      get_prec_1_iters ()
      {
        return base_t::prec ? base_t::prec->get_prop ()->get_iters () : 0;
      }

      /**
       * \brief  Returns first preconditioner
       * \return First preconditioner instance
       * */
      sp_base_t 
      get_prec_1 () const
        {
          return base_t::prec;
        }

      /**
       * \brief  Returns second preconditioner
       * \return Second preconditioner instance
       * */
      sp_base_t 
      get_prec_2 () const
        {
          return prec_2;
        }

      ///**
      // * \brief set sub node by their name
      // *
      // * \param[in] name  The name of subnode
      // * \param[in] obj   The smart pointer to new object
      // */
      //virtual void set_subnode (const std::string &name, sp_obj obj)
      //{
      //  if (name == "prec_2")
      //    {
      //      base_t::set_subnode_internal (name, obj, prec_2);
      //    }
      //  else
      //    {
      //      base_t::set_subnode (name, obj);
      //    }
      //}

      //! blue-sky type declaration
      BLUE_SKY_TYPE_DECL (two_stage_preconditioner);

    private:
      //-----------------------------------------
      //  VARIABLES
      //-----------------------------------------

    private:
      sp_base_t       prec_2;   //!< Second preconditioner

      item_array_t    r_array;  //!< Temporary array
      item_array_t    w_array;  //!< Temporary array
    };

  /**
   * \brief  Registers types in blue-sky kernel
   * \param  pd plugin_descriptor
   * \return True if all types registered successfully
   * */
  bool
  two_stage_prec_register_type (const blue_sky::plugin_descriptor &pd);

} // namespace blue_sky
#endif // TWO_STAGE_PREC_H__
