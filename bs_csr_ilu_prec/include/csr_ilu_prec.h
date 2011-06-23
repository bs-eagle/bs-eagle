/*!
 * \file csr_ilu_prec.h
 * \brief class declaration for ILU preconditioner for CSR matrix
 * \author Borschuk Oleg
 * \date 2006-11-07
 */
#ifndef CSR_ILU__PREC__H__
#define CSR_ILU__PREC__H__

#include BS_FORCE_PLUGIN_IMPORT ()
#include "linear_solvers.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {

  /**
   * \brief ILU preconditioner for CSR matrix
   */
  class BS_API_PLUGIN csr_ilu_prec: public linear_solver_base
    {
      //-----------------------------------------
      // TYPES
      //-----------------------------------------
    public:
      typedef strategy_t::matrix_t       matrix_t;       ///< short name to matrix type
      typedef strategy_t::item_array_t   item_array_t;   ///< short name to array type
      typedef strategy_t::index_array_t  index_array_t;
      typedef strategy_t::item_t         item_t;         ///< short name to array item type
      typedef strategy_t::index_t        index_t;        ///< short name to matrix's index type

      typedef strategy_t::rhs_item_t	  rhs_item_t;     ///< short name for type of rhs
      typedef strategy_t::rhs_item_array_t	rhs_item_array_t;   ///< short name for rhs array type

      typedef linear_solver_base      this_t;         ///< typedef to this type
      typedef linear_solver_base      base_t;         ///< typedef to this type. in child classes used as a short name of base class
      typedef smart_ptr<this_t, true>             sp_this_t;      ///< short name to smart pointer to this class
      typedef smart_ptr<linear_solver_prop, true> sp_prop_t;      ///< short name to smart pointer to properties holder class

      typedef smart_ptr<matrix_t, true>           sp_matrix_t;    ///< short name to smart pointer on matrix

      //-----------------------------------------
      //  TYPES
      //-----------------------------------------
    private:
      typedef bcsr_matrix<rhs_item_array_t, index_array_t>        bcsr_matrix_t;        ///< short name for used matrix
      typedef smart_ptr<bcsr_matrix_t, true>                      sp_bcsr_matrix_t;     ///< short name for smart_pointer on used matrix
      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
    public:
      //! destructor
      ~csr_ilu_prec ();

      //! solve
      virtual int solve (matrix_t *matrix, rhs_item_array_t &rhs, item_array_t &sol);

      virtual int solve_prec (matrix_t *matrix, item_array_t &rhs, item_array_t &sol);

      //! setup via merging matrices extracted from matrix
      virtual int setup (matrix_t *matrix);


      sp_bcsr_matrix_t get_ilu_matrix () const;

    protected:
      //! setup
      int setup_internal (bcsr_matrix_t *matrix);

      template <class rhs_t>
      int
      templ_solve (matrix_t *matrix, rhs_t &rhs, item_array_t &sol);


    private:

      /**
       * \brief Check internal state of preconditioner

       * \return If ILU matrix is valid return True
       */
      bool check_state_internal ();

    public:
      BLUE_SKY_TYPE_DECL (csr_ilu_prec);
      //-----------------------------------------
      //  VARIABLES
      //-----------------------------------------
    public:

    private:

      sp_bcsr_matrix_t sp_ilu;                    //!< storage for ILU matrix
    };

  //! register types in kernel
  bool csr_ilu_prec_register_type (const blue_sky::plugin_descriptor &pd);

} // namespace blue_sky
#endif // CSR_ILU__PREC__H__
