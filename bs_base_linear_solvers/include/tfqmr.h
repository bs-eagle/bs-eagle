/**
 * \file tfqmr.h
 * \brief TFQMR linear solver
 * \author Salimgareeva E.M.
 * \date
 * */
#ifndef BS_TFQMR_LINEAR_SOLVER_H_
#define BS_TFQMR_LINEAR_SOLVER_H_

#include "linear_solvers.h"

namespace blue_sky
  {
  /**
  * @brief TFQMR linear solver
  */
  class BS_API_PLUGIN tfqmr_solver : public linear_solver_base
    {

      //-----------------------------------------
      // TYPES
      //-----------------------------------------
    public:
      typedef strategy_t::matrix_t       matrix_t;       ///< short name to matrix type
      typedef strategy_t::item_array_t   item_array_t;   ///< short name to array type
      typedef strategy_t::item_t         item_t;         ///< short name to array item type
      typedef strategy_t::index_t        index_t;        ///< short name to matrix's index type
      typedef strategy_t::rhs_item_t	  rhs_item_t;     ///< short name for type of rhs
      typedef strategy_t::rhs_item_array_t	rhs_item_array_t;   ///< short name for rhs array type

      typedef strategy_t::index_array_t  index_array_t;

      typedef linear_solver_base      this_t;         ///< typedef to this type
      typedef linear_solver_base      base_t;         ///< typedef to this type. in child classes used as a short name of base class
      typedef smart_ptr<this_t, true>             sp_this_t;      ///< short name to smart pointer to this class
      typedef smart_ptr<linear_solver_prop, true> sp_prop_t;      ///< short name to smart pointer to properties holder class

      typedef smart_ptr<matrix_t, true>           sp_matrix_t;    ///< short name to smart pointer on matrix_t

      typedef bcsr_matrix<rhs_item_array_t, index_array_t>    bcsr_matrix_t;        ///< short name for used matrix
      typedef smart_ptr<bcsr_matrix_t, true>                  sp_bcsr_matrix_t;     ///< short name for smart_pointer on used matrix

      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
    public:
      //! destructor
      virtual ~tfqmr_solver ();

      //! solve
      virtual int solve (matrix_t *matrix, rhs_item_array_t &rhs, item_array_t &sol);

      virtual int solve_prec (matrix_t *matrix, item_array_t &rhs, item_array_t &sol);

      template <class rhs_t>
      int
      templ_solve (matrix_t *matrix, rhs_t &rhs, item_array_t &sol);

      //! setup
      virtual int setup (matrix_t *matrix);

    private:
      //-----------------------------------------
      //  VARIABLES
      //-----------------------------------------
    public:

    public:
      BLUE_SKY_TYPE_DECL (tfqmr_solver);
    };

  }	// namespace blue_sky

#endif // #ifndef BS_TFQMR_LINEAR_SOLVER_H_
