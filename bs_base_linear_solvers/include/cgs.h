/**
 * \file cgs.h
 * \brief CGS linear solver
 * \author SalimgareevaEM
 * \date
 * */
#ifndef BS_CGS_LINEAR_SOLVER_H_
#define BS_CGS_LINEAR_SOLVER_H_

#include "linear_solvers.h"

namespace blue_sky
  {
  /**
  * @brief CGS linear solver
  */
  template <class strategy_t>
  class BS_API_PLUGIN cgs_solver : public linear_solver_base<strategy_t>
    {

      //-----------------------------------------
      // TYPES
      //-----------------------------------------
    public:
      typedef typename strategy_t::matrix_t       matrix_t;       ///< short name to matrix type
      typedef typename strategy_t::item_array_t   item_array_t;   ///< short name to array type
      typedef typename strategy_t::item_t         item_t;         ///< short name to array item type
      typedef typename strategy_t::index_t        index_t;        ///< short name to matrix's index type
      typedef typename strategy_t::rhs_item_t	  rhs_item_t;     ///< short name for type of rhs
      typedef typename strategy_t::rhs_item_array_t	rhs_item_array_t;   ///< short name for rhs array type

      typedef typename strategy_t::index_array_t  index_array_t;

      typedef linear_solver_base<strategy_t>      this_t;         ///< typedef to this type
      typedef linear_solver_base<strategy_t>      base_t;         ///< typedef to this type. in child classes used as a short name of base class
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
      virtual ~cgs_solver ();

      //! solve
      virtual int solve (matrix_t *matrix, rhs_item_array_t &rhs, item_array_t &sol);

      virtual int solve_prec (matrix_t *matrix, item_array_t &rhs, item_array_t &sol);

      //! setup
      virtual int setup (matrix_t *matrix);

      //template <class mx_t>
      //int init_by_matrix (const mx_t *mx);

      template <class rhs_t>
      int
      templ_solve (matrix_t *matrix, rhs_t &rhs, item_array_t &sol);

    private:
      //-----------------------------------------
      //  VARIABLES
      //-----------------------------------------
    public:

    public:
      BLUE_SKY_TYPE_DECL (cgs_solver);
    };

  }	// namespace blue_sky

#endif // #ifndef BS_CGS_LINEAR_SOLVER_H_

