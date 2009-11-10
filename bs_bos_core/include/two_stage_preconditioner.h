#ifndef TWO_STAGE_PREC_H__
#define TWO_STAGE_PREC_H__
/*!
 * \file two_stage_preconditioner.h
 * \brief class declaration for two stage preconditioner
 * \author Borschuk Oleg
 * \date 2006-07-26
 */
#ifdef _MPI
#include "mpi_csr_matrix.h"
#endif //_MPI
#include "memory_macroses.h"

#include BS_FORCE_PLUGIN_IMPORT ()
#include "linear_solvers.h"
#include BS_STOP_PLUGIN_IMPORT ()

namespace blue_sky
  {

  template <class strategy_t>
  class BS_API_PLUGIN two_stage_preconditioner: public linear_solver_base<strategy_t>
    {
    public:
      typedef linear_solver_base<strategy_t>          base_t;         ///< short name for base class
      typedef two_stage_preconditioner<strategy_t>    this_t;         ///< short name for this class
      typedef typename strategy_t::matrix_t           matrix_t;       ///< short name for matrix type
      typedef typename strategy_t::item_array_t       item_array_t;        ///< short name for array type
      typedef typename strategy_t::item_t             item_t;   ///< short name for array item type
      typedef typename strategy_t::rhs_item_t	      rhs_item_t;     ///< short name for type of rhs
      typedef typename strategy_t::rhs_item_array_t	  rhs_item_array_t;   ///< short name for rhs array type
      typedef typename strategy_t::index_t            index_t;        ///< short name for index type
      typedef item_t                                  fp_type;        ///< short name for array item type
      typedef index_t                                 i_type;         ///< short name for index type
      typedef smart_ptr <base_t, true>                sp_base_t;      ///< short name for smart pointer on base class
      typedef smart_ptr <matrix_t, true>              sp_matrix_t;    ///< short name for smart pointer on matrix
      //-----------------------------------------
      //  METHODS
      //-----------------------------------------
    public:
      //! destructor
      virtual ~two_stage_preconditioner ();

      ////! solve preconditioner
      //virtual int solve (matrix_t *matrix, item_array_t &rhs, item_array_t &sol);
      //! solve
      virtual int solve (matrix_t *matrix, rhs_item_array_t &rhs, item_array_t &sol);

      virtual int solve_prec (matrix_t *matrix, item_array_t &rhs, item_array_t &sol);

      template <class rhs_t>
      int
      templ_solve (matrix_t *matrix, rhs_t &rhs, item_array_t &sol);


      //! setup preconditioner
      virtual int setup (matrix_t *matrix);

      //! set first preconditioner
      void set_prec_1 (const base_t *prec)
      {
        base_t::prec = prec;
        set_subnode_in_tree ("prec", sp_base_t (prec));
      }

      //! set second preconditioner
      void set_prec_2 (const base_t *prec)
      {
        prec_2 = prec;
        set_subnode_in_tree ("prec_2", sp_base_t (prec));
      }

      //! get iteration from first preconditioner
      int get_prec_1_iters ()
      {
        return base_t::prec ? base_t::prec->get_prop ()->get_iters () : 0;
      }

      //! get first preconditioner
      sp_base_t get_prec_1 () const
        {
          return base_t::prec;
        }

      //! get second preconditioner
      sp_base_t get_prec_2 () const
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

      BLUE_SKY_TYPE_DECL (two_stage_preconditioner);

    private:
      //-----------------------------------------
      //  VARIABLES
      //-----------------------------------------
      //#ifdef _MPI
      //    mpi_vector <double> mpi_r, mpi_w;
      //#endif //_MPI

    public:

    private:
      sp_base_t       prec_2;

      item_array_t    r_array;
      item_array_t    w_array;
    };

  //! register types into kernel
  bool
  two_stage_prec_register_type (const blue_sky::plugin_descriptor &pd);

/*
  template <class solver_t>
  int
  solve_helper (solver_t *solver, typename solver_t::matrix_t *mx, seq_vector <float> &rhs, typename solver_t::item_array_t &sol)
    {
      solver->solve (mx, rhs, sol);
    }

  template <class solver_t>
  int
  solve_helper (solver_t *solver, typename solver_t::matrix_t *mx, seq_vector <double> &rhs, typename solver_t::item_array_t &sol)
    {
      solver->solve_prec (mx, rhs, sol);
    }

*/
} // namespace blue_sky
#endif // TWO_STAGE_PREC_H__
