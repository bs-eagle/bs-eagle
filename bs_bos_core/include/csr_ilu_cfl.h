/**
 *       \file  csr_ilu_cfl.h
 *      \brief  Builds matrix based on CFL for csr_ilu preconditioner
 *              and solves this matrix
 *     \author  Elmira Salimgareeva
 *       \date  11.05.2009
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef BS_CSR_ILU_CFL_PREC_H_
#define BS_CSR_ILU_CFL_PREC_H_


#include BS_FORCE_PLUGIN_IMPORT ()
#include "lsolver_iface.h"
#include "csr_ilu_prec.h"
#include BS_STOP_PLUGIN_IMPORT ()

#include "calc_model.h"


namespace blue_sky
{

  /**
   * \class csr_ilu_cfl_prec
   * \brief Builds matrix based on CFL for csr_ilu prec
   *        and solves this matrix
   * */
  template <class strategy_t>
  class BS_API_PLUGIN csr_ilu_cfl_prec: public csr_ilu_prec<strategy_t>
    {
      //-------------------------------------
      // TYPES
      //-------------------------------------

    public:
      typedef typename strategy_t::matrix_t                 matrix_t;        //< short name for matrix type
      typedef typename strategy_t::item_array_t             item_array_t;    //< short name for array
      typedef typename strategy_t::rhs_item_array_t         rhs_item_array_t;
      typedef typename strategy_t::item_t                   item_t;          //< short name for array item type
      typedef typename strategy_t::rhs_item_t               rhs_item_t;
      typedef typename strategy_t::index_t                  index_t;         //< short name for matrix'es index type
      typedef typename strategy_t::index_array_t            index_array_t;
      typedef calc_model <strategy_t>                       calc_model_t;
      typedef jacobian_matrix <strategy_t>                  jmatrix_t;
      typedef rs_mesh_iface <strategy_t>                    mesh_iface_t;
      typedef typename calc_model_t::data_t                 data_t;
      typedef typename calc_model_t::data_array_t           data_array_t;

      typedef csr_ilu_prec <strategy_t>                     csr_ilu_prec_t;

    private:
      typedef bcsr_matrix<rhs_item_array_t, index_array_t>  bcsr_matrix_t;    //< short name for used matrix
      typedef smart_ptr<bcsr_matrix_t, true>                sp_bcsr_matrix_t; //< short name for smart pointer for matrix

      //------------------------------------
      // METHODS
      //------------------------------------
    public:

      //! destructor
      ~csr_ilu_cfl_prec ();

      //! solve
      virtual int 
      solve (matrix_t *matrix, rhs_item_array_t &rhs, item_array_t &sol);

      virtual int 
      solve_prec (matrix_t *matrix, item_array_t &rhs, item_array_t &sol);

      //! setup via merging matrices extracted from matrix
      virtual int 
      setup (matrix_t *matrix);

    private:

      sp_bcsr_matrix_t create_merged_matrix_for_ilu (jmatrix_t *matrix, const rhs_item_array_t &cfl_vector_);

    public:
      BLUE_SKY_TYPE_DECL (csr_ilu_cfl_prec);


      //---------------------------------------
      // VARIABLES
      //---------------------------------------
    private:
      sp_bcsr_matrix_t  sp_ilu_cfl_;

  };

}


#endif // #ifdef BS_CSR_ILU_CFL_PREC_H_
