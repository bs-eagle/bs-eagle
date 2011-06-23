/**
* @file jacobian.h
* @brief implemetation of jacobian matrix
* @author Morozov Andrew
* @date 2008-05-29
*/

#ifndef JACOBIAN_MATRIX_H_
#define JACOBIAN_MATRIX_H_

#include "bcsr_matrix.h"

namespace blue_sky
{
  class BS_API_PLUGIN jacobian_matrix : public strategy_t::matrix_t, public setup_preconditioner <strategy_t::matrix_t>
    {
      //////////////////////////////////////////////////////////////////////////
      // TYPES
      //////////////////////////////////////////////////////////////////////////
    public:

      typedef strategy_t::matrix_t                   matrix_t;
      typedef strategy_t::csr_matrix_t               csr_matrix_t;

      typedef strategy_t::item_array_t               item_array_t;
      typedef strategy_t::rhs_item_array_t           rhs_item_array_t;
      typedef strategy_t::index_array_t              index_array_t;

      typedef strategy_t::index_t                    index_t;
      typedef strategy_t::item_t                     item_t;
      typedef strategy_t::rhs_item_t                 rhs_item_t;

      typedef item_array_t::array <float>::type   float_array_t;
      typedef item_array_t::array <double>::type  double_array_t;

      typedef bcsr_matrix<item_array_t, index_array_t>        csr_matrix_2_t;       ///< short name for csr_matrix type (it's fail way today)
      typedef smart_ptr <csr_matrix_2_t, true>                sp_csr_matrix_2_t;

      typedef smart_ptr <matrix_t, true>                      sp_matrix_t;
      typedef smart_ptr <csr_matrix_t, true>                  sp_csr_matrix_t;

      typedef matrix_t                                        base_t;
      typedef jacobian_matrix this_t;
      typedef matrix_t                                        matrix_base_t;

      //////////////////////////////////////////////////////////////////////////
      // METHODS
      //////////////////////////////////////////////////////////////////////////
    public:
      BLUE_SKY_TYPE_DECL (jacobian_matrix);

      // initialize matrix, return 0 if success
      void init (int N_blocks, int N_block_size, int N_of_diag_bands, const int, const int n_sec);

      //! read regular matrix from csr file
      void read_regular_from_csr_matrix (const char *filename);

      //! read regular matrix from b_matrix file
      void read_regular_from_b_matrix (const char *filename, bool summ_diag);

      //! read rhs vector from section
      void read_rhs_vector (const char *filename);

      //! read rhs_flux vector from section
      void read_rhs_flux_vector (const char *filename);

      //! read solution vector form file
      void read_solution_vector (const char *filename);

      //! read accumulative diag of regular matrix from b_matrix file
      void read_acc_diag_from_b_matrix (const char *filename);

      //! read irregular matrix from csr file
      void read_irregular_csr_matrix (const char *filename);

      const sp_csr_matrix_t &get_regular_matrix () const;
      sp_csr_matrix_t get_regular_matrix ();

      const sp_csr_matrix_t &get_irregular_matrix () const;
      sp_csr_matrix_t get_irregular_matrix ();

      const item_array_t &get_solution () const;
      item_array_t &get_solution ();

      //! return pointer for the sec_rhs
      const rhs_item_array_t &get_sec_rhs () const;
      rhs_item_array_t &get_sec_rhs ();

      //! return pointer for the sec_solution
      const item_array_t &get_sec_solution () const;
      item_array_t &get_sec_solution ();

      //! return pointer for the ss_diagonal
      const rhs_item_array_t &get_ss_diagonal () const;
      rhs_item_array_t &get_ss_diagonal ();

      //! return pointer for the sp_diagonal
      const rhs_item_array_t &get_sp_diagonal () const;
      rhs_item_array_t &get_sp_diagonal ();

      const rhs_item_array_t &get_rhs () const;
      rhs_item_array_t &get_rhs ();

      const rhs_item_array_t &get_rhs_flux () const;
      rhs_item_array_t &get_rhs_flux ();

      const rhs_item_array_t &get_regular_acc_diag () const;
      rhs_item_array_t &get_regular_acc_diag ();

      //! prepare matrix
      void prepare_matrix ();

      void
      summ_main_acc_diag ();

      //! return merged BCSR matrix
      sp_matrix_t get_prepared_matrix () const;

      //! same as previous but with out upcast
      sp_csr_matrix_t get_prepared_bcsr_matrix () const;

      const rhs_item_array_t &
      get_prepared_values_wo_acc_diag () const;

      rhs_item_array_t &
      get_cfl_vector ();

      void restore_sec_solution ();

      void summ_diag ();

      int mult_flux_part (const item_t mult);

      void reset ();

      void clear_solution ();
      void summ_rhs ();

      virtual int matrix_vector_product (const double_array_t &v, float_array_t &r) const;
      virtual int matrix_vector_product (const float_array_t &v, float_array_t &r) const;
      virtual int matrix_vector_product (const double_array_t &v, double_array_t &r) const;

      //////////////////////////////////////////////////////////////////////////
      // VARIABLES
      //////////////////////////////////////////////////////////////////////////
      index_array_t         m_array;
      index_array_t         p_array;
      sp_csr_matrix_t       trns_matrix;

    private:

      sp_csr_matrix_t       regular_matrix;
      sp_csr_matrix_t       irregular_matrix;
      sp_csr_matrix_t       merged_reg_irreg_matrix;
      item_array_t          solution;

      rhs_item_array_t      rhs;
      rhs_item_array_t      rhs_flux;
      rhs_item_array_t      regular_accumulative_diag;
      rhs_item_array_t      prepared_values_wo_acc_diag_;

      rhs_item_array_t      sec_rhs;                            //!< [n_sec_vars * n_rows] right hand side vector for secondary variables
      item_array_t          sec_solution;                       //!< [n_sec_vars * n_rows] solution vector for secondary variables
      rhs_item_array_t      ss_diagonal;                        //!< [n_sec_vars * n_sec_vars * n_rows] SS part of matrix diagonal
      rhs_item_array_t      sp_diagonal;                        //!< [n_sec_vars * n_block_size * n_rows] SP part of matrix diagonal

      rhs_item_array_t      cfl_vector_;

      index_t               n_sec_vars;
    };

  bool
  jacobian_matrix_register_type (const blue_sky::plugin_descriptor &pd);
 
} // namespace blue_sky
#endif // JACOBIAN_MATRIX_H_

