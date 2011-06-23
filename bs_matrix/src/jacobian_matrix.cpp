/**
 * \file jacobian_matrix.cpp
 * \brief impl of jacobian matrix
 * \author Sergey Miryanov
 * \date 11.06.2008
 * */
#include "bs_matrix_stdafx.h"

#include "strategies.h"
#include "jacobian_matrix.h"
#include "shared_vector.h"
#include "naive_file_reader.h"
#include "read_b_matrix.h"
//#include "mesh_rs.h"

#include "save_seq_vector.h"

namespace blue_sky
  {

  //! constructor
  jacobian_matrix::jacobian_matrix (bs_type_ctor_param /*param*/)
      : regular_matrix (BS_KERNEL.create_object (csr_matrix_t::bs_type ()))
      , irregular_matrix (BS_KERNEL.create_object (csr_matrix_t::bs_type ()))
      , merged_reg_irreg_matrix (BS_KERNEL.create_object (csr_matrix_t::bs_type ()))
  {

  }

  //! copy constructor
  jacobian_matrix::jacobian_matrix (const jacobian_matrix &matrix) 
        : bs_refcounter (matrix), base_t (matrix)
  {
    BS_ASSERT (false && "TEST ME");

    if (this != &matrix)
      {
        *(regular_matrix->lock ())        = *matrix.regular_matrix;
        *(irregular_matrix)               = *matrix.irregular_matrix;

        solution.assign                   (matrix.solution.begin (), matrix.solution.end ());
        rhs.assign                        (matrix.rhs.begin (), matrix.rhs.end ());
        rhs_flux.assign                   (matrix.rhs_flux.begin (), matrix.rhs_flux.end ());
        regular_accumulative_diag.assign  (matrix.regular_accumulative_diag.begin (), matrix.regular_accumulative_diag.end ());
      }

  }

  /*!
  \brief initialize matrix
  \param N_blocks -- matrix size in blocks
  \param N_block_size -- block size
  \param N_of_diag_bands -- number of diagonal bands
  \return 0 if success
  */
  void
  jacobian_matrix::init (int N_blocks, int N_block_size, int N_of_diag_bands,
                                      const int /*flag_is_line_search*/, const int n_sec)
  {
    BS_ASSERT (N_blocks >= 1 && N_block_size >= 1 && N_of_diag_bands >= 0) (N_blocks) (N_block_size) (N_of_diag_bands);
    if (N_blocks < 1 || N_block_size < 1 || N_of_diag_bands < 0)
      throw bs_exception ("invalid input params", "");

    base_t::n_rows = N_blocks;
    base_t::n_cols = base_t::n_rows; // TODO: BUG
    base_t::n_block_size = N_block_size;

#ifdef _MPI
    BS_ASSERT(false&&"MPI init in jacobian matrix hasn't been iplemented yet");
    //mpi_sol.init (row_starts, comm, n_block_size);
    //mpi_rhs.init (row_starts, comm, n_block_size);
    //mpi_flux_rhs.init (row_starts, comm, n_block_size);

    //solution = mpi_sol.get_local_part ();
    //rhs = mpi_rhs.get_local_part ();
    //flux_rhs = mpi_flux_rhs.get_local_part ();
#else //_MPI
    assign (rhs, N_blocks * N_block_size, 0);
    assign (rhs_flux, N_blocks * N_block_size, 0);
    assign (regular_accumulative_diag, N_blocks * N_block_size * N_block_size, 0);

#endif //_MPI

    // TODO: BUG
    if (sec_rhs.size () == 0)
      {
        if (n_sec > 0)
          {
            n_sec_vars = n_sec;
            assign (sec_rhs, n_sec * N_blocks, 0);
            assign (sec_solution, n_sec * N_blocks, 0);
            assign (ss_diagonal, n_sec * n_sec * N_blocks, 0);
            assign (sp_diagonal, n_sec * N_block_size * N_blocks, 0);
          }
        else
          {
            assign (sec_rhs, 1, 0);
            assign (sec_solution, 1, 0);
            assign (ss_diagonal, 1, 0);
            assign (sp_diagonal, 1, 0);
          }
      }

    const sp_csr_matrix_t &locked_irreg_mx (irregular_matrix);
    const sp_csr_matrix_t &locked_reg_mx (regular_matrix);

    assign (locked_reg_mx->get_values (), locked_reg_mx->get_values ().size (), 0);
    if (locked_reg_mx->get_rows_ptr ().size ())
      {
        item_t value = 0;
        if (locked_irreg_mx->get_rows_ptr ().size ())
          value = locked_irreg_mx->get_rows_ptr ().back ();

        locked_irreg_mx->get_rows_ptr ().resize (locked_reg_mx->get_rows_ptr ().size (), value);
      }
    assign (locked_irreg_mx->get_values (), locked_reg_mx->get_values ().size (), 0);
  }


  //template <class strategy_t>
  //void
  //jacobian_matrix::init_regular_matrix(int N_block_size, const sp_mesh_t &mesh)
  //{
  //  const sp_mesh_t &locked_mesh (mesh);

  //  typename mesh_t::i_vector_type boundary_array;
  //  typedef typename mesh_t::flux_conn_t   flux_conn_t;
  //  typename mesh_t::sp_flux_conn flux_conn = BS_KERNEL.create_object(flux_conn_t::bs_type(), true);
  //
  //  locked_mesh->build_jacobian_and_flux_connections (N_block_size, this->get_regular_matrix(), flux_conn, boundary_array);
  //  assign (regular_matrix->get_values(), N_block_size * N_block_size * (2* locked_mesh->get_n_connections() + locked_mesh->get_n_active_elements()), item_t (0));
  //  m_array = flux_conn->matrix_block_idx_minus;
  //  p_array = flux_conn->matrix_block_idx_plus;
  //  trns_matrix = flux_conn->conn_trans;
  //}

  //template <typename strategy_t>
  //void
  //jacobian_matrix::init_irregular_matrix (index_t N_blocks, index_t N_block_size)
  //{
  //  // TODO: throw exception
  //  irregular_matrix->init (N_blocks, N_blocks, N_block_size, 0);
  //}


  void
  jacobian_matrix::read_regular_from_csr_matrix (const char *filename)
  {
    regular_matrix->ascii_read_from_csr_format (filename);
  }

  void
  jacobian_matrix::read_regular_from_b_matrix (const char *filename, bool summ_diag)
  {
    regular_matrix = tools::read_b_matrix_from_file <csr_matrix_t> (filename, summ_diag);
  }

  void
  jacobian_matrix::read_rhs_vector (const char *filename)
  {
    naive_file_reader (filename)
    .locate_section ("RIGHT HAND SIDE")
    .read_list (rhs);
  }

  void
  jacobian_matrix::read_rhs_flux_vector (const char *filename)
  {
    naive_file_reader (filename)
    .locate_section ("RIGHT HAND SIDE FLUX")
    .read_list (rhs_flux);
  }

  void
  jacobian_matrix::read_solution_vector (const char *filename)
  {
    naive_file_reader (filename)
    .read_list (solution);
  }

  void
  jacobian_matrix::read_acc_diag_from_b_matrix (const char *filename)
  {
    naive_file_reader (filename)
    .locate_section ("MAIN DIAGONAL ACCUMULATIVE")
    .read_list (regular_accumulative_diag);
  }

  void
  jacobian_matrix::read_irregular_csr_matrix (const char *filename)
  {
    BS_ASSERT (!irregular_matrix);

    irregular_matrix = BS_KERNEL.create_object (csr_matrix_t::bs_type ());
    irregular_matrix->ascii_read_from_csr_format (filename);
  }

  const jacobian_matrix::sp_csr_matrix_t &
  jacobian_matrix::get_regular_matrix () const
  {
    return regular_matrix;
  }
  jacobian_matrix::sp_csr_matrix_t
  jacobian_matrix::get_regular_matrix ()
  {
    return regular_matrix;
  }

  const jacobian_matrix::sp_csr_matrix_t &
  jacobian_matrix::get_irregular_matrix () const
  {
    return irregular_matrix;
  }
  jacobian_matrix::sp_csr_matrix_t
  jacobian_matrix::get_irregular_matrix ()
  {
    return irregular_matrix;
  }

  const jacobian_matrix::item_array_t &
  jacobian_matrix::get_solution () const
    {
      return solution;
    }
  jacobian_matrix::item_array_t &
  jacobian_matrix::get_solution ()
  {
    return solution;
  }

  const jacobian_matrix::rhs_item_array_t &
  jacobian_matrix::get_sec_rhs () const
    {
      return sec_rhs;
    }
  jacobian_matrix::rhs_item_array_t &
  jacobian_matrix::get_sec_rhs ()
  {
    return sec_rhs;
  }

  const jacobian_matrix::item_array_t &
  jacobian_matrix::get_sec_solution () const
    {
      return sec_solution;
    }
  jacobian_matrix::item_array_t &
  jacobian_matrix::get_sec_solution ()
  {
    return sec_solution;
  }

  const jacobian_matrix::rhs_item_array_t &
  jacobian_matrix::get_ss_diagonal () const
    {
      return ss_diagonal;
    }
  jacobian_matrix::rhs_item_array_t &
  jacobian_matrix::get_ss_diagonal ()
  {
    return ss_diagonal;
  }

  const jacobian_matrix::rhs_item_array_t &
  jacobian_matrix::get_sp_diagonal () const
    {
      return sp_diagonal;
    }
  jacobian_matrix::rhs_item_array_t &
  jacobian_matrix::get_sp_diagonal ()
  {
    return sp_diagonal;
  }


  const jacobian_matrix::rhs_item_array_t &
  jacobian_matrix::get_rhs () const
    {
      return rhs;
    }
  jacobian_matrix::rhs_item_array_t &
  jacobian_matrix::get_rhs ()
  {
    return rhs;
  }

  const jacobian_matrix::rhs_item_array_t &
  jacobian_matrix::get_rhs_flux () const
    {
      return rhs_flux;
    }
  jacobian_matrix::rhs_item_array_t &
  jacobian_matrix::get_rhs_flux ()
  {
    return rhs_flux;
  }

  const jacobian_matrix::rhs_item_array_t &
  jacobian_matrix::get_regular_acc_diag () const
    {
      return regular_accumulative_diag;
    }
  jacobian_matrix::rhs_item_array_t &
  jacobian_matrix::get_regular_acc_diag ()
  {
    return regular_accumulative_diag;
  }

  //! prepare matrix
  void
  jacobian_matrix::prepare_matrix ()
  {
    BS_ASSERT (merged_reg_irreg_matrix);
    if (!merged_reg_irreg_matrix)
      throw bs_exception ("jacobian_matrix::prepare_matrix", "merged_reg_irreg_matrix is null");

    const sp_csr_matrix_t &lmx (merged_reg_irreg_matrix);

    if (lmx->init (1, 1, regular_matrix->n_block_size, 1))
      throw bs_exception ("jacobian_matrix::prepare_matrix", "can't init merged_reg_irreg_matrix");

    // here we want to use usual merge (not ACPR) thats why we set key=0
    if (lmx->merge (regular_matrix, irregular_matrix, 0, cfl_vector_))
      throw bs_exception ("jacobian_matrix::prepare_matrix", "can't merge reg and irreg matrixes");

    if (lmx->init_diag_ind ())
      throw bs_exception ("jacobian_matrix::prepare_matrix", "can't init merged_reg_irreg_matrix's diag_ind");

    index_t b_sqr = regular_matrix->n_block_size;
    b_sqr = b_sqr * b_sqr;
    rhs_item_array_t &values = lmx->get_values ();

	prepared_values_wo_acc_diag_.assign (values.begin (), values.end ());

    const index_array_t &rows = lmx->get_rows_ptr ();
    const index_array_t &cols = lmx->get_cols_ind ();
    for (index_t i = 0, n_rows = regular_matrix->n_rows; i < n_rows; ++i)
      {
        index_t l = rows[i];
        BS_ASSERT (cols[l] == i) (cols[l]) (i);
        for (index_t k = 0; k < b_sqr; ++k)
          {
            values[l * b_sqr + k] += regular_accumulative_diag[i * b_sqr + k];
          }
      }

    lmx->n_rows       = regular_matrix->n_rows;
    lmx->n_block_size = regular_matrix->n_block_size;
    lmx->n_cols       = regular_matrix->n_cols;
  }

  const strategy_t::rhs_item_array_t &
  jacobian_matrix::get_prepared_values_wo_acc_diag () const
  {
    return prepared_values_wo_acc_diag_;
  }

  //! return merged BCSR matrix
  jacobian_matrix::sp_matrix_t
  jacobian_matrix::get_prepared_matrix () const
  {
    return merged_reg_irreg_matrix;
  }
  //! return merged BCSR matrix
  jacobian_matrix::sp_csr_matrix_t
  jacobian_matrix::get_prepared_bcsr_matrix () const
  {
    return merged_reg_irreg_matrix;
  }


  void
  jacobian_matrix::summ_diag ()
  {
    index_t b_sqr = regular_matrix->n_block_size;
    b_sqr = b_sqr * b_sqr;
    rhs_item_array_t &values = regular_matrix->get_values ();
    const index_array_t &rows = regular_matrix->get_rows_ptr ();
    const index_array_t &cols = regular_matrix->get_cols_ind ();
    for (index_t i = 0, n_rows = regular_matrix->n_rows; i < n_rows; ++i)
      {
        index_t l = rows[i];
        BS_ASSERT (cols[l] == i) (cols[l]) (i);
        for (index_t k = 0; k < b_sqr; ++k)
          {
            values[l * b_sqr + k] += regular_accumulative_diag[i * b_sqr + k];
          }
      }
  }

  /**
  * @brief restore solution for secondary variables
  *        Xs = Dss * Bs - Dss * Asp * Xp
  *        Xs -- (OUTPUT) solution vector for secondary variables
  *        Dss -- (Ass)^(-1)
  *        Bs -- rhs vector for secondary variables
  *        Xp -- solution vector for primary variables
  *
  * @return 0 if success
  */
  void
  jacobian_matrix::restore_sec_solution ()
  {
    rhs_item_t *sp_block = 0;
    item_t *sol_block = 0;

    if (n_sec_vars == 1 && regular_matrix->get_rows_ptr ().size ())
      {
        index_t nb = regular_matrix->n_block_size;
        for (index_t i = 0, cnt = (index_t)regular_matrix->get_rows_ptr ().size (); i < cnt - 1; ++i)
          {
            sp_block = &sp_diagonal[n_sec_vars * nb * i];
            sol_block = &solution[nb * i];
            sec_solution[i] = sec_rhs[i];
            item_t &sec_sol = sec_solution[i];

            VV_PROD_M (nb, sp_block, sol_block, sec_sol);
          }
      }
  }

  // multiply flux part of matrix
  int
  jacobian_matrix::mult_flux_part (const item_t mult)
  {
    static int coef_count = 0;
    ++coef_count;

    const sp_csr_matrix_t &l_irr_mat (irregular_matrix);
    const sp_csr_matrix_t &l_reg_mat (regular_matrix);

    rhs_item_array_t &irr_values = l_irr_mat->get_values ();
    rhs_item_array_t &reg_values = l_reg_mat->get_values ();
    for (size_t i = 0, cnt = irr_values.size (); i < cnt; ++i)
      {
        irr_values[i] *= mult;
      }
    for (size_t i = 0, cnt = reg_values.size (); i < cnt; ++i)
      {
        reg_values[i] *= mult;
      }

    for (size_t i = 0, cnt = (index_t)rhs_flux.size (); i < cnt; ++i)
      {
        rhs_flux[i] *= mult;
      }

    return coef_count;
  }

  void
  jacobian_matrix::reset ()
  {
    BS_ASSERT (this->n_block_size);
    BS_ASSERT (this->n_rows);

    assign (solution, this->n_block_size * this->n_rows, 0);
    assign (rhs, this->n_block_size * this->n_rows, 0);
    assign (rhs_flux, this->n_block_size * this->n_rows, 0);
    assign (regular_accumulative_diag, this->n_block_size * this->n_rows, 0);
  }

  void
  jacobian_matrix::clear_solution ()
  {
    BS_ASSERT (solution.size () >= (size_t)base_t::n_block_size * base_t::n_rows) (solution.size ()) (base_t::n_block_size) (base_t::n_rows);
    assign (solution, solution.size (), 0);
  }

  void
  jacobian_matrix::summ_rhs ()
  {
    BS_ASSERT (rhs.size () >= (size_t)base_t::n_block_size * base_t::n_rows) (rhs.size ()) (base_t::n_block_size) (base_t::n_rows);
    BS_ASSERT (rhs_flux.size () >= (size_t)base_t::n_block_size * base_t::n_rows) (rhs_flux.size ()) (base_t::n_block_size) (base_t::n_rows);
    BS_ASSERT (rhs.size () == rhs_flux.size ()) (rhs.size ()) (rhs_flux.size ());

    size_t i = 0, cnt = rhs.size ();
    for (size_t cnt1 = cnt - (cnt % 4); i < cnt1; i+=4)
      {
        rhs[i] += rhs_flux[i];
        rhs[i+1] += rhs_flux[i+1];
        rhs[i+2] += rhs_flux[i+2];
        rhs[i+3] += rhs_flux[i+3];
      }

    for (; i < cnt; ++i)
      {
        rhs[i] += rhs_flux[i];
      }
  }

  template <typename diag_t, typename v_t, typename r_t>
  void
  matrix_vector_product_free (size_t n_block_size, size_t n_rows, const diag_t &diag, const v_t &v, r_t &r)
  {
    size_t b_sqr = n_block_size * n_block_size;
    for (size_t i = 0, idx = 0, diag_idx = 0; i < n_rows; ++i, idx += n_block_size, diag_idx += b_sqr)
      {
        const typename diag_t::value_type *m_block = &diag[diag_idx];
        const typename v_t::value_type    *v_block = &v[idx];
        typename r_t::value_type          *r_block = &r[idx];
        MV_PROD (n_block_size, m_block, v_block, r_block);
      }
  }

  int
  jacobian_matrix::matrix_vector_product (const double_array_t &v, float_array_t &r) const
  {
    //matrix_vector_product_free (regular_matrix->n_block_size, regular_matrix->n_rows, regular_diag, v, r);
    regular_matrix->matrix_vector_product (v, r);
    irregular_matrix->matrix_vector_product (v, r);

    return 0;
  }
  int
  jacobian_matrix::matrix_vector_product (const float_array_t &v, float_array_t &r) const
  {
    //matrix_vector_product_free (regular_matrix->n_block_size, regular_matrix->n_rows, regular_diag, v, r);
    regular_matrix->matrix_vector_product (v, r);
    irregular_matrix->matrix_vector_product (v, r);

    return 0;
  }
  int
  jacobian_matrix::matrix_vector_product (const double_array_t &v, double_array_t &r) const
  {
    //matrix_vector_product_free (regular_matrix->n_block_size, regular_matrix->n_rows, regular_diag, v, r);
    regular_matrix->matrix_vector_product (v, r);
    irregular_matrix->matrix_vector_product (v, r);

    return 0;
  }

  jacobian_matrix::rhs_item_array_t &
  jacobian_matrix::get_cfl_vector ()
  {
    return cfl_vector_;
  }

  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE (jacobian_matrix);
  BLUE_SKY_TYPE_STD_COPY (jacobian_matrix);

  BLUE_SKY_TYPE_IMPL (jacobian_matrix, strategy_t::matrix_t, "jacobian_matrix", "jacobian_matrix", "jacobian_matrix");


  //////////////////////////////////////////////////////////////////////////
  bool
  jacobian_matrix_register_type (const blue_sky::plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, jacobian_matrix::bs_type ()); BS_ASSERT (res);

    return res;
  }

} // namespace blue_sky
