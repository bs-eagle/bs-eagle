/*!
 * \file csr_ilu_prec.cpp
 * \brief implementation of incomplit LU precondition for CSR matrix
 * \author Borschuk Oleg
 * \date 2006-11-07
 */
#include "bs_csr_ilu_prec_stdafx.h"
#include "csr_ilu_prec.h"


namespace blue_sky
  {
  /*!
   * \brief constructor
   */
  csr_ilu_prec::csr_ilu_prec (bs_type_ctor_param param)
      : linear_solver_base (param)
      , sp_ilu(BS_KERNEL.create_object(bcsr_matrix_t::bs_type()))
  {
  }

  csr_ilu_prec::csr_ilu_prec(const csr_ilu_prec& solver)
      : bs_refcounter (solver), linear_solver_base (solver)
  {
    if (&solver != this)
      *this = solver;
  }

  /*!
   * \brief destructor
   */
  csr_ilu_prec::~csr_ilu_prec ()
  {
  }

  /**
   * \brief check internal state of preconditioner
   */
  bool
  csr_ilu_prec::check_state_internal ()
  {
    BS_ASSERT (sp_ilu);
    if (!sp_ilu)
      {
        return false;
      }
    bcsr_matrix_t &ilu = *sp_ilu;

    BS_ASSERT (!ilu.get_values ().empty ());
    if (ilu.get_values ().empty ())
      return false;
    BS_ASSERT (!ilu.get_rows_ptr ().empty ());
    if (ilu.get_rows_ptr ().empty ())
      return false;
    BS_ASSERT (!ilu.get_cols_ind ().empty ());
    if (ilu.get_cols_ind ().empty ())
      return false;
    BS_ASSERT (!ilu.get_diag_ind ().empty ());
    if (ilu.get_diag_ind ().empty ())
      return false;


    return true;
  }


  int csr_ilu_prec::solve(matrix_t *matrix, rhs_item_array_t &rhs, item_array_t &sol)
  {
    return templ_solve (matrix, rhs, sol);
  }

  int csr_ilu_prec::solve_prec(matrix_t *matrix, item_array_t &rhs, item_array_t &sol)
  {
     return templ_solve (matrix, rhs, sol);
  }
  /*!
   * \brief
   *
   * \param v
   * \param r
   *
   * \return
   */
  //template <class strategy_t> int
  //csr_ilu_prec::solve (matrix_t *matrix, rhs_item_array_t &rhs, item_array_t &sol)
  template <class rhs_t> int
  csr_ilu_prec::templ_solve (matrix_t *matrix, rhs_t &rhs, item_array_t &sol)
  {
    typedef item_t fp_type;

    BS_ASSERT (matrix);
    BS_ASSERT (rhs.size ());
    BS_ASSERT (sol.size ());
    BS_ASSERT (rhs.size () == sol.size ()) (rhs.size ()) (sol.size ());
    BS_ASSERT (sp_ilu);
    BS_ASSERT (!sp_ilu->get_values ().empty ());
    BS_ASSERT (!sp_ilu->get_rows_ptr ().empty ());
    BS_ASSERT (!sp_ilu->get_cols_ind ().empty ());

    bcsr_matrix_t &ilu = *sp_ilu;

    int j, k, l;
    int j1, j2;
    //const fp_type *D_block, *M_block;
    typedef typename strategy_t::rhs_item_t rhs_item_t;
    const rhs_item_t *D_block, *M_block;

    fp_type *v, *r;

    index_t b_sqr     = ilu.n_block_size * ilu.n_block_size;
    index_t n         = ilu.n_rows;
    index_t nb        = ilu.n_block_size;

    //const item_array_t &values    = ilu.get_values   ();
    const rhs_item_array_t &values    = ilu.get_values   ();
    const index_array_t &rows     = ilu.get_rows_ptr ();
    const index_array_t &cols     = ilu.get_cols_ind ();
    //const index_array_t &diag_ind = ilu.get_diag_ind ();

    //sol.assign (rhs.begin (), rhs.end ());
    memcpy (&sol[0], &rhs[0], rhs.size () * sizeof (rhs[0]));

#ifdef _DEBUG
    item_t *sol_ = &sol[0];
    sol_;
#endif
    // solve Ly = b
    r = &sol[0];
    for (k = 0; k < n; ++k, r += nb)
      {
        // pointer to matrix row
        j1          = rows[k];
        j2          = j1;//diag_ind[k];
        index_t j3  = rows[k + 1];
        for (j = j1; j < j3; ++j)
          {
            l = cols[j];
            if (l < k)
              {
                v = &sol[l * nb];
                M_block = &values[j * b_sqr];
                LU_FIND_ROOT_UPDATE (nb, M_block, v, r);
              }
          }
        BS_ASSERT (j1 == j2);
        BS_ASSERT (cols[j2] == k) (cols[j2]) (k);
        D_block = &values[j2 * b_sqr];
        uLU_FIND_ROOT_L (nb, D_block, r);
      }
    // solve Ux = y
    r = &sol [(n - 1) * nb];
    for (k = n - 1; k >= 0; --k, r -= nb)
      {
        // find element in k column
        index_t j0  = rows[k];
        j1          = j0;//diag_ind[k];
        j2          = rows[k + 1];
        for (j = j0; j < j2; ++j)
          {
            l = cols[j];
            if (l > k)
              {
                v = &sol[l * nb];
                M_block = &values[j * b_sqr];
                LU_FIND_ROOT_UPDATE (nb, M_block, v, r);
              }
          }
        BS_ASSERT (j0 == j1);
        BS_ASSERT (cols[j1] == k) (cols[j1]) (k);
        D_block = &values[j1 * b_sqr];
        uLU_FIND_ROOT_U (nb, D_block, r);
      }

    return 0;
  }


  /**
   * \brief setup preconditioner (no matrix merging)
   *
   * \param matrix The BCSR matrix
   * \return 0 is success
   */
  int
  csr_ilu_prec::setup_internal (bcsr_matrix_t *matrix)
  {
    BS_ASSERT (matrix);

    typedef item_t fp_type;
    typedef index_t i_type;
    typedef typename strategy_t::rhs_item_array_t rhs_item_array_t;
    typedef typename strategy_t::rhs_item_t rhs_item_t;

#ifdef _DEBUG
    //matrix->ascii_write_in_csr_format ("my_csr_ilu_setup.csr");
#endif

    BS_ASSERT (sp_ilu);
    bcsr_matrix_t &ilu = *sp_ilu;

    // common
    i_type i, j, j1, j2, j3, cl, k;

    rhs_item_t *block;


    if (ilu.init (*matrix))
      {
        bs_throw_exception ("CSR_ILU: Can't init ilu");
      }

    if (!check_state_internal ())
      {
        bs_throw_exception ("CSR_ILU: Internal state of CSR_ILU_PREC is invalid");
      }

    BS_ASSERT (ilu.n_cols == matrix->n_cols) (ilu.n_cols) (matrix->n_cols);

    i_type n              = ilu.n_rows;
    i_type nb             = ilu.n_block_size = matrix->n_block_size;
    i_type b_sqr          = nb * nb;
    ilu.n_cols            = matrix->n_cols;

    const index_array_t &ilu_rows     = ilu.get_rows_ptr ();
    const index_array_t &ilu_cols     = ilu.get_cols_ind ();
    rhs_item_array_t  &ilu_values         = ilu.get_values ();

    //ilu.ascii_write_in_csr_format ("ilu_set_internal");

    //  ilu.write_matrix_to_file ("matrix.out");
    // ----------------------
    // STAGE 3: build ilu factorization
    int i_str;
    //fp_type *d_block, *dd_block;
    rhs_item_t *d_block, *dd_block;
    fp_type d;
    int jj, jj1, jj2;

    // loop through strings
    for (i = 0; i < n; ++i)
      {
        // update i-th string by 0..i-1 strings
        j1 = ilu_rows[i];
        j2 = j1;//ilu_diag_ind[i];
        j3 = ilu_rows[i + 1];

        BS_ASSERT (j1 == j2) (j1) (j2);

        // for all nonzero elements in i-th string (in L part)
        // update it by diagonal element and update all other nonzero elements
        // in i-th string
        for (j = j1; j < j3; ++j)
          {
            // find corresponding diagonal element
            i_str = ilu_cols[j];
            if (i_str < i)
              {
                // update by diagonal element
                block   = &ilu_values[j * b_sqr];
                d_block = &ilu_values[ilu_rows[i_str]/*ilu_diag_ind[i_str]*/ * b_sqr];
                //BS_ASSERT (ilu_diag_ind[i_str] == ilu_rows[i_str]) (ilu_diag_ind[i_str]) (ilu_rows[i_str]);
                uLU_SEEK_L (nb, d_block, block, d);

                jj1 = ilu_rows[i_str];
                jj2 = ilu_rows[i_str + 1];
                k = j + 1;

                for (jj = jj1; jj < jj2; ++jj)
                  {
                    cl = ilu_cols[jj];
                    if (cl <= i_str)
                      continue;

                    for (k = j1; k < j3 && ilu_cols[k] < cl; ++k)
                      ;

                    if (k < j3 && ilu_cols[k] == cl)
                      {
                        // upgrade by corresponding element
                        d_block   = &ilu_values[jj * b_sqr];
                        dd_block  = &ilu_values[k * b_sqr];
                        LU_UPGRADE(nb, block, d_block, dd_block);
                      }
                  }
              }
          }

        // factorize i-th string
        block = &ilu_values[j2 * b_sqr];
        uLU (nb, block, d);

        for (j = j1; j < j3; ++j)
          {
            index_t cl = ilu_cols[j];
            if (cl > i)
              {
                d_block = &ilu_values[j * b_sqr];
                uLU_SEEK_U (nb, block, d_block);
              }
          }
      }
    // ----------------------
    return 0;
  }

  /**
   * \brief setup preconditioner (merge matrices if needed)
   *
   * \param matrix Various matrix
   * \return 0 if success
   */
  int
  csr_ilu_prec::setup (matrix_t *matrix_)
  {
    BS_ASSERT (matrix_);
    if (!matrix_)
      {
        bs_throw_exception ("CSR_ILU: Passed matrix is null");
      }

    jacobian_matrix *jmx = dynamic_cast <jacobian_matrix *> (matrix_);
    BS_ASSERT (jmx) (bs::type_name (*matrix_));
    if (jmx)
      {
        const sp_bcsr_matrix_t &prepared_mx = jmx->get_prepared_bcsr_matrix ();
        BS_ASSERT (prepared_mx);

        return setup_internal (prepared_mx);
      }
    else
      {
        bcsr_matrix_t *matrix = dynamic_cast <bcsr_matrix_t *> (matrix_);
        if (!matrix)
          {
            bs_throw_exception ("Passed matrix is not a bcsr matrix");
          }
        return setup_internal (matrix);
      }
  }

  csr_ilu_prec::sp_bcsr_matrix_t
  csr_ilu_prec::get_ilu_matrix () const
  {
    return sp_ilu;
  }

  //////////////////////////////////////////////////////////////////////////
  BLUE_SKY_TYPE_STD_CREATE (csr_ilu_prec);
  BLUE_SKY_TYPE_STD_COPY (csr_ilu_prec);
  BLUE_SKY_TYPE_IMPL (csr_ilu_prec, linear_solver_base, "csr_ilu_prec", "csr_ilu_prec", "csr_ilu_prec");

  //////////////////////////////////////////////////////////////////////////
  //! register types in kernel
  bool csr_ilu_prec_register_type (const blue_sky::plugin_descriptor &pd)
  {
    bool res = true;

    res &= BS_KERNEL.register_type (pd, csr_ilu_prec::bs_type ());

    return res;
  }

} // namespace blue_sky
